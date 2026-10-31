"""Linkwright web service (FastAPI).

Jobs run in a small thread pool and are stored on disk under DATA_DIR, so a restart does not
lose finished results. For heavier public traffic, move `_execute` into a Celery/RQ worker.
"""
from __future__ import annotations

import json
import os
import hmac
import re
import threading
import time
import uuid
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor

from fastapi import FastAPI, File, Form, Header, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, RedirectResponse, HTMLResponse

from . import hybrid, pipeline

DATA_DIR = os.environ.get("LINKWRIGHT_DATA", "/tmp/linkwright-jobs")
FRONTEND = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "frontend"))
if not os.path.isdir(FRONTEND):
    FRONTEND = "/app/frontend"
MAX_JOBS_PER_HOUR = int(os.environ.get("LINKWRIGHT_RATE", "20"))
MAX_UPLOAD_BYTES = int(os.environ.get("LINKWRIGHT_MAX_UPLOAD_BYTES", str(20 * 1024 * 1024)))
RETENTION_DAYS = int(os.environ.get("LINKWRIGHT_RETENTION_DAYS", "7"))
PLAN_LIMITS = {"pilot": 10, "researcher_plus": 20, "lab": 50}
ALLOWED_FILES = {
    "report.docx": "application/vnd.openxmlformats-officedocument.wordprocessingml.document",
    "report.md": "text/markdown",
    "hybrids.sdf": "chemical/x-mdl-sdfile",
    "gromacs_package.zip": "application/zip",
    "results.json": "application/json",
    "docking_package.zip": "application/zip",
    "docking_results.json": "application/json",
}
UUID_RE = re.compile(r"^[0-9a-f]{32}$")

os.makedirs(DATA_DIR, exist_ok=True)
app = FastAPI(title="Linkwright", version=pipeline.VERSION)
from . import auth
auth.init_db()
SESSION_COOKIE = "linkwright_session"

@app.middleware("http")
async def account_and_ownership_guard(request: Request, call_next):
    path = request.url.path
    public_paths = {"/api/health", "/api/auth/register", "/api/auth/login", "/api/auth/me", "/api/auth/verify", "/api/auth/verify-otp", "/api/auth/resend-otp", "/verify"}
    # Reject cross-origin browser mutations, including auth endpoints. SameSite cookies
    # provide a second barrier; non-browser clients may omit Origin.
    origin = request.headers.get("origin")
    if path.startswith("/api/") and request.method not in ("GET", "HEAD", "OPTIONS") and origin:
        from urllib.parse import urlsplit
        if urlsplit(origin).netloc.lower() != request.headers.get("host", "").lower():
            return JSONResponse({"detail": "Cross-origin request rejected."}, status_code=403)
    if path.startswith("/api/") and path not in public_paths:
        user = auth.get_session(request.cookies.get(SESSION_COOKIE))
        if not user:
            return JSONResponse({"detail": "Please sign in to use Linkwright."}, status_code=401)
        if not user.get("approved"):
            return JSONResponse({"detail": "Your researcher account is not active."}, status_code=403)
        if user.get("trial_expired") and user.get("role") != "admin":
            return JSONResponse({"detail": "Your 7-day trial has expired. Contact Linkwright to extend access."}, status_code=403)
        request.state.user = user
    response = await call_next(request)
    response.headers.setdefault("X-Content-Type-Options", "nosniff")
    response.headers.setdefault("Referrer-Policy", "same-origin")
    response.headers.setdefault("X-Frame-Options", "DENY")
    response.headers.setdefault("Permissions-Policy", "camera=(), microphone=(), geolocation=()")
    return response

pool = ThreadPoolExecutor(max_workers=int(os.environ.get("LINKWRIGHT_WORKERS", "2")))
_hits: dict[str, deque] = defaultdict(deque)
_auth_hits: dict[str, deque] = defaultdict(deque)
_lock = threading.Lock()


def _enforce_rate_limit(request: Request):
    ip = (getattr(request.state, "user", {}) or {}).get("id") or (request.client.host if request.client else "unknown")
    now = time.time()
    q = _hits[ip]
    while q and now - q[0] > 3600:
        q.popleft()
    user = getattr(request.state, "user", None)
    limit = PLAN_LIMITS.get((user or {}).get("plan", "pilot"), MAX_JOBS_PER_HOUR)
    if len(q) >= limit:
        raise HTTPException(429, f"Hourly job limit reached for the {(user or {}).get('plan', 'pilot')} plan ({limit}/hour). Try again later.")
    q.append(now)


def _require_api_key(provided: str | None):
    expected = os.environ.get("LINKWRIGHT_API_KEY", "")
    # Researcher endpoints are protected by the account session middleware.
    # Preserve optional key validation for operators who explicitly send a key,
    # but never require a shared browser-side secret for ordinary researcher use.
    if expected and provided and not hmac.compare_digest(provided, expected):
        raise HTTPException(401, "Invalid X-API-Key header.")


def _status_path(job_id):
    return os.path.join(DATA_DIR, job_id, "status.json")


def _write_status(job_id, **kw):
    with _lock:
        path = _status_path(job_id)
        cur = {}
        if os.path.exists(path):
            with open(path) as f:
                cur = json.load(f)
        cur.update(kw)
        with open(path, "w") as f:
            json.dump(cur, f)


def _execute(job_id, inputs, params):
    job_dir = os.path.join(DATA_DIR, job_id)

    def progress(stage, state, message):
        stages = _read_status(job_id).get("stages", {})
        stages[stage] = dict(state=state, message=message)
        _write_status(job_id, stages=stages)

    try:
        res = pipeline.run(job_dir, inputs, params, progress)
        _write_status(job_id, state="done", warnings=res.get("warnings", []), finished=time.time())
    except Exception as e:  # noqa: BLE001
        _write_status(job_id, state="failed", error=str(e)[:500], finished=time.time())



def _execute_docking(job_id, receptor_path, ligand_paths, params):
    """Run a Vina docking batch in the background and package poses + provenance."""
    import zipfile
    from .vina import run_vina
    job_dir = os.path.join(DATA_DIR, job_id)
    results = {"job_id": job_id, "engine": "AutoDock Vina", "parameters": params, "ligands": [], "warnings": []}
    package_path = os.path.join(job_dir, "docking_package.zip")
    try:
        with zipfile.ZipFile(package_path, "w", zipfile.ZIP_DEFLATED) as z:
            for i, ligand_path in enumerate(ligand_paths, 1):
                ligand_name = os.path.basename(ligand_path)
                stem = re.sub(r"[^A-Za-z0-9_-]+", "_", os.path.splitext(ligand_name)[0])[:48] or f"ligand_{i}"
                out_path = os.path.join(job_dir, f"{stem}_vina_out.pdbqt")
                _write_status(job_id, stages={"docking": {"state": "running", "message": f"Docking {i}/{len(ligand_paths)}: {stem}"}})
                try:
                    result = run_vina(
                        receptor_path, ligand_path, out_path,
                        center=tuple(params["center"]), size=tuple(params["size"]),
                        exhaustiveness=params["exhaustiveness"], num_modes=params["num_modes"],
                        seed=params["seed"], timeout=params["timeout"],
                    )
                    result["ligand"] = ligand_name
                    results["ligands"].append(result)
                    z.write(out_path, f"poses/{stem}_vina_out.pdbqt")
                except Exception as e:
                    results["warnings"].append(f"{ligand_name}: {type(e).__name__}: {e}")
                    results["ligands"].append({"ligand": ligand_name, "status": "failed", "error": str(e)[:1000]})
            z.writestr("README.txt",
                "AutoDock Vina docking batch. Receptor and ligand preparation are researcher responsibilities.\\n"
                "Scores are approximate ranking estimates, not measured binding affinities. Validate poses.\\n")
        results["completed"] = sum(1 for x in results["ligands"] if x.get("scores") is not None)
        with open(os.path.join(job_dir, "docking_results.json"), "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        _write_status(job_id, state="done", warnings=results["warnings"], finished=time.time(),
                      stages={"docking": {"state": "done", "message": "Docking batch finished"}})
    except Exception as e:
        _write_status(job_id, state="failed", error=f"{type(e).__name__}: {e}"[:500], finished=time.time())

def _read_status(job_id):
    path = _status_path(job_id)
    if not os.path.exists(path):
        raise HTTPException(404, "Job not found")
    with open(path) as f:
        return json.load(f)


def _cleanup():
    cutoff = time.time() - RETENTION_DAYS * 86400
    for d in os.listdir(DATA_DIR):
        p = os.path.join(DATA_DIR, d)
        if os.path.isdir(p) and os.path.getmtime(p) < cutoff:
            for root, _, files in os.walk(p, topdown=False):
                for fn in files:
                    os.remove(os.path.join(root, fn))
                os.rmdir(root)


def _parse_sites(s):
    if not s or not s.strip():
        return None
    try:
        return [int(x) for x in re.split(r"[,\s]+", s.strip()) if x != ""]
    except ValueError:
        raise HTTPException(400, "Atom indices must be whole numbers separated by commas.")


def _client_ip(request: Request) -> str:
    """Real client IP. Behind Render's proxy request.client is the proxy, so every
    visitor would share one rate-limit bucket. Enable only on a trusted proxy."""
    if os.environ.get("LINKWRIGHT_TRUST_PROXY", "0") == "1":
        fwd = request.headers.get("x-forwarded-for", "")
        if fwd:
            return fwd.split(",")[0].strip()[:64]
    return request.client.host if request.client else "unknown"

def _limit_email(email: str, action: str, maximum: int):
    """Per-address cap so nobody can use signup/resend to flood a mailbox."""
    key = f"{action}:{(email or '').strip().lower()}"
    now = time.time(); q = _auth_hits[key]
    while q and now - q[0] > 3600: q.popleft()
    if len(q) >= maximum: raise HTTPException(429, "Too many codes requested for this address. Try again in an hour.")
    q.append(now)

def _limit_auth(request: Request, action: str, maximum: int = 10):
    ip = _client_ip(request)
    key = f"{action}:{ip}"
    now = time.time()
    q = _auth_hits[key]
    while q and now - q[0] > 3600:
        q.popleft()
    if len(q) >= maximum:
        raise HTTPException(429, "Too many account attempts. Try again later.")
    q.append(now)

def _set_session_cookie(response, token: str):
    response.set_cookie(SESSION_COOKIE, token, httponly=True,
        secure=os.environ.get("LINKWRIGHT_COOKIE_SECURE", "0") == "1",
        samesite="strict", max_age=auth.SESSION_SECONDS, path="/")

@app.get("/api/health")
def health():
    return {"ok": True, "version": pipeline.VERSION}


def _send_otp_email(email: str, code: str):
    """Deliver a short-lived six-digit OTP. Console mode is for local testing only."""
    import smtplib
    from email.message import EmailMessage
    text=f"Your Linkwright verification code is {code}. It expires in 10 minutes. If you did not request this code, ignore this email."
    # Prefer HTTPS email APIs: Render Free blocks outbound SMTP ports.
    import json, urllib.request
    sender=os.environ.get("LINKWRIGHT_MAIL_FROM", "")
    resend_key=os.environ.get("RESEND_API_KEY", "")
    if resend_key:
        if not sender: raise RuntimeError("LINKWRIGHT_MAIL_FROM is not set to a verified sender.")
        req=urllib.request.Request("https://api.resend.com/emails", method="POST",
            data=json.dumps({"from":f"Linkwright <{sender}>","to":[email],
                             "subject":"Your Linkwright verification code","text":text}).encode(),
            headers={"Authorization":f"Bearer {resend_key}","content-type":"application/json","accept":"application/json"})
        try:
            with urllib.request.urlopen(req, timeout=20) as response:
                if response.status >= 300: raise RuntimeError(f"Resend API returned HTTP {response.status}")
        except urllib.error.HTTPError as exc:
            detail=exc.read(1000).decode("utf-8","replace")
            raise RuntimeError(f"Resend API returned HTTP {exc.code}: {detail}") from exc
        return None
    brevo_key=os.environ.get("LINKWRIGHT_BREVO_API_KEY", "")
    if brevo_key:
        if not sender: raise RuntimeError("LINKWRIGHT_MAIL_FROM is not set to a verified sender.")
        req=urllib.request.Request("https://api.brevo.com/v3/smtp/email", method="POST",
            data=json.dumps({"sender":{"name":"Linkwright","email":sender},"to":[{"email":email}],
                             "subject":"Your Linkwright verification code","textContent":text}).encode(),
            headers={"api-key":brevo_key,"content-type":"application/json","accept":"application/json"})
        with urllib.request.urlopen(req, timeout=20) as response:
            if response.status >= 300: raise RuntimeError(f"Brevo API returned HTTP {response.status}")
        return None
    host=os.environ.get("LINKWRIGHT_SMTP_HOST", "")
    if not host:
        if os.environ.get("LINKWRIGHT_EMAIL_MODE", "").lower()=="console":
            print(f"[LOCAL LINKWRIGHT OTP] recipient={email} code={code} expires_in=600s", flush=True)
            return code
        raise RuntimeError("Email delivery is not configured. Configure SMTP before accepting public registrations.")
    msg=EmailMessage()
    msg["Subject"]="Your Linkwright verification code"
    msg["From"]=os.environ.get("LINKWRIGHT_SMTP_FROM", os.environ.get("LINKWRIGHT_SMTP_USER", ""))
    msg["To"]=email
    msg.set_content(f"Your Linkwright verification code is {code}. It expires in 10 minutes. If you did not request this code, ignore this email.")
    port=int(os.environ.get("LINKWRIGHT_SMTP_PORT", "587"))
    user=os.environ.get("LINKWRIGHT_SMTP_USER", "")
    password=os.environ.get("LINKWRIGHT_SMTP_PASSWORD", "")
    with smtplib.SMTP(host,port,timeout=20) as server:
        if os.environ.get("LINKWRIGHT_SMTP_STARTTLS", "1")=="1": server.starttls()
        if user: server.login(user,password)
        server.send_message(msg)
    return None

@app.post("/api/auth/register")
async def register_account(request: Request):
    _limit_auth(request, "register", 5)
    try: payload=await request.json()
    except Exception: raise HTTPException(400,"Request body must be valid JSON.")
    if not isinstance(payload,dict): raise HTTPException(400,"Request body must be a JSON object.")
    try:
        _limit_email(str(payload.get("email","")), "otp_send", 5)
        user=auth.create_user(str(payload.get("email","")),str(payload.get("password","")),str(payload.get("admin_token","")))
        if user.get("role")=="admin":
            token=auth.new_session(user["id"])
            response=JSONResponse({"user":user,"message":"Administrator account created."})
            _set_session_cookie(response,token)
            return response
        code=auth.create_email_otp(user["id"])
        dev_code=_send_otp_email(user["email"],code)
        result={"verification_required":True,"message":"We sent a six-digit code to your email. Enter it within 10 minutes to activate your 7-day trial."}
        if dev_code: result["development_otp"]=dev_code
        return JSONResponse(result,status_code=202)
    except ValueError as exc: raise HTTPException(400,str(exc))
    except Exception as exc: raise HTTPException(503,f"Registration could not send the verification code. Check email delivery settings and retry. ({type(exc).__name__})")

@app.post("/api/auth/verify-otp")
async def verify_otp_account(request: Request):
    _limit_auth(request, "verify_otp", 10)
    try: payload=await request.json()
    except Exception: raise HTTPException(400,"Request body must be valid JSON.")
    user=auth.verify_email_otp(str(payload.get("email","")),str(payload.get("otp","")))
    if not user: raise HTTPException(400,"The code is invalid, expired, or has reached its attempt limit. Request a new code.")
    session=auth.new_session(user["id"])
    response=JSONResponse({"verified":True,"user":user,"message":"Email verified. Your 7-day Linkwright trial is active."})
    _set_session_cookie(response,session)
    return response

@app.post("/api/auth/resend-otp")
async def resend_otp_account(request: Request):
    _limit_auth(request, "resend_otp", 3)
    try: payload=await request.json()
    except Exception: raise HTTPException(400,"Request body must be valid JSON.")
    email=str(payload.get("email","")).strip().lower()
    pending=auth.get_unverified_user(email)
    # Use the same outward response for unknown and already-verified addresses.
    if pending and not pending.get("email_verified"):
        _limit_email(email, "otp_send", 5)
        code=auth.create_email_otp(pending["id"])
        dev_code=_send_otp_email(email,code)
        result={"message":"If this address has a pending registration, a new verification code has been sent."}
        if dev_code: result["development_otp"]=dev_code
        return result
    return {"message":"If this address has a pending registration, a new verification code has been sent."}

@app.post("/api/auth/login")
async def login_account(request: Request):
    _limit_auth(request, "login", 20)
    try:
        payload = await request.json()
    except Exception:
        raise HTTPException(400, "Request body must be valid JSON.")
    if not isinstance(payload, dict):
        raise HTTPException(400, "Request body must be a JSON object.")
    user = auth.authenticate(str(payload.get("email", "")), str(payload.get("password", "")))
    if not user:
        raise HTTPException(401, "Email or password is incorrect.")
    if user.get("needs_verification"):
        raise HTTPException(403, "Please verify your email using the six-digit code we sent before signing in.")
    if user.get("trial_expired"):
        raise HTTPException(403, "Your 7-day trial has expired. Contact Linkwright to request an extension.")
    if not user.get("approved"):
        raise HTTPException(403, "Your researcher account is awaiting administrator approval.")
    token = auth.new_session(user["id"])
    response = JSONResponse({"user": user})
    _set_session_cookie(response, token)
    return response

@app.post("/api/auth/logout")
async def logout_account(request: Request):
    auth.revoke_session(request.cookies.get(SESSION_COOKIE))
    response = JSONResponse({"ok": True})
    response.delete_cookie(SESSION_COOKIE, path="/", httponly=True,
                           secure=os.environ.get("LINKWRIGHT_COOKIE_SECURE", "0") == "1", samesite="strict")
    return response

@app.get("/api/auth/me")
def current_account(request: Request):
    user = auth.get_session(request.cookies.get(SESSION_COOKIE))
    return {"authenticated": bool(user and user.get("approved")), "pending": bool(user and not user.get("approved")), "user": user}

@app.get("/api/admin/users")
def admin_users(request: Request):
    if request.state.user.get("role") != "admin":
        raise HTTPException(403, "Administrator access required.")
    return {"users": auth.list_users()}

@app.post("/api/admin/users/{user_id}/approve")
def admin_approve_user(user_id: str, request: Request, approved: bool = True):
    if request.state.user.get("role") != "admin":
        raise HTTPException(403, "Administrator access required.")
    if not re.fullmatch(r"[0-9a-f]{32}", user_id):
        raise HTTPException(404, "User not found.")
    if not auth.set_approved(user_id, approved):
        raise HTTPException(404, "User not found or cannot modify administrator.")
    return {"ok": True, "user_id": user_id, "approved": approved}

@app.post("/api/admin/users/{user_id}/plan")
def admin_set_user_plan(user_id: str, request: Request, plan: str):
    if request.state.user.get("role") != "admin":
        raise HTTPException(403, "Administrator access required.")
    if not re.fullmatch(r"[0-9a-f]{32}", user_id):
        raise HTTPException(404, "User not found.")
    try:
        changed = auth.set_plan(user_id, plan)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    if not changed:
        raise HTTPException(404, "User not found or cannot modify administrator.")
    return {"ok": True, "user_id": user_id, "plan": plan}

@app.get("/api/dashboard")
def dashboard_data(request: Request):
    user = request.state.user
    jobs = []
    for jid in os.listdir(DATA_DIR):
        if not UUID_RE.match(jid):
            continue
        try:
            st = _read_status(jid)
            if st.get("owner_id") == user["id"]:
                jobs.append({"job_id": jid, "state": st.get("state"), "created": st.get("created"), "finished": st.get("finished"), "error": st.get("error")})
        except (HTTPException, OSError, ValueError):
            continue
    jobs.sort(key=lambda x: x.get("created") or 0, reverse=True)
    return {"user": user, "jobs": jobs[:30], "limits": {"jobs_per_hour": PLAN_LIMITS.get(user.get("plan", "pilot"), MAX_JOBS_PER_HOUR), "retention_days": RETENTION_DAYS}, "billing": "Pilot access; payment processing is not configured."}

@app.post("/api/feedback")
async def submit_feedback(request: Request):
    try: p=await request.json()
    except Exception: raise HTTPException(400,"Invalid feedback payload.")
    try: rating=int(p.get("rating",0))
    except Exception: rating=0
    ease=str(p.get("ease","")); feature=str(p.get("feature","")); comments=str(p.get("comments",""))
    if rating not in range(1,6): raise HTTPException(400,"Rating must be from 1 to 5.")
    if ease not in {"Very easy","Easy","Difficult","Very difficult"}: raise HTTPException(400,"Select an ease-of-use response.")
    if feature not in {"Hybrid Molecule Design","Reaction Discovery","Structure Analysis","Molecular Docking","Molecular Dynamics"}: raise HTTPException(400,"Select a feature.")
    fid=auth.add_feedback(request.state.user["id"],rating,ease,feature,comments)
    return {"ok":True,"feedback_id":fid,"message":"Thank you. Your feedback has been recorded."}

@app.get("/api/admin/analytics")
def admin_analytics(request: Request):
    if request.state.user.get("role")!="admin": raise HTTPException(403,"Administrator access required.")
    return auth.feedback_stats()

@app.get("/api/options")
def options():
    return {
        "rules": [dict(id=r["id"], name=r["name"], roles=r["roles"]) for r in hybrid.RULES],
        "linkers": [dict(id=l["id"], name=l["name"]) for l in hybrid.LINKERS],
    }

@app.get("/api/reactions/status")
def reaction_database_status(x_api_key: str | None = Header(None)):
    """Report whether the local RetroRules source and searchable index are available."""
    _require_api_key(x_api_key)
    from .reaction_db import ReactionTemplateStore
    return ReactionTemplateStore().status()


@app.post("/api/reactions/index")
def index_reaction_database(force: bool = False, x_api_key: str | None = Header(None)):
    """Build a persistent SQLite index from the locally installed RetroRules CSV.GZ."""
    _require_api_key(x_api_key)
    from .reaction_db import ReactionTemplateStore
    try:
        return ReactionTemplateStore().build_index(force=force)
    except FileNotFoundError as exc:
        raise HTTPException(404, str(exc))
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, f"Reaction database indexing failed: {type(exc).__name__}: {exc}")


@app.get("/api/reactions/search")
def search_reaction_templates(
    q: str | None = None,
    template_id: str | None = None,
    dataset: str | None = None,
    valid: bool | None = None,
    limit: int = 20,
    offset: int = 0,
    validate_smarts: bool = False,
    x_api_key: str | None = Header(None),
):
    """Search indexed reaction templates; SMARTS execution is intentionally not performed here."""
    _require_api_key(x_api_key)
    from .reaction_db import ReactionTemplateStore
    try:
        return ReactionTemplateStore().search(q=q, template_id=template_id, dataset=dataset,
                                              valid=valid, limit=limit, offset=offset,
                                              validate_smarts=validate_smarts)
    except RuntimeError as exc:
        raise HTTPException(409, str(exc))
    except ValueError as exc:
        raise HTTPException(400, str(exc))


@app.get("/api/reactions/template/{template_id}")
def get_reaction_template(template_id: str, validate_smarts: bool = True,
                          x_api_key: str | None = Header(None)):
    """Retrieve one template and optionally check whether RDKit parses its reaction SMARTS."""
    _require_api_key(x_api_key)
    from .reaction_db import ReactionTemplateStore
    try:
        result = ReactionTemplateStore().get(template_id, validate_smarts=validate_smarts)
    except RuntimeError as exc:
        raise HTTPException(409, str(exc))
    if result is None:
        raise HTTPException(404, "Reaction template not found")
    return result


@app.post("/api/reactions/analyze")
async def analyze_reaction_pair(request: Request, x_api_key: str | None = Header(None)):
    """Detect parent functional groups and evaluate curated two-parent forward reactions."""
    _require_api_key(x_api_key)
    try:
        payload = await request.json()
        from .reaction_engine import analyze_pair
        return analyze_pair(
            str(payload.get("smiles_a", "")), str(payload.get("smiles_b", "")),
            str(payload.get("name_a", "Parent A")), str(payload.get("name_b", "Parent B")),
        )
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, f"Reaction analysis failed: {type(exc).__name__}: {exc}")


@app.post("/api/reactions/generate")
async def generate_forward_reactions(request: Request, x_api_key: str | None = Header(None)):
    """Generate products only with Linkwright's curated forward reaction templates."""
    _require_api_key(x_api_key)
    try:
        payload = await request.json()
        from .reaction_engine import analyze_pair
        result = analyze_pair(
            str(payload.get("smiles_a", "")), str(payload.get("smiles_b", "")),
            str(payload.get("name_a", "Parent A")), str(payload.get("name_b", "Parent B")),
        )
        return {"products": result["compatible_forward_reactions"],
                "reaction_count": result["reaction_count"],
                "parents": result["parents"], "database": result["retrorules_database"],
                "limitations": result["limitations"]}
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, f"Forward reaction generation failed: {type(exc).__name__}: {exc}")


@app.get("/api/reactions/experimental-validation-template")
def experimental_validation_template(x_api_key: str | None = Header(None)):
    """Return a structured bench-validation checklist; this endpoint does not perform experiments."""
    _require_api_key(x_api_key)
    return {
        "status": "protocol_template_only",
        "notice": "Linkwright cannot experimentally synthesize or verify a compound. A qualified laboratory must perform and document the work.",
        "required_records": [
            {"field": "reaction_id", "purpose": "Unique link to the computational prediction"},
            {"field": "starting_material_identity", "purpose": "Supplier, lot, purity and identity confirmation"},
            {"field": "reaction_conditions", "purpose": "Actual equivalents, concentration, solvent, temperature, time, atmosphere and catalyst/reagents"},
            {"field": "workup_and_purification", "purpose": "Quench, extraction, chromatography/crystallization and isolated mass"},
            {"field": "yield_percent", "purpose": "Isolated yield with calculation basis"},
            {"field": "analytical_evidence", "purpose": "Raw and processed 1H/13C NMR, HRMS/LC-MS, HPLC/UPLC and other relevant data"},
            {"field": "identity_and_purity_conclusion", "purpose": "Chemist-reviewed structural assignment and measured purity"},
            {"field": "deviations_and_failures", "purpose": "Failed runs, side products, safety observations and protocol deviations"},
            {"field": "reviewer_and_date", "purpose": "Responsible researcher, review date and institutional record reference"}
        ],
        "minimum_evidence_for_experimental_claim": ["independent analytical evidence supporting identity", "measured purity", "documented isolated yield", "traceable experimental procedure"],
        "safety": "Perform only under institutional risk assessment, approved SOPs, trained supervision and appropriate containment/PPE."
    }


@app.post("/api/jobs")
async def create_job(
    request: Request,
    fileA: UploadFile | None = File(None),
    fileB: UploadFile | None = File(None),
    smilesA: str = Form(""),
    smilesB: str = Form(""),
    mode: str = Form("auto"),
    max_hybrids: int = Form(12),
    xref_top: int = Form(5),
    md_top: int = Form(3),
    md_ns: int = Form(100),
    temp: float = Form(310.15),
    linkers: str = Form(""),
    sitesA: str = Form(""),
    sitesB: str = Form(""),
    x_api_key: str | None = Header(None),
):
    _require_api_key(x_api_key)
    _enforce_rate_limit(request)

    if mode not in ("auto", "reactions", "linkers", "both"):
        raise HTTPException(400, "mode must be auto, reactions, linkers or both")
    if not (1 <= max_hybrids <= 50 and 0 <= xref_top <= 10 and 0 <= md_top <= 5 and 1 <= md_ns <= 1000):
        raise HTTPException(400, "A setting is outside its allowed range.")
    valid_linkers = {l["id"] for l in hybrid.LINKERS}
    chosen = [x for x in re.split(r"[,\s]+", linkers.strip()) if x] or None
    if chosen and not set(chosen) <= valid_linkers:
        raise HTTPException(400, "Unknown linker id.")

    inputs = {}
    for key, up, smi in (("A", fileA, smilesA), ("B", fileB, smilesB)):
        data = await up.read(MAX_UPLOAD_BYTES + 1) if up is not None and up.filename else None
        if data is not None and len(data) > MAX_UPLOAD_BYTES:
            raise HTTPException(413, f"Parent {key} upload exceeds {MAX_UPLOAD_BYTES} bytes.")
        inputs[key] = (up.filename if up is not None and up.filename else None, data, smi)
        if data is None and not smi.strip():
            raise HTTPException(400, f"Parent {key}: upload a file or paste a SMILES string.")

    params = dict(mode=mode, max_hybrids=max_hybrids, xref_top=xref_top, md_top=md_top, md_ns=md_ns,
                  temp=temp, linkers=chosen, sites_a=_parse_sites(sitesA), sites_b=_parse_sites(sitesB))
    job_id = uuid.uuid4().hex
    os.makedirs(os.path.join(DATA_DIR, job_id))
    _write_status(job_id, state="running", created=time.time(), owner_id=request.state.user["id"], stages={})
    pool.submit(_execute, job_id, inputs, params)
    if len(os.listdir(DATA_DIR)) % 25 == 0:
        threading.Thread(target=_cleanup, daemon=True).start()
    return {"job_id": job_id}


@app.post("/api/jobs/{job_id}/admet-csv")
async def import_admet_csv(job_id: str, request: Request, file: UploadFile = File(...), source: str = Form("User-imported ADMET CSV"),
                          x_api_key: str | None = Header(None)):
    _require_api_key(x_api_key)
    """Merge externally generated ADMET endpoints into a completed job report."""
    if not UUID_RE.match(job_id):
        raise HTTPException(404, "Job not found")
    st = _read_status(job_id)
    if st.get("owner_id") != request.state.user["id"]:
        raise HTTPException(404, "Job not found")
    if st.get("state") != "done":
        raise HTTPException(409, "ADMET CSV can only be imported after the job is complete.")
    if not file.filename or not file.filename.lower().endswith(".csv"):
        raise HTTPException(400, "Upload a .csv file.")
    raw = await file.read()
    if len(raw) > 2 * 1024 * 1024:
        raise HTTPException(413, "CSV exceeds the 2 MB limit.")
    from . import admet_csv, report
    job_dir = os.path.join(DATA_DIR, job_id)
    results_path = os.path.join(job_dir, "results.json")
    if not os.path.isfile(results_path):
        raise HTTPException(409, "This job does not contain hybrid-design results.")
    with open(results_path, encoding="utf-8") as f:
        res = json.load(f)
    try:
        imported = admet_csv.import_csv(raw, res.get("hybrids", []), source)
    except ValueError as e:
        raise HTTPException(400, str(e))
    for hb in res.get("hybrids", []):
        if hb["id"] in imported["matched"]:
            hb.setdefault("admet_imports", []).append({
                "source": imported["source"], "queried": report.now(),
                "endpoints": imported["matched"][hb["id"]],
            })
    res.setdefault("admet_import_history", []).append({
        "source": imported["source"], "queried": report.now(),
        "matched_count": imported["matched_count"],
        "unmatched_count": imported["unmatched_count"],
        "endpoints": imported["endpoints"],
    })
    with open(results_path, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(job_dir, "report.md"), "w", encoding="utf-8") as f:
        f.write(report.build_markdown(res))
    report.build_docx(res, os.path.join(job_dir, "report.docx"))
    return imported


@app.post("/api/ligand/prepare")
async def prepare_ligand_endpoint(
    request: Request,
    file: UploadFile | None = File(None),
    smiles: str = Form(""),
    seed: int = Form(20261030),
    x_api_key: str | None = Header(None),
):
    """Generate a reproducible, minimized 3D SDF; this is not PDBQT preparation."""
    _require_api_key(x_api_key)
    _enforce_rate_limit(request)
    if not 0 <= seed <= 2_147_483_647:
        raise HTTPException(400, "Seed must be between 0 and 2147483647.")
    raw = await file.read(5 * 1024 * 1024 + 1) if file is not None and file.filename else None
    if raw is not None and len(raw) > 5 * 1024 * 1024:
        raise HTTPException(413, "Ligand upload exceeds 5 MB.")
    if raw is None and not smiles.strip():
        raise HTTPException(400, "Upload a ligand file or provide SMILES.")
    from . import structure_qc
    try:
        return structure_qc.prepare_ligand(file.filename if file else "ligand.smi", raw, smiles, seed)
    except (ValueError, RuntimeError) as e:
        raise HTTPException(422, str(e))


@app.post("/api/protein/validate")
async def validate_protein_endpoint(
    request: Request,
    file: UploadFile = File(...),
    x_api_key: str | None = Header(None),
):
    """Inventory PDB coordinate records and flag common review requirements."""
    _require_api_key(x_api_key)
    _enforce_rate_limit(request)
    if not file.filename or not file.filename.lower().endswith((".pdb", ".ent")):
        raise HTTPException(400, "Upload a PDB or ENT file.")
    raw = await file.read(100 * 1024 * 1024 + 1)
    if len(raw) > 100 * 1024 * 1024:
        raise HTTPException(413, "Protein PDB exceeds the 100 MB limit.")
    from . import structure_qc
    try:
        return structure_qc.validate_protein_pdb(raw, file.filename)
    except ValueError as e:
        raise HTTPException(422, str(e))


@app.post("/api/docking")
async def create_docking_job(
    request: Request,
    receptor: UploadFile = File(...),
    ligands: list[UploadFile] = File(...),
    center_x: float = Form(...), center_y: float = Form(...), center_z: float = Form(...),
    size_x: float = Form(...), size_y: float = Form(...), size_z: float = Form(...),
    exhaustiveness: int = Form(8), num_modes: int = Form(9), seed: int = Form(2026),
    timeout: int = Form(3600),
    x_api_key: str | None = Header(None),
):
    _require_api_key(x_api_key)
    _enforce_rate_limit(request)
    """Queue a Vina batch using pre-prepared PDBQT files. No chemistry preparation is implied."""
    if not ligands or len(ligands) > 20:
        raise HTTPException(400, "Upload between 1 and 20 ligand PDBQT files.")
    center = [center_x, center_y, center_z]
    size = [size_x, size_y, size_z]
    if any(not (-10000 <= float(x) <= 10000) for x in center):
        raise HTTPException(400, "Docking center coordinate is outside allowed bounds.")
    if any(not (1 <= float(x) <= 100) for x in size):
        raise HTTPException(400, "Each docking box dimension must be between 1 and 100 Å.")
    if not (1 <= exhaustiveness <= 128 and 1 <= num_modes <= 20 and 1 <= timeout <= 86400):
        raise HTTPException(400, "Invalid Vina exhaustiveness, mode count or timeout.")
    job_id = uuid.uuid4().hex
    job_dir = os.path.join(DATA_DIR, job_id)
    os.makedirs(job_dir)
    receptor_bytes = await receptor.read(100 * 1024 * 1024 + 1)
    if len(receptor_bytes) > 100 * 1024 * 1024 or not receptor_bytes:
        raise HTTPException(413, "Receptor file is empty or exceeds 100 MB.")
    receptor_path = os.path.join(job_dir, "receptor.pdbqt")
    with open(receptor_path, "wb") as f:
        f.write(receptor_bytes)
    ligand_paths = []
    for i, ligand in enumerate(ligands, 1):
        if not ligand.filename or not ligand.filename.lower().endswith(".pdbqt"):
            raise HTTPException(400, "All ligand uploads must be .pdbqt files.")
        data = await ligand.read(20 * 1024 * 1024 + 1)
        if not data or len(data) > 20 * 1024 * 1024:
            raise HTTPException(413, f"Ligand {ligand.filename} is empty or exceeds 20 MB.")
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", os.path.basename(ligand.filename))[:80] or f"ligand_{i}.pdbqt"
        path = os.path.join(job_dir, f"ligand_{i}_{safe}")
        with open(path, "wb") as f:
            f.write(data)
        ligand_paths.append(path)
    params = dict(center=center, size=size, exhaustiveness=exhaustiveness, num_modes=num_modes,
                  seed=seed, timeout=timeout)
    _write_status(job_id, state="running", created=time.time(), owner_id=request.state.user["id"], stages={
        "docking": {"state": "queued", "message": f"{len(ligand_paths)} ligand(s) queued"}})
    pool.submit(_execute_docking, job_id, receptor_path, ligand_paths, params)
    return {"job_id": job_id, "status_url": f"/api/jobs/{job_id}",
            "results_url": f"/api/jobs/{job_id}/files/docking_results.json",
            "package_url": f"/api/jobs/{job_id}/files/docking_package.zip"}


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str, request: Request):
    if not UUID_RE.match(job_id):
        raise HTTPException(404, "Job not found")
    st = _read_status(job_id)
    if st.get("owner_id") != request.state.user["id"]:
        raise HTTPException(404, "Job not found")
    if st.get("state") == "done":
        regular = os.path.join(DATA_DIR, job_id, "results.json")
        docking = os.path.join(DATA_DIR, job_id, "docking_results.json")
        if os.path.isfile(regular):
            with open(regular, encoding="utf-8") as f:
                st["result"] = json.load(f)
        elif os.path.isfile(docking):
            with open(docking, encoding="utf-8") as f:
                st["result"] = json.load(f)
    return JSONResponse(st)


@app.get("/api/jobs/{job_id}/files/{name}")
def job_file(job_id: str, name: str, request: Request):
    if not UUID_RE.match(job_id) or name not in ALLOWED_FILES:
        raise HTTPException(404, "File not found")
    st = _read_status(job_id)
    if st.get("owner_id") != request.state.user["id"]:
        raise HTTPException(404, "File not found")
    path = os.path.join(DATA_DIR, job_id, name)
    if not os.path.exists(path):
        raise HTTPException(404, "File not ready")
    return FileResponse(path, media_type=ALLOWED_FILES[name], filename=f"linkwright_{job_id[:8]}_{name}")


def _require_page_user(request: Request):
    user = auth.get_session(request.cookies.get(SESSION_COOKIE))
    if not user or not user.get("approved"):
        return RedirectResponse(url="/", status_code=303)
    return None

@app.get("/docking")
def docking_page(request: Request):
    denied = _require_page_user(request)
    return denied or FileResponse(os.path.join(FRONTEND, "docking.html"))

@app.get("/structure-tools")
def structure_tools_page(request: Request):
    denied = _require_page_user(request)
    return denied or FileResponse(os.path.join(FRONTEND, "structure.html"))

@app.get("/hybrid")
def hybrid_workspace(request: Request):
    denied = _require_page_user(request)
    return denied or FileResponse(os.path.join(FRONTEND, "index.html"))

@app.get("/signup")
def signup_page():
    """Public entry point for self-service researcher registration."""
    return FileResponse(os.path.join(FRONTEND, "dashboard.html"))

@app.get("/")
def index():
    return FileResponse(os.path.join(FRONTEND, "dashboard.html"))
