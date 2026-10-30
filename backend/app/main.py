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
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware

from . import hybrid, pipeline

DATA_DIR = os.environ.get("LINKWRIGHT_DATA", "/tmp/linkwright-jobs")
FRONTEND = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "frontend"))
if not os.path.isdir(FRONTEND):
    FRONTEND = "/app/frontend"
MAX_JOBS_PER_HOUR = int(os.environ.get("LINKWRIGHT_RATE", "20"))
MAX_UPLOAD_BYTES = int(os.environ.get("LINKWRIGHT_MAX_UPLOAD_BYTES", str(20 * 1024 * 1024)))
RETENTION_DAYS = int(os.environ.get("LINKWRIGHT_RETENTION_DAYS", "7"))
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
app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://singhadd25-sys.github.io"],
    allow_credentials=False,
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Authorization", "Content-Type", "X-API-Key"],
)

# Supabase Auth protects every API endpoint except the lightweight health check.
# The publishable key is safe for client use; never configure a service-role key here.
SUPABASE_URL = os.environ.get("SUPABASE_URL", "").rstrip("/")
SUPABASE_PUBLISHABLE_KEY = os.environ.get("SUPABASE_PUBLISHABLE_KEY", "")

@app.middleware("http")
async def require_supabase_session(request: Request, call_next):
    if request.method != "OPTIONS" and request.url.path.startswith("/api/") and request.url.path != "/api/health":
        if not SUPABASE_URL or not SUPABASE_PUBLISHABLE_KEY:
            return JSONResponse({"detail": "Backend authentication is not configured."}, status_code=503)
        auth = request.headers.get("authorization", "")
        scheme, _, token = auth.partition(" ")
        if scheme.lower() != "bearer" or not token:
            return JSONResponse({"detail": "Sign in to use the Linkwright research workspace."}, status_code=401)
        try:
            import httpx
            async with httpx.AsyncClient(timeout=8.0) as client:
                response = await client.get(
                    f"{SUPABASE_URL}/auth/v1/user",
                    headers={"apikey": SUPABASE_PUBLISHABLE_KEY, "Authorization": f"Bearer {token}"},
                )
            if response.status_code != 200:
                return JSONResponse({"detail": "Your session is invalid or expired. Please sign in again."}, status_code=401)
            request.state.supabase_user = response.json()
        except Exception:
            return JSONResponse({"detail": "Could not verify your session. Please retry."}, status_code=503)
    return await call_next(request)
pool = ThreadPoolExecutor(max_workers=int(os.environ.get("LINKWRIGHT_WORKERS", "2")))
_hits: dict[str, deque] = defaultdict(deque)
_lock = threading.Lock()


def _enforce_rate_limit(request: Request):
    ip = request.client.host if request.client else "unknown"
    now = time.time()
    q = _hits[ip]
    while q and now - q[0] > 3600:
        q.popleft()
    if len(q) >= MAX_JOBS_PER_HOUR:
        raise HTTPException(429, "Hourly job limit reached. Try again later.")
    q.append(now)


def _require_api_key(provided: str | None):
    expected = os.environ.get("LINKWRIGHT_API_KEY", "")
    if expected and not hmac.compare_digest(provided or "", expected):
        raise HTTPException(401, "A valid X-API-Key header is required.")


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


def _require_job_owner(job_id: str, request: Request):
    status = _read_status(job_id)
    user = getattr(request.state, "supabase_user", {}) or {}
    if not user.get("id") or status.get("owner_id") != user.get("id"):
        raise HTTPException(404, "Job not found")
    return status


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


@app.get("/api/health")
def health():
    return {"ok": True, "version": pipeline.VERSION}


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


@app.post("/api/jobs")
async def create_job(
    request: Request,
    fileA: UploadFile | None = File(None),
    fileB: UploadFile | None = File(None),
    smilesA: str = Form(""),
    smilesB: str = Form(""),
    mode: str = Form("both"),
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

    if mode not in ("reactions", "linkers", "both"):
        raise HTTPException(400, "mode must be reactions, linkers or both")
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
    _write_status(job_id, state="running", created=time.time(), owner_id=request.state.supabase_user["id"], stages={})
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
    st = _require_job_owner(job_id, request)
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
    receptor_bytes = await receptor.read()
    if len(receptor_bytes) > 100 * 1024 * 1024 or not receptor_bytes:
        raise HTTPException(413, "Receptor file is empty or exceeds 100 MB.")
    receptor_path = os.path.join(job_dir, "receptor.pdbqt")
    with open(receptor_path, "wb") as f:
        f.write(receptor_bytes)
    ligand_paths = []
    for i, ligand in enumerate(ligands, 1):
        if not ligand.filename or not ligand.filename.lower().endswith(".pdbqt"):
            raise HTTPException(400, "All ligand uploads must be .pdbqt files.")
        data = await ligand.read()
        if not data or len(data) > 20 * 1024 * 1024:
            raise HTTPException(413, f"Ligand {ligand.filename} is empty or exceeds 20 MB.")
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", os.path.basename(ligand.filename))[:80] or f"ligand_{i}.pdbqt"
        path = os.path.join(job_dir, f"ligand_{i}_{safe}")
        with open(path, "wb") as f:
            f.write(data)
        ligand_paths.append(path)
    params = dict(center=center, size=size, exhaustiveness=exhaustiveness, num_modes=num_modes,
                  seed=seed, timeout=timeout)
    _write_status(job_id, state="running", created=time.time(), owner_id=request.state.supabase_user["id"], stages={
        "docking": {"state": "queued", "message": f"{len(ligand_paths)} ligand(s) queued"}})
    pool.submit(_execute_docking, job_id, receptor_path, ligand_paths, params)
    return {"job_id": job_id, "status_url": f"/api/jobs/{job_id}",
            "results_url": f"/api/jobs/{job_id}/files/docking_results.json",
            "package_url": f"/api/jobs/{job_id}/files/docking_package.zip"}


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str, request: Request):
    if not UUID_RE.match(job_id):
        raise HTTPException(404, "Job not found")
    st = _require_job_owner(job_id, request)
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
    _require_job_owner(job_id, request)
    path = os.path.join(DATA_DIR, job_id, name)
    if not os.path.exists(path):
        raise HTTPException(404, "File not ready")
    return FileResponse(path, media_type=ALLOWED_FILES[name], filename=f"linkwright_{job_id[:8]}_{name}")


@app.get("/docking")
def docking_page():
    return FileResponse(os.path.join(FRONTEND, "docking.html"))


@app.get("/structure-tools")
def structure_tools_page():
    return FileResponse(os.path.join(FRONTEND, "structure.html"))


@app.get("/")
def index():
    return FileResponse(os.path.join(FRONTEND, "index.html"))
