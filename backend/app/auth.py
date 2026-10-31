"""SQLite accounts, email verification, seven-day trials and feedback."""
from __future__ import annotations
import hashlib, hmac, os, re, secrets, sqlite3, time
from pathlib import Path

EMAIL_RE = re.compile(r"^[^\s@]+@[^\s@]+\.[^\s@]+$")
ITERATIONS = 310_000
SESSION_SECONDS = int(os.environ.get("LINKWRIGHT_SESSION_DAYS", "7")) * 86400
TRIAL_DAYS = int(os.environ.get("LINKWRIGHT_TRIAL_DAYS", "7"))

def db_path() -> str:
    root = Path(os.environ.get("LINKWRIGHT_DATA", "/tmp/linkwright-jobs")); root.mkdir(parents=True, exist_ok=True)
    return str(root / "accounts.sqlite3")
DATABASE_URL=os.environ.get("LINKWRIGHT_DATABASE_URL","").strip()
IntegrityErrors=(sqlite3.IntegrityError,)
if DATABASE_URL:
    import psycopg
    IntegrityErrors=(sqlite3.IntegrityError,psycopg.errors.IntegrityError)

class _Row(dict):
    """Row addressable by column name or position, like sqlite3.Row."""
    def __getitem__(self,k): return list(self.values())[k] if isinstance(k,int) else dict.__getitem__(self,k)

def _pg_sql(sql):
    s=sql.strip()
    if s.upper().startswith(("PRAGMA","BEGIN")): return None
    return re.sub(r"\bREAL\b","DOUBLE PRECISION",s.replace("?","%s"))

class _PG:
    """Thin wrapper so the SQLite-style code below runs unchanged on Postgres."""
    def __init__(self):
        from psycopg.rows import dict_row
        self.c=psycopg.connect(DATABASE_URL,row_factory=lambda cur:(lambda vals:_Row(zip([d.name for d in cur.description],vals))),prepare_threshold=None,connect_timeout=15)
    def execute(self,sql,params=()):
        q=_pg_sql(sql)
        if q is None: return self.c.cursor()
        return self.c.execute(q,params)
    def __enter__(self): return self
    def __exit__(self,et,ev,tb):
        try:
            (self.c.rollback() if et else self.c.commit())
        finally: self.c.close()

def connect():
    if DATABASE_URL: return _PG()
    c=sqlite3.connect(db_path(), timeout=20); c.row_factory=sqlite3.Row; c.execute("PRAGMA foreign_keys=ON"); return c

def _columns(c,table):
    if DATABASE_URL:
        return {r["column_name"] for r in c.execute("SELECT column_name FROM information_schema.columns WHERE table_schema=current_schema() AND table_name=?",(table,)).fetchall()}
    return {r[1] for r in c.execute(f"PRAGMA table_info({table})")}
def init_db():
    with connect() as c:
        c.execute("""CREATE TABLE IF NOT EXISTS users(id TEXT PRIMARY KEY,email TEXT UNIQUE NOT NULL,password_hash TEXT NOT NULL,salt TEXT NOT NULL,plan TEXT NOT NULL DEFAULT 'pilot',created REAL NOT NULL,role TEXT NOT NULL DEFAULT 'researcher',approved INTEGER NOT NULL DEFAULT 1)""")
        cols=_columns(c,"users")
        migrations={"role":"TEXT NOT NULL DEFAULT 'researcher'","approved":"INTEGER NOT NULL DEFAULT 1","email_verified":"INTEGER NOT NULL DEFAULT 1","trial_started":"REAL","trial_expires":"REAL","last_login":"REAL"}
        for k,v in migrations.items():
            if k not in cols: c.execute(f"ALTER TABLE users ADD COLUMN {k} {v}")
        c.execute("""CREATE TABLE IF NOT EXISTS sessions(token_hash TEXT PRIMARY KEY,user_id TEXT NOT NULL REFERENCES users(id) ON DELETE CASCADE,expires REAL NOT NULL,created REAL NOT NULL)""")
        c.execute("CREATE INDEX IF NOT EXISTS idx_sessions_user ON sessions(user_id)")
        c.execute("CREATE TABLE IF NOT EXISTS app_settings(key TEXT PRIMARY KEY,value TEXT NOT NULL)")
        c.execute("""CREATE TABLE IF NOT EXISTS email_verifications(token_hash TEXT PRIMARY KEY,user_id TEXT NOT NULL REFERENCES users(id) ON DELETE CASCADE,expires REAL NOT NULL,created REAL NOT NULL)""")
        c.execute("""CREATE TABLE IF NOT EXISTS email_otps(user_id TEXT PRIMARY KEY REFERENCES users(id) ON DELETE CASCADE,otp_hash TEXT NOT NULL,expires REAL NOT NULL,attempts INTEGER NOT NULL DEFAULT 0,created REAL NOT NULL)""")
        c.execute("""CREATE TABLE IF NOT EXISTS feedback(id TEXT PRIMARY KEY,user_id TEXT NOT NULL REFERENCES users(id),rating INTEGER NOT NULL,ease TEXT NOT NULL,feature TEXT NOT NULL,comments TEXT NOT NULL DEFAULT '',created REAL NOT NULL)""")
        admin_email=os.environ.get("LINKWRIGHT_ADMIN_EMAIL","").strip().lower()
        if admin_email and EMAIL_RE.match(admin_email):
            c.execute("UPDATE users SET role='admin',approved=1,email_verified=1 WHERE email=?",(admin_email,))
            c.execute("UPDATE users SET trial_started=COALESCE(trial_started,created),trial_expires=NULL WHERE email=?",(admin_email,))
            if c.execute("SELECT 1 FROM users WHERE role='admin' LIMIT 1").fetchone(): c.execute("INSERT INTO app_settings VALUES('admin_bootstrap_used','1') ON CONFLICT(key) DO UPDATE SET value='1'")
        if os.environ.get("LINKWRIGHT_REQUIRE_APPROVAL","0")=="1":
            if not admin_email or not EMAIL_RE.match(admin_email): raise RuntimeError("A valid LINKWRIGHT_ADMIN_EMAIL is required.")
            if not c.execute("SELECT 1 FROM users WHERE role='admin' LIMIT 1").fetchone() and not os.environ.get("LINKWRIGHT_ADMIN_BOOTSTRAP_TOKEN"): raise RuntimeError("Set LINKWRIGHT_ADMIN_BOOTSTRAP_TOKEN before administrator registration.")

def _hash(p,s): return hashlib.pbkdf2_hmac('sha256',p.encode(),s,ITERATIONS,dklen=32)
def _public_user(r):
    now=time.time(); expiry=r['trial_expires']
    return {"id":r['id'],"email":r['email'],"plan":r['plan'],"role":r['role'],"approved":bool(r['approved']),"email_verified":bool(r['email_verified']),"trial_started":r['trial_started'],"trial_expires":expiry,"trial_expired":bool(expiry and expiry<=now)}
def create_user(email,password,admin_token=""):
    email=(email or '').strip().lower()
    if len(email)>254 or not EMAIL_RE.match(email): raise ValueError("Enter a valid email address.")
    if len(password or '')<12 or len(password)>256: raise ValueError("Password must be 12–256 characters long.")
    admin_email=os.environ.get('LINKWRIGHT_ADMIN_EMAIL','').strip().lower(); is_admin=bool(admin_email and email==admin_email)
    uid,salt=secrets.token_hex(16),secrets.token_bytes(16); now=time.time(); digest=_hash(password,salt)
    try:
        with connect() as c:
            c.execute('BEGIN IMMEDIATE')
            if is_admin:
                expected=os.environ.get('LINKWRIGHT_ADMIN_BOOTSTRAP_TOKEN',''); used=c.execute("SELECT value FROM app_settings WHERE key='admin_bootstrap_used'").fetchone()
                if not expected or used or not hmac.compare_digest(admin_token or '',expected): raise ValueError("Administrator setup requires the configured one-time bootstrap token.")
            # Admin bootstrap is trusted and immediate. Researcher accounts are inactive until email verification.
            prior=c.execute("SELECT id,role,email_verified FROM users WHERE email=?",(email,)).fetchone()
            if prior and not prior["email_verified"] and prior["role"]!="admin" and not is_admin:
                # Unverified address: replace the pending password so the real owner can finish signup.
                c.execute("UPDATE users SET password_hash=?,salt=?,created=? WHERE id=?",(digest.hex(),salt.hex(),now,prior["id"]))
                return {"id":prior["id"],"email":email,"plan":"pilot","role":"researcher","approved":False,"email_verified":False,"trial_started":None,"trial_expires":None}
            c.execute("INSERT INTO users(id,email,password_hash,salt,created,role,approved,email_verified,trial_started,trial_expires) VALUES(?,?,?,?,?,?,?,?,?,?)",(uid,email,digest.hex(),salt.hex(),now,'admin' if is_admin else 'researcher',int(is_admin),int(is_admin),now if is_admin else None,None))
            if is_admin: c.execute("INSERT INTO app_settings VALUES('admin_bootstrap_used','1') ON CONFLICT(key) DO UPDATE SET value='1'")
    except IntegrityErrors: raise ValueError("An account with this email already exists.")
    return {"id":uid,"email":email,"plan":"pilot","role":"admin" if is_admin else "researcher","approved":is_admin,"email_verified":is_admin,"trial_started":now if is_admin else None,"trial_expires":None}
def create_verification(user_id):
    token=secrets.token_urlsafe(36); digest=hashlib.sha256(token.encode()).hexdigest(); now=time.time()
    with connect() as c:
        c.execute("DELETE FROM email_verifications WHERE user_id=?",(user_id,)); c.execute("INSERT INTO email_verifications VALUES(?,?,?,?)",(digest,user_id,now+86400,now))
    return token
def verify_email(token):
    digest=hashlib.sha256((token or '').encode()).hexdigest(); now=time.time()
    with connect() as c:
        row=c.execute("SELECT user_id FROM email_verifications WHERE token_hash=? AND expires>?",(digest,now)).fetchone()
        if not row: return None
        c.execute("UPDATE users SET email_verified=1,approved=1,trial_started=?,trial_expires=? WHERE id=?",(now,now+TRIAL_DAYS*86400,row['user_id']))
        c.execute("DELETE FROM email_verifications WHERE user_id=?",(row['user_id'],))
        u=c.execute("SELECT * FROM users WHERE id=?",(row['user_id'],)).fetchone()
    return _public_user(u)
def authenticate(email,password):
    with connect() as c: row=c.execute("SELECT * FROM users WHERE email=?",((email or '').strip().lower(),)).fetchone()
    if not row:
        _hash(password or '',b'linkwright-salt'); return None
    try: ok=hmac.compare_digest(bytes.fromhex(row['password_hash']),_hash(password or '',bytes.fromhex(row['salt'])))
    except (ValueError,TypeError): return None
    if not ok: return None
    if row['role']!='admin' and not row['email_verified']: return {**_public_user(row),"needs_verification":True}
    if row['trial_expires'] and row['trial_expires']<=time.time(): return {**_public_user(row),"trial_expired":True}
    with connect() as c: c.execute("UPDATE users SET last_login=? WHERE id=?",(time.time(),row['id']))
    return _public_user(row)
def new_session(user_id):
    token=secrets.token_urlsafe(32); digest=hashlib.sha256(token.encode()).hexdigest(); now=time.time()
    with connect() as c: c.execute("DELETE FROM sessions WHERE expires<?",(now,)); c.execute("INSERT INTO sessions VALUES(?,?,?,?)",(digest,user_id,now+SESSION_SECONDS,now))
    return token
def get_session(token):
    if not token:return None
    digest=hashlib.sha256(token.encode()).hexdigest()
    with connect() as c: row=c.execute("SELECT users.* FROM users JOIN sessions ON users.id=sessions.user_id WHERE sessions.token_hash=? AND sessions.expires>?",(digest,time.time())).fetchone()
    if not row:return None
    u=_public_user(row)
    if u['trial_expired'] and u['role']!='admin': return None
    return u
def revoke_session(token):
    if token:
        with connect() as c:c.execute("DELETE FROM sessions WHERE token_hash=?",(hashlib.sha256(token.encode()).hexdigest(),))
def list_users():
    with connect() as c: rows=c.execute("SELECT * FROM users ORDER BY created DESC").fetchall()
    return [{**_public_user(r),"created":r['created'],"last_login":r['last_login']} for r in rows]
def set_approved(uid,approved):
    with connect() as c:
        cur=c.execute("UPDATE users SET approved=? WHERE id=? AND role!='admin' AND email_verified=1",(int(approved),uid))
        if cur.rowcount:c.execute("DELETE FROM sessions WHERE user_id=?",(uid,))
    return cur.rowcount>0
def set_plan(uid,plan):
    if plan not in {'pilot','researcher_plus','lab'}:raise ValueError('Invalid plan.')
    with connect() as c:cur=c.execute("UPDATE users SET plan=? WHERE id=? AND role!='admin'",(plan,uid))
    return cur.rowcount>0
def add_feedback(uid,rating,ease,feature,comments):
    fid=secrets.token_hex(16)
    with connect() as c:c.execute("INSERT INTO feedback VALUES(?,?,?,?,?,?,?)",(fid,uid,int(rating),ease,feature,(comments or '')[:4000],time.time()))
    return fid
def feedback_stats():
    with connect() as c:
        total=c.execute("SELECT COUNT(*) n,AVG(rating) avg FROM feedback").fetchone(); rows=c.execute("SELECT f.*,u.email FROM feedback f JOIN users u ON u.id=f.user_id ORDER BY f.created DESC LIMIT 100").fetchall()
        users=c.execute("SELECT COUNT(*) n FROM users WHERE role='researcher'").fetchone()[0]; active=c.execute("SELECT COUNT(*) n FROM users WHERE role='researcher' AND email_verified=1 AND approved=1 AND (trial_expires IS NULL OR trial_expires>?)",(time.time(),)).fetchone()[0]
    return {"registered_researchers":users,"active_trials":active,"feedback_count":total['n'],"average_rating":round(total['avg'],2) if total['avg'] is not None else None,"feedback":[dict(r) for r in rows]}


def create_email_otp(user_id):
    """Create a six-digit, single-use email OTP valid for ten minutes."""
    code=f"{secrets.randbelow(1_000_000):06d}"
    digest=hashlib.sha256(code.encode("ascii")).hexdigest()
    now=time.time()
    with connect() as c:
        c.execute("INSERT INTO email_otps(user_id,otp_hash,expires,attempts,created) VALUES(?,?,?,?,?) "
                  "ON CONFLICT(user_id) DO UPDATE SET otp_hash=excluded.otp_hash,expires=excluded.expires,attempts=0,created=excluded.created",
                  (user_id,digest,now+600,0,now))
    return code

def verify_email_otp(email, code):
    """Verify OTP with expiry and five-attempt limit; activate the seven-day trial."""
    email=(email or "").strip().lower()
    code=(code or "").strip()
    if not re.fullmatch(r"\d{6}",code):
        return None
    now=time.time()
    with connect() as c:
        row=c.execute("SELECT u.id,o.otp_hash,o.expires,o.attempts FROM users u JOIN email_otps o ON o.user_id=u.id WHERE u.email=? AND u.email_verified=0",(email,)).fetchone()
        if not row or row["expires"]<=now or row["attempts"]>=5:
            return None
        digest=hashlib.sha256(code.encode("ascii")).hexdigest()
        if not hmac.compare_digest(row["otp_hash"],digest):
            c.execute("UPDATE email_otps SET attempts=attempts+1 WHERE user_id=?",(row["id"],))
            return None
        c.execute("UPDATE users SET email_verified=1,approved=1,trial_started=?,trial_expires=? WHERE id=?",(now,now+TRIAL_DAYS*86400,row["id"]))
        c.execute("DELETE FROM email_otps WHERE user_id=?",(row["id"],))
        c.execute("DELETE FROM email_verifications WHERE user_id=?",(row["id"],))
        user=c.execute("SELECT * FROM users WHERE id=?",(row["id"],)).fetchone()
    return _public_user(user)

def get_unverified_user(email):
    with connect() as c:
        row=c.execute("SELECT id,email,email_verified FROM users WHERE email=?",(email or "").strip().lower()).fetchone()
    return dict(row) if row else None
