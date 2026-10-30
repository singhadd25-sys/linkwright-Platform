"""Local, searchable index for RetroRules reaction-template CSV archives.

This module deliberately indexes reaction templates without assuming that every
RetroRules transformation is a forward, two-parent synthesis reaction. Applying
or ranking a transformation is a separate chemistry-validation step.
"""
from __future__ import annotations

import argparse
from contextlib import closing
import csv
import gzip
import hashlib
import os
import sqlite3
import tempfile
import threading
import time
from pathlib import Path
from typing import Any

APP_DIR = Path(__file__).resolve().parent
DEFAULT_DATASET = APP_DIR / "data" / "retrorules" / "templates.csv.gz"
DEFAULT_INDEX = Path(os.environ.get("LINKWRIGHT_DATA", "/tmp/linkwright-jobs")) / "retrorules.sqlite3"
EXPECTED_COLUMNS = {
    "TEMPLATE_ID", "SCORE", "TEMPLATE", "REACTIONS", "RADII", "REACTIONS_COUNT",
    "ECS", "ECS_COUNT", "RADIUS_MIN", "RADIUS_MAX", "VALID", "DATASETS",
}
_DB_LOCK = threading.RLock()


def dataset_path() -> Path:
    return Path(os.environ.get("LINKWRIGHT_RETRORULES_PATH", str(DEFAULT_DATASET)))


def index_path() -> Path:
    return Path(os.environ.get("LINKWRIGHT_RETRORULES_INDEX_PATH", str(DEFAULT_INDEX)))


def _fingerprint(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {"source_size_bytes": st.st_size, "source_mtime_ns": st.st_mtime_ns}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _connect(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(path), timeout=60)
    conn.row_factory = sqlite3.Row
    return conn


class ReactionTemplateStore:
    def __init__(self, source: str | Path | None = None, database: str | Path | None = None):
        self.source = Path(source) if source else dataset_path()
        self.database = Path(database) if database else index_path()

    def _metadata(self) -> dict[str, str]:
        if not self.database.is_file():
            return {}
        try:
            with closing(_connect(self.database)) as conn:
                rows = conn.execute("SELECT key, value FROM metadata").fetchall()
                return {row["key"]: row["value"] for row in rows}
        except sqlite3.Error:
            return {}

    def status(self) -> dict[str, Any]:
        exists = self.source.is_file()
        meta = self._metadata()
        indexed = False
        current_fp = _fingerprint(self.source) if exists else {}
        if meta and exists:
            indexed = (
                meta.get("source_size_bytes") == str(current_fp["source_size_bytes"])
                and meta.get("source_mtime_ns") == str(current_fp["source_mtime_ns"])
            )
        result: dict[str, Any] = {
            "dataset_exists": exists,
            "dataset_path": str(self.source),
            "index_path": str(self.database),
            "indexed": indexed,
            "index_current": indexed,
        }
        if exists:
            result["dataset_size_bytes"] = current_fp["source_size_bytes"]
        if meta:
            for key in ("template_count", "valid_count", "invalid_count", "source_sha256", "indexed_at_utc"):
                if key in meta:
                    result[key] = meta[key]
        if indexed:
            with closing(_connect(self.database)) as conn:
                result["source_distribution"] = [dict(row) for row in conn.execute(
                    "SELECT dataset, template_count, valid_count FROM source_counts ORDER BY template_count DESC"
                )]
        return result

    def build_index(self, force: bool = False, batch_size: int = 5000) -> dict[str, Any]:
        with _DB_LOCK:
            if not self.source.is_file():
                raise FileNotFoundError(f"RetroRules dataset not found: {self.source}")
            if batch_size < 100 or batch_size > 50000:
                raise ValueError("batch_size must be between 100 and 50000")
            current_fp = _fingerprint(self.source)
            existing = self.status()
            if existing.get("indexed") and not force:
                existing["message"] = "Index already matches the source dataset; no rebuild needed."
                return existing

            self.database.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(prefix="retrorules-", suffix=".sqlite3", dir=str(self.database.parent))
            os.close(fd)
            tmp_path = Path(tmp_name)
            total = valid_count = invalid_count = 0
            source_counts: dict[str, list[int]] = {}
            try:
                conn = _connect(tmp_path)
                conn.execute("PRAGMA journal_mode=OFF")
                conn.execute("PRAGMA synchronous=OFF")
                conn.execute("PRAGMA temp_store=MEMORY")
                conn.execute("""CREATE TABLE templates (
                    template_id TEXT PRIMARY KEY,
                    score TEXT, template TEXT NOT NULL, reactions TEXT, radii TEXT,
                    reactions_count INTEGER, ecs TEXT, ecs_count INTEGER,
                    radius_min INTEGER, radius_max INTEGER, valid INTEGER NOT NULL,
                    datasets TEXT NOT NULL
                )""")
                conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
                conn.execute("CREATE TABLE source_counts (dataset TEXT PRIMARY KEY, template_count INTEGER, valid_count INTEGER)")
                insert_sql = """INSERT OR REPLACE INTO templates
                    (template_id, score, template, reactions, radii, reactions_count, ecs, ecs_count,
                     radius_min, radius_max, valid, datasets)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""
                batch = []
                with gzip.open(self.source, "rt", encoding="utf-8-sig", newline="") as handle:
                    reader = csv.DictReader(handle)
                    columns = set(reader.fieldnames or [])
                    missing = sorted(EXPECTED_COLUMNS - columns)
                    if missing:
                        raise ValueError("RetroRules CSV is missing required columns: " + ", ".join(missing))
                    for row in reader:
                        template_id = (row.get("TEMPLATE_ID") or "").strip()
                        template = (row.get("TEMPLATE") or "").strip()
                        if not template_id or not template:
                            continue
                        is_valid = (row.get("VALID") or "").strip().lower() == "true"
                        total += 1
                        valid_count += int(is_valid)
                        invalid_count += int(not is_valid)
                        sources = [s.strip().lower() for s in (row.get("DATASETS") or "").split(";") if s.strip()]
                        for source_name in sources:
                            counts = source_counts.setdefault(source_name, [0, 0])
                            counts[0] += 1
                            counts[1] += int(is_valid)
                        def as_int(value):
                            try:
                                return int(value) if value not in (None, "") else None
                            except (TypeError, ValueError):
                                return None
                        batch.append((
                            template_id, row.get("SCORE", ""), template, row.get("REACTIONS", ""),
                            row.get("RADII", ""), as_int(row.get("REACTIONS_COUNT")), row.get("ECS", ""),
                            as_int(row.get("ECS_COUNT")), as_int(row.get("RADIUS_MIN")),
                            as_int(row.get("RADIUS_MAX")), int(is_valid), ";".join(sources),
                        ))
                        if len(batch) >= batch_size:
                            conn.executemany(insert_sql, batch)
                            batch.clear()
                    if batch:
                        conn.executemany(insert_sql, batch)
                conn.execute("CREATE INDEX idx_templates_dataset ON templates(datasets)")
                conn.execute("CREATE INDEX idx_templates_valid ON templates(valid)")
                conn.execute("CREATE INDEX idx_templates_reactions_count ON templates(reactions_count)")
                conn.executemany("INSERT INTO source_counts VALUES (?, ?, ?)",
                                 [(name, counts[0], counts[1]) for name, counts in source_counts.items()])
                metadata = {
                    **current_fp,
                    "template_count": str(total),
                    "valid_count": str(valid_count),
                    "invalid_count": str(invalid_count),
                    "source_sha256": _sha256(self.source),
                    "indexed_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                }
                conn.executemany("INSERT INTO metadata VALUES (?, ?)", metadata.items())
                conn.commit()
                conn.close()
                os.replace(tmp_path, self.database)
            except Exception:
                try:
                    conn.close()
                except Exception:
                    pass
                tmp_path.unlink(missing_ok=True)
                raise
            return self.status()

    @staticmethod
    def _row_dict(row: sqlite3.Row, validate_smarts: bool = False) -> dict[str, Any]:
        result = dict(row)
        result["valid"] = bool(result["valid"])
        if validate_smarts:
            try:
                from rdkit.Chem import AllChem
                reaction = AllChem.ReactionFromSmarts(result["template"])
                if reaction is None:
                    raise ValueError("RDKit returned no reaction object")
                warnings, errors = reaction.Validate(silent=True)
                result["smarts_parse_ok"] = errors == 0
                result["smarts_validation_warnings"] = warnings
                result["smarts_validation_errors"] = errors
            except Exception as exc:  # noqa: BLE001
                result["smarts_parse_ok"] = False
                result["smarts_validation_error"] = f"{type(exc).__name__}: {exc}"
        return result

    def search(self, q: str | None = None, template_id: str | None = None, dataset: str | None = None,
               valid: bool | None = None, limit: int = 20, offset: int = 0,
               validate_smarts: bool = False) -> dict[str, Any]:
        if not self.status().get("indexed"):
            raise RuntimeError("RetroRules dataset has not been indexed. POST /api/reactions/index first.")
        if not 1 <= limit <= 100 or offset < 0:
            raise ValueError("limit must be 1-100 and offset must be non-negative")
        clauses, params = [], []
        if template_id:
            clauses.append("template_id = ?")
            params.append(template_id)
        if dataset:
            if len(dataset) > 50:
                raise ValueError("dataset must be 50 characters or fewer")
            clauses.append("instr(';' || datasets || ';', ';' || lower(?) || ';') > 0")
            params.append(dataset.strip())
        if valid is not None:
            clauses.append("valid = ?")
            params.append(int(valid))
        if q:
            if len(q) > 200:
                raise ValueError("q must be 200 characters or fewer")
            clauses.append("(template_id LIKE ? OR template LIKE ? OR reactions LIKE ?)")
            needle = f"%{q.strip()}%"
            params.extend([needle, needle, needle])
        where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
        with closing(_connect(self.database)) as conn:
            total = conn.execute("SELECT COUNT(*) FROM templates" + where, params).fetchone()[0]
            rows = conn.execute(
                "SELECT * FROM templates" + where + " ORDER BY template_id LIMIT ? OFFSET ?",
                [*params, limit, offset],
            ).fetchall()
        return {"total": total, "limit": limit, "offset": offset,
                "results": [self._row_dict(row, validate_smarts) for row in rows]}

    def get(self, template_id: str, validate_smarts: bool = True) -> dict[str, Any] | None:
        result = self.search(template_id=template_id, limit=1, validate_smarts=validate_smarts)
        return result["results"][0] if result["results"] else None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build or inspect the Linkwright RetroRules index")
    parser.add_argument("action", choices=["index", "status"])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    store = ReactionTemplateStore()
    result = store.build_index(force=args.force) if args.action == "index" else store.status()
    import json
    print(json.dumps(result, indent=2))
