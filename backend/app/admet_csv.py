"""Import researcher-generated ADMET prediction CSVs without scraping third-party services."""
from __future__ import annotations
import csv
import io
import re

ID_KEYS = {"id", "hybrid_id", "compound_id", "name", "molecule", "molecule_id"}
SMILES_KEYS = {"smiles", "canonical_smiles", "isomeric_smiles"}
META_KEYS = ID_KEYS | SMILES_KEYS | {"source", "tool", "model", "date", "timestamp"}

def _key(s):
    return re.sub(r"[^a-z0-9]+", "_", (s or "").strip().lower()).strip("_")

def import_csv(content: bytes, hybrids: list[dict], source: str | None = None) -> dict:
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError as e:
        raise ValueError("CSV must be UTF-8 encoded.") from e
    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames:
        raise ValueError("CSV has no header row.")
    headers = {_key(h): h for h in reader.fieldnames}
    id_col = next((headers[k] for k in ID_KEYS if k in headers), None)
    smi_col = next((headers[k] for k in SMILES_KEYS if k in headers), None)
    if not id_col and not smi_col:
        raise ValueError("CSV needs a hybrid ID/name column or a SMILES column.")
    numeric_cols = [h for k, h in headers.items() if k not in META_KEYS]
    if not numeric_cols:
        raise ValueError("No ADMET endpoint columns were found.")
    by_id = {str(h.get("id", "")).strip().lower(): h for h in hybrids}
    by_smiles = {str(h.get("smiles", "")).strip(): h for h in hybrids}
    matched, unmatched, errors = {}, [], []
    for line_no, row in enumerate(reader, 2):
        ident = (row.get(id_col) or "").strip().lower() if id_col else ""
        smi = (row.get(smi_col) or "").strip() if smi_col else ""
        target = by_id.get(ident) if ident else None
        target = target or (by_smiles.get(smi) if smi else None)
        if target is None:
            unmatched.append(line_no)
            continue
        values = {}
        for col in numeric_cols:
            raw = (row.get(col) or "").strip()
            if not raw:
                continue
            try:
                values[_key(col)] = float(raw)
            except ValueError:
                values[_key(col)] = raw
        if values:
            matched[target["id"]] = values
    return {
        "source": (source or "User-imported CSV").strip()[:120],
        "matched": matched,
        "matched_count": len(matched),
        "unmatched_rows": unmatched,
        "unmatched_count": len(unmatched),
        "endpoints": sorted({k for v in matched.values() for k in v}),
    }
