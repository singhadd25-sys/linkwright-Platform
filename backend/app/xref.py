"""Cross-referencing against public chemistry databases (PubChem PUG-REST, ChEMBL web services).

Every query records its source URL, UTC timestamp and outcome, so a report can show exactly
what was checked. A failed or rate-limited query is reported as 'unavailable' and is never
treated as 'not found'. 'No exact match' means no match at query time, not proof of novelty.
"""
from __future__ import annotations

import time
from datetime import datetime, timezone
from urllib.parse import quote

TIMEOUT = 25.0
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS = {"User-Agent": "Linkwright/0.1 (research tool; contact: set-in-config)"}


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def _client():
    import httpx

    return httpx.Client(timeout=TIMEOUT, headers=HEADERS, follow_redirects=True)


def _pubchem_cids(client, path: str, smiles: str, params=None):
    """POST the SMILES (safe for long strings), follow PubChem's ListKey if it answers 'Waiting'."""
    url = f"{PUBCHEM}/{path}"
    r = client.post(url, data={"smiles": smiles}, params=params)
    if r.status_code == 404:
        return []
    r.raise_for_status()
    js = r.json()
    for _ in range(6):
        if "Waiting" not in js:
            break
        key = js["Waiting"]["ListKey"]
        time.sleep(2)
        r = client.get(f"{PUBCHEM}/compound/listkey/{key}/cids/JSON")
        r.raise_for_status()
        js = r.json()
    return js.get("IdentifierList", {}).get("CID", [])


def pubchem_lookup(smiles: str, similarity=90, limit=5) -> dict:
    out = dict(source="PubChem PUG-REST", queried=_now(), status="ok")
    try:
        with _client() as c:
            ident = _pubchem_cids(c, "compound/fastidentity/smiles/cids/JSON", smiles,
                                  {"identity_type": "same_connectivity"})
            sim = _pubchem_cids(c, "compound/fastsimilarity_2d/smiles/cids/JSON", smiles,
                                {"Threshold": similarity, "MaxRecords": limit + 1})
            sim = [x for x in sim if x not in ident][:limit]
            out["identity_cids"] = ident[:limit]
            out["similar_cids"] = sim
            out["identity_urls"] = [f"https://pubchem.ncbi.nlm.nih.gov/compound/{x}" for x in ident[:limit]]
            out["similar_urls"] = [f"https://pubchem.ncbi.nlm.nih.gov/compound/{x}" for x in sim]
            out["similarity_threshold"] = similarity
            out["verdict"] = (
                "Exact connectivity match found in PubChem." if ident
                else "No exact connectivity match in PubChem at query time."
            )
    except Exception as e:  # noqa: BLE001
        out.update(status="unavailable", error=f"{type(e).__name__}: {e}"[:200],
                   verdict="PubChem could not be reached; result unknown.")
    return out


def pubchem_title(cid: int) -> dict:
    try:
        with _client() as c:
            r = c.get(f"{PUBCHEM}/compound/cid/{cid}/property/Title,MolecularFormula,InChIKey/JSON")
            r.raise_for_status()
            return r.json()["PropertyTable"]["Properties"][0]
    except Exception:  # noqa: BLE001
        return {}


def chembl_similar(smiles: str, similarity=85, limit=5) -> dict:
    out = dict(source="ChEMBL web services", queried=_now(), status="ok", similarity_threshold=similarity)
    try:
        with _client() as c:
            url = f"{CHEMBL}/similarity/{quote(smiles, safe='')}/{similarity}.json"
            r = c.get(url, params={"limit": limit})
            r.raise_for_status()
            hits = []
            for m in r.json().get("molecules", [])[:limit]:
                hits.append(dict(
                    chembl_id=m.get("molecule_chembl_id"),
                    name=m.get("pref_name"),
                    similarity=m.get("similarity"),
                    url=f"https://www.ebi.ac.uk/chembl/compound_report_card/{m.get('molecule_chembl_id')}/",
                ))
            out["hits"] = hits
            out["verdict"] = (
                f"{len(hits)} ChEMBL compound(s) at >= {similarity}% similarity." if hits
                else f"No ChEMBL compound at >= {similarity}% similarity at query time."
            )
    except Exception as e:  # noqa: BLE001
        out.update(status="unavailable", error=f"{type(e).__name__}: {e}"[:200],
                   verdict="ChEMBL could not be reached; result unknown.")
    return out


def check_parent(smiles: str) -> dict:
    """Confirm a parent ligand exists in PubChem and fetch its public name."""
    res = pubchem_lookup(smiles, limit=3)
    if res.get("identity_cids"):
        res["record"] = pubchem_title(res["identity_cids"][0])
    return res
