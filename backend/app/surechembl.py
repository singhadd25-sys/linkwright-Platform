"""SureChEMBL 2.0 structure-search adapter.

Uses the published REST API. Results are time-stamped and never interpreted as proof of novelty.
The service is evolving; unexpected response shapes are reported as unavailable/unknown.
"""
from __future__ import annotations
import os, time
from datetime import datetime, timezone
from urllib.parse import quote

BASE = os.environ.get("LINKWRIGHT_SURECHEMBL_API", "https://www.surechembl.org/api").rstrip("/")
TIMEOUT = float(os.environ.get("LINKWRIGHT_SURECHEMBL_TIMEOUT", "30"))

def _now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

def _find_key(obj, names):
    if isinstance(obj, dict):
        for k,v in obj.items():
            if k.lower() in names and isinstance(v,(str,int)):
                return str(v)
        for v in obj.values():
            found=_find_key(v,names)
            if found: return found
    elif isinstance(obj,list):
        for v in obj:
            found=_find_key(v,names)
            if found: return found
    return None

def _extract_records(obj):
    if isinstance(obj, list) and (not obj or all(isinstance(x, dict) for x in obj)):
        return obj
    if isinstance(obj,dict):
        for k,v in obj.items():
            if k.lower() in {"results","compounds","hits","items","content","chemicals","searchresults","documents"} and isinstance(v,list):
                return v
        for v in obj.values():
            got=_extract_records(v)
            if got is not None: return got
    return None

def _chemical_ids(obj):
    found=set()
    if isinstance(obj,dict):
        for k,v in obj.items():
            if k.lower() in {"chemicalid","chemical_id","chem_id"} and isinstance(v,(str,int)):
                found.add(str(v))
            elif isinstance(v,(dict,list)):
                found.update(_chemical_ids(v))
    elif isinstance(obj,list):
        for v in obj: found.update(_chemical_ids(v))
    return sorted(found)


def structure_search(smiles: str, search_type="similarity", max_results=50, poll_seconds=1.5, max_polls=20):
    """Run identity/connectivity/similarity/substructure search, returning raw compact hits."""
    out={"source":"SureChEMBL 2.0 REST API","queried":_now(),"status":"unavailable",
         "search_type":search_type,"max_results":max_results,
         "interpretation":"Patent-mined chemical matches are leads for manual review, not a legal novelty opinion."}
    if search_type not in {"similarity","identity","connectivity","substructure"}:
        out["error"]="Unsupported SureChEMBL search type."
        return out
    try:
        import httpx
        headers={"User-Agent":"Linkwright/0.3 research prototype (contact configured by operator)"}
        with httpx.Client(timeout=TIMEOUT,headers=headers,follow_redirects=True) as c:
            payload={"StructureSearchRequest":{"struct":smiles,"structSearchType":search_type.upper(),
                      "maxResults":int(max_results),"saveSearch":"true"}}
            r=c.post(f"{BASE}/search/structure",json=payload)
            r.raise_for_status()
            first=r.json()
            search_hash=_find_key(first,{"hash","searchhash","search_hash"})
            if not search_hash:
                # Some API deployments may return a synchronous result directly.
                records=_extract_records(first)
                if records is not None:
                    out.update(status="ok",hits=records[:max_results],hit_count=len(records),
                               response_shape="synchronous")
                    return out
                out.update(status="unknown",error="SureChEMBL response did not expose a search hash or result list.")
                return out
            out["search_hash"]=search_hash
            final=None
            for _ in range(max_polls):
                st=c.get(f"{BASE}/search/{quote(search_hash,safe='')}/status")
                st.raise_for_status()
                status_json=st.json()
                status_text=str(_find_key(status_json,{"status","state"}) or "").lower()
                if any(x in status_text for x in ("complete","completed","finished","success","done")):
                    final=status_json
                    break
                if any(x in status_text for x in ("fail","error")):
                    out.update(status="unavailable",error="SureChEMBL search reported failure.")
                    return out
                time.sleep(poll_seconds)
            if final is None:
                out.update(status="pending",error="Search did not finish within the polling window.")
                return out
            rr=c.get(f"{BASE}/search/{quote(search_hash,safe='')}/results",
                     params={"page":0,"max_results":max_results})
            rr.raise_for_status()
            result_json=rr.json()
            records=_extract_records(result_json)
            if records is None:
                out.update(status="unknown",error="SureChEMBL completed but result records were not recognized.",
                           status_response=final)
                return out
            hits=records[:max_results]
            out.update(status="ok",hits=hits,hit_count=len(records),
                       result_metadata=result_json.get("data",{}) if isinstance(result_json,dict) else {})
            ids=_chemical_ids(hits)
            if ids:
                docs=c.post(f"{BASE}/search/documents_for_structures",
                            params=[("chemicalIds",x) for x in ids[:100]]+[("page","0"),("itemsPerPage",str(max_results))])
                docs.raise_for_status()
                out["patent_documents_response"]=docs.json()
                out["chemical_ids"]=ids[:100]
            out["verdict"]=f"{len(records)} patent-linked chemical record(s) returned by SureChEMBL."
            return out
    except Exception as e:
        out.update(status="unavailable",error=f"{type(e).__name__}: {e}"[:250],
                   verdict="SureChEMBL could not be queried; patent-search result unknown.")
        return out
