"""Job pipeline: parse -> generate -> properties -> cross-reference -> MD package -> report."""
from __future__ import annotations

import json
import os
import zipfile
from datetime import datetime, timezone

VERSION = "0.4.2"
STAGES = ["parse", "generate", "properties", "xref", "md", "report"]


def _props_block(mol, smiles=None):
    from . import properties as P

    ph = P.physchem(mol)
    return dict(props=ph, rules=P.drug_likeness(ph, mol), alerts=P.alerts(mol))


def run(job_dir: str, inputs: dict, params: dict, progress) -> dict:
    """inputs: {A: (filename, bytes|None, smiles|None), B: ...}. progress(stage, state, message)."""
    import rdkit
    from rdkit import Chem

    from . import chem_io, gromacs, hybrid, properties as P, report, xref, surechembl

    res: dict = dict(
        job_id=os.path.basename(job_dir),
        created=report.now(),
        params=params,
        software=dict(linkwright=VERSION, rdkit=rdkit.__version__),
        warnings=[],
    )

    # 1. parse ---------------------------------------------------------------------
    progress("parse", "running", "Reading parent ligands")
    mols, names = {}, {}
    for k in ("A", "B"):
        fn, data, smi = inputs[k]
        mol, name, warns = chem_io.parse_ligand(fn or f"parent{k}", data, smi)
        mols[k], names[k] = mol, name
        res["warnings"] += warns
    res["parents"] = {}
    for k in ("A", "B"):
        blk = _props_block(mols[k])
        res["parents"][k] = dict(name=names[k], smiles=Chem.MolToSmiles(mols[k]),
                                 svg=chem_io.mol_svg(mols[k]), **blk)
    progress("parse", "done", "Parents read")

    # 2. generate ------------------------------------------------------------------
    progress("generate", "running", "Building hybrids")
    cand = hybrid.generate(
        mols["A"], mols["B"], params["mode"], names["A"], names["B"],
        params.get("linkers"), params.get("sites_a"), params.get("sites_b"),
    )
    if not cand:
        progress("generate", "failed", "No hybrid could be built")
        res["hybrids"] = []
        res["warnings"].append(
            "No hybrid could be built. In reaction mode the parents need complementary functional groups "
            "(for example a carboxylic acid and an amine). Try linker mode, which joins at any C-H or N-H site."
        )
        raise RuntimeError(res["warnings"][-1])
    progress("generate", "done", f"{len(cand)} candidate hybrids")

    # 3. properties + ranking ------------------------------------------------------
    progress("properties", "running", "Computing properties and drug-likeness")
    scored = []
    for h in cand:
        blk = _props_block(h.mol)
        n_pass = sum(1 for v in blk["rules"].values() if v["passed"])
        score = blk["props"]["qed"] + 0.05 * n_pass - 0.1 * len(blk["alerts"]["PAINS"])
        scored.append((score, h, blk))
    scored.sort(key=lambda t: t[0], reverse=True)
    keep = scored[: params["max_hybrids"]]
    hybrids = []
    for rank, (score, h, blk) in enumerate(keep, 1):
        try:
            ik = Chem.MolToInchiKey(h.mol)
        except Exception:  # noqa: BLE001
            ik = ""
        hybrids.append(dict(
            rank=rank, id=f"H{rank:02d}", smiles=h.smiles, inchikey=ik, method=h.method,
            rule_id=h.rule_id, rule_name=h.rule_name, explanation=h.explanation,
            details=h.details, svg=chem_io.mol_svg(h.mol, h.highlight_atoms),
            score=round(score, 3), admet_ai=None, linkouts=P.link_outs(h.smiles),
            xref=None, **blk,
        ))
    mol_by_id = {hb["id"]: keep[i][1].mol for i, hb in enumerate(hybrids)}
    progress("properties", "done", f"{len(hybrids)} hybrids ranked")

    # 4. cross-reference -----------------------------------------------------------
    progress("xref", "running", "Querying PubChem and ChEMBL")
    for k in ("A", "B"):
        res["parents"][k]["xref"] = xref.check_parent(res["parents"][k]["smiles"])
    for hb in hybrids[: params["xref_top"]]:
        hb["xref"] = dict(pubchem=xref.pubchem_lookup(hb["smiles"]), chembl=xref.chembl_similar(hb["smiles"]))
        if os.environ.get("LINKWRIGHT_SURECHEMBL", "1").strip().lower() not in {"0", "false", "off"}:
            hb["xref"]["surechembl_identity"] = surechembl.structure_search(hb["smiles"], "identity", 20)
            hb["xref"]["surechembl_similarity"] = surechembl.structure_search(hb["smiles"], "similarity", 20)
        hb["admet_ai"] = P.admet_ai_predict(hb["smiles"])
    unavailable = any(
        hb["xref"] and (
            hb["xref"]["pubchem"]["status"] != "ok" or hb["xref"]["chembl"]["status"] != "ok"
            or any(v.get("status") not in {"ok"} for k, v in hb["xref"].items() if k.startswith("surechembl_"))
        )
        for hb in hybrids if hb["xref"]
    )
    if unavailable:
        res["warnings"].append("One or more database queries failed. Re-run later; unreachable is not 'not found'.")
    progress("xref", "done", "Database checks finished")
    res["hybrids"] = hybrids

    # 5. MD package ---------------------------------------------------------------
    progress("md", "running", "Preparing GROMACS packages")
    top_n = min(params["md_top"], len(hybrids))
    pkg_path = os.path.join(job_dir, "gromacs_package.zip")
    selected = []
    with zipfile.ZipFile(pkg_path, "w", zipfile.ZIP_DEFLATED) as outer:
        for hb in hybrids[:top_n]:
            try:
                blob = gromacs.build_package(mol_by_id[hb["id"]], hb["id"], hb["smiles"], params["md_ns"], params["temp"])
                outer.writestr(f"{hb['id']}.zip", blob)
                selected.append(hb["id"])
            except Exception as e:  # noqa: BLE001
                res["warnings"].append(f"{hb['id']}: GROMACS package skipped ({e}).")
    res["gromacs"] = dict(selected=selected, ns=params["md_ns"], temp=params["temp"])
    progress("md", "done", f"{len(selected)} package(s)")

    # 6. report + exports ----------------------------------------------------------
    progress("report", "running", "Writing report")
    sdf = os.path.join(job_dir, "hybrids.sdf")
    w = Chem.SDWriter(sdf)
    for hb in hybrids:
        m = Chem.Mol(mol_by_id[hb["id"]])
        m.SetProp("_Name", hb["id"])
        m.SetProp("Rule", hb["rule_name"])
        m.SetProp("Explanation", hb["explanation"])
        w.write(m)
    w.close()
    with open(os.path.join(job_dir, "report.md"), "w", encoding="utf-8") as f:
        f.write(report.build_markdown(res))
    report.build_docx(res, os.path.join(job_dir, "report.docx"))
    with open(os.path.join(job_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    progress("report", "done", "Report ready")
    return res
