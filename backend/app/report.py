"""Methods-and-results report in DOCX and Markdown, generated from the job result dict only.

Nothing is written that the pipeline did not compute or query; claims that cannot be verified
automatically (synthetic feasibility, biological activity) are stated as limitations.
"""
from __future__ import annotations

from datetime import datetime, timezone

REFERENCES = [
    "Landrum G. RDKit: Open-source cheminformatics. https://www.rdkit.org",
    "Kim S. et al. PubChem 2023 update. Nucleic Acids Res. 2023, 51, D1373-D1380.",
    "Zdrazil B. et al. The ChEMBL Database in 2023. Nucleic Acids Res. 2024, 52, D1180-D1192.",
    "Lipinski C.A. et al. Adv. Drug Deliv. Rev. 2001, 46, 3-26.",
    "Veber D.F. et al. J. Med. Chem. 2002, 45, 2615-2623.",
    "Baell J.B., Holloway G.A. J. Med. Chem. 2010, 53, 2719-2740 (PAINS).",
    "Bickerton G.R. et al. Nat. Chem. 2012, 4, 90-98 (QED).",
    "Abraham M.J. et al. GROMACS: High performance molecular simulations. SoftwareX 2015, 1-2, 19-25.",
    "Wang J. et al. Development and testing of a general Amber force field. J. Comput. Chem. 2004, 25, 1157-1174 (GAFF).",
    "Sousa da Silva A.W., Vranken W.F. ACPYPE. BMC Res. Notes 2012, 5, 367.",
    "SureChEMBL 2.0: an open database of chemical and biomedical annotations from patent literature. Nucleic Acids Res. 2026. https://doi.org/10.1093/nar/gkag1057.",
]


def _bits(h):
    p = h["props"]
    return (f"MW {p['mw']}, cLogP {p['logp']}, HBD {p['hbd']}, HBA {p['hba']}, TPSA {p['tpsa']}, "
            f"RotB {p['rotb']}, QED {p['qed']}")


def build_markdown(res: dict) -> str:
    L = []
    L.append(f"# Hybrid ligand design report\n\nJob {res['job_id']} - generated {res['created']}\n")
    L.append("## 1. Parent ligands\n")
    for k in ("A", "B"):
        p = res["parents"][k]
        L.append(f"- **Parent {k}: {p['name']}** - `{p['smiles']}` ({_bits(p)})")
        if p.get("xref"):
            L.append(f"  - PubChem: {p['xref'].get('verdict', 'n/a')} ({p['xref'].get('queried', '')})")
    L.append("\n## 2. Methods\n")
    L.append(
        f"Hybrids were enumerated with Linkwright {res['software']['linkwright']} using RDKit "
        f"{res['software']['rdkit']}. Mode: {res['params']['mode']}. Named coupling rules (amide, ester, "
        "ether, N-alkylation, reductive amination, sulfonamide, urea, CuAAC triazole, Suzuki-Miyaura) "
        "were applied as reaction SMARTS in both parent orientations, and/or a hydrogen on each parent was "
        "replaced by a linker. Products were sanitised, restricted to a single fragment and de-duplicated by "
        "canonical SMILES. Properties and drug-likeness rules (Lipinski, Veber, Egan, Ghose, Muegge) and "
        "PAINS/Brenk alerts were computed with RDKit. Identity and similarity searches used PubChem "
        "PUG-REST and ChEMBL web services; each query is time-stamped below.\n"
    )
    L.append("## 3. Results\n")
    for h in res["hybrids"]:
        L.append(f"### Hybrid {h['rank']}: {h['rule_name']}\n")
        L.append(f"- SMILES: `{h['smiles']}`\n- InChIKey: {h.get('inchikey', 'n/a')}\n- Properties: {_bits(h)}")
        L.append(f"- How it was formed: {h['explanation']}")
        if h["details"].get("reference"):
            L.append(f"- Method reference: {h['details']['reference']}")
        fails = [k for k, v in h["rules"].items() if not v["passed"]]
        L.append(f"- Drug-likeness: {'all rule sets passed' if not fails else 'failed: ' + ', '.join(fails)}")
        L.append(f"- Alerts: PAINS {h['alerts']['PAINS'] or 'none'}; Brenk {h['alerts']['Brenk'] or 'none'}")
        x = h.get("xref")
        if x:
            L.append(f"- PubChem: {x['pubchem']['verdict']} ({x['pubchem']['queried']})")
            L.append(f"- ChEMBL: {x['chembl']['verdict']} ({x['chembl']['queried']})")
            for key in ("surechembl_identity", "surechembl_similarity"):
                if x.get(key):
                    sx=x[key]
                    L.append(f"- SureChEMBL ({sx.get('search_type')}): {sx.get('verdict', sx.get('error', sx.get('status')))} ({sx.get('queried')})")
        if h.get("admet_ai"):
            L.append("- ADMET-AI predictions (local model): " + ", ".join(f"{k}={v}" for k, v in list(h["admet_ai"].items())[:12]))
        for imported in h.get("admet_imports", []):
            L.append(f"- Imported ADMET results from {imported.get('source')} ({imported.get('queried')}): " +
                     ", ".join(f"{k}={v}" for k, v in imported.get("endpoints", {}).items()))
        L.append("")
    L.append("## 4. MD preparation\n")
    g = res.get("gromacs", {})
    L.append(f"GROMACS input packages were prepared for hybrids {g.get('selected', [])} "
             f"({g.get('ns')} ns, {g.get('temp')} K). No simulation was performed by Linkwright.\n")
    L.append("## 5. Limitations\n")
    for t in LIMITS:
        L.append(f"- {t}")
    L.append("\n## References\n")
    for i, r in enumerate(REFERENCES, 1):
        L.append(f"{i}. {r}")
    for h in res["hybrids"]:
        if h["details"].get("reference"):
            L.append(f"- Reaction method: {h['details']['reference']}")
    return "\n".join(L) + "\n"


LIMITS = [
    "Hybrids are design-level constructions. Synthetic feasibility, stability and selectivity are not verified.",
    "'No exact match' in PubChem, ChEMBL or SureChEMBL means no match at query time; it does not prove novelty.",
    "SureChEMBL patent hits require manual review of claims, family, priority dates, legal status and disclosure context; this is not a freedom-to-operate or legal opinion.",
    "Drug-likeness rules and ML ADMET predictions are screening aids, not experimental measurements.",
    "Docking poses and molecular dynamics must be produced and validated separately by the researcher.",
    "Reaction reference citations are general method references; verify them before publication.",
]


def build_docx(res: dict, path: str):
    from docx import Document
    from docx.shared import Pt

    d = Document()
    d.styles["Normal"].font.name = "Calibri"
    d.styles["Normal"].font.size = Pt(11)
    d.add_heading("Hybrid ligand design report", 0)
    d.add_paragraph(f"Job {res['job_id']} - generated {res['created']}")

    d.add_heading("1. Parent ligands", 1)
    t = d.add_table(rows=1, cols=4)
    t.style = "Light Grid Accent 1"
    for c, txt in zip(t.rows[0].cells, ("Parent", "Name", "SMILES", "Properties")):
        c.text = txt
    for k in ("A", "B"):
        p = res["parents"][k]
        r = t.add_row().cells
        r[0].text, r[1].text, r[2].text, r[3].text = k, p["name"], p["smiles"], _bits(p)

    d.add_heading("2. Methods", 1)
    md = build_markdown(res)
    methods = md.split("## 2. Methods\n")[1].split("## 3. Results")[0].strip()
    d.add_paragraph(methods)

    d.add_heading("3. Results", 1)
    for h in res["hybrids"]:
        d.add_heading(f"Hybrid {h['rank']}: {h['rule_name']}", 2)
        d.add_paragraph(f"SMILES: {h['smiles']}")
        d.add_paragraph(f"InChIKey: {h.get('inchikey', 'n/a')}")
        d.add_paragraph(f"Properties: {_bits(h)}")
        d.add_paragraph(f"How it was formed: {h['explanation']}")
        if h["details"].get("reference"):
            d.add_paragraph(f"Method reference: {h['details']['reference']}")
        rt = d.add_table(rows=1, cols=3)
        rt.style = "Light Grid Accent 1"
        for c, txt in zip(rt.rows[0].cells, ("Rule set", "Result", "Detail")):
            c.text = txt
        for name, v in h["rules"].items():
            r = rt.add_row().cells
            r[0].text, r[1].text, r[2].text = name, "pass" if v["passed"] else "fail", v["detail"]
        d.add_paragraph(f"PAINS alerts: {', '.join(h['alerts']['PAINS']) or 'none'}. "
                        f"Brenk alerts: {', '.join(h['alerts']['Brenk']) or 'none'}.")
        x = h.get("xref")
        if x:
            d.add_paragraph(f"PubChem ({x['pubchem']['queried']}): {x['pubchem']['verdict']}")
            d.add_paragraph(f"ChEMBL ({x['chembl']['queried']}): {x['chembl']['verdict']}")
            for key in ("surechembl_identity", "surechembl_similarity"):
                if x.get(key):
                    sx=x[key]
                    d.add_paragraph(f"SureChEMBL {sx.get('search_type')} ({sx.get('queried')}): "
                                    f"{sx.get('verdict', sx.get('error', sx.get('status')))}")
        if h.get("admet_ai"):
            d.add_paragraph("ADMET-AI (local model): " + ", ".join(f"{k}={v}" for k, v in list(h["admet_ai"].items())[:12]))
        for imported in h.get("admet_imports", []):
            d.add_paragraph(f"Imported ADMET results from {imported.get('source')} ({imported.get('queried')}): " +
                            ", ".join(f"{k}={v}" for k, v in imported.get("endpoints", {}).items()))

    if res.get("admet_import_history"):
        d.add_heading("Imported ADMET data provenance", 1)
        for entry in res["admet_import_history"]:
            d.add_paragraph(f"{entry.get('source')} ({entry.get('queried')}): {entry.get('matched_count')} hybrid(s) matched; "
                            f"{entry.get('unmatched_count')} row(s) unmatched. Endpoints: {', '.join(entry.get('endpoints', []))}")

    d.add_heading("4. MD preparation", 1)
    g = res.get("gromacs", {})
    d.add_paragraph(f"GROMACS input packages were prepared for hybrids {g.get('selected', [])} "
                    f"({g.get('ns')} ns, {g.get('temp')} K). No simulation was performed by Linkwright.")
    d.add_heading("5. Limitations", 1)
    for t_ in LIMITS:
        d.add_paragraph(t_, style="List Bullet")
    d.add_heading("References", 1)
    for r in REFERENCES:
        d.add_paragraph(r, style="List Number")
    d.save(path)


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
