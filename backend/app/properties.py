"""Physicochemical properties, drug-likeness rules, structural alerts, optional ML ADMET.

Everything here is computed locally with RDKit (open source). Optional ADMET-AI predictions are
computed locally too if the `admet_ai` package is installed. External web servers
(SwissADME, pkCSM, ADMETlab 3.0, ProTox 3.0) have no stable public APIs and restrict
automated use, so Linkwright offers hand-off links instead of scraping them.
"""
from __future__ import annotations

import os
import sys
from urllib.parse import quote

_ADMET_MODEL = None
_ADMET_TRIED = False

LINK_OUTS = [
    dict(name="SwissADME", url="http://www.swissadme.ch/", note="Paste the SMILES in the input box."),
    dict(name="pkCSM", url="https://biosig.lab.uq.edu.au/pkcsm/prediction", note="Choose a prediction mode, paste SMILES."),
    dict(name="ADMETlab 3.0", url="https://admetlab3.scbdd.com/", note="Use the single-molecule or batch tool."),
    dict(name="ProTox 3.0", url="https://tox.charite.de/protox3/", note="Toxicity endpoints; paste SMILES."),
]


def _sascore(mol):
    try:
        from rdkit.Chem import RDConfig

        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer  # type: ignore

        return round(float(sascorer.calculateScore(mol)), 2)
    except Exception:  # noqa: BLE001
        return None


def physchem(mol) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors as rd

    return dict(
        formula=rd.CalcMolFormula(mol),
        mw=round(Descriptors.MolWt(mol), 2),
        logp=round(Crippen.MolLogP(mol), 2),
        hbd=rd.CalcNumHBD(mol),
        hba=rd.CalcNumHBA(mol),
        tpsa=round(rd.CalcTPSA(mol), 2),
        rotb=rd.CalcNumRotatableBonds(mol),
        fsp3=round(rd.CalcFractionCSP3(mol), 2),
        rings=rd.CalcNumRings(mol),
        aromatic_rings=rd.CalcNumAromaticRings(mol),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        mr=round(Crippen.MolMR(mol), 2),
        charge=Chem.GetFormalCharge(mol),
        qed=round(QED.qed(mol), 3),
        sa_score=_sascore(mol),
    )


def drug_likeness(p: dict, mol) -> dict:
    """Published rule sets, evaluated on the computed properties."""
    carbons = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    hetero = sum(1 for a in mol.GetAtoms() if a.GetSymbol() not in ("C", "H"))
    lip_viol = sum([p["mw"] > 500, p["logp"] > 5, p["hbd"] > 5, p["hba"] > 10])
    rules = {
        "Lipinski (<=1 violation)": (lip_viol <= 1, f"{lip_viol} violation(s)"),
        "Veber (RotB<=10, TPSA<=140)": (p["rotb"] <= 10 and p["tpsa"] <= 140, f"RotB {p['rotb']}, TPSA {p['tpsa']}"),
        "Egan (logP<=5.88, TPSA<=131.6)": (p["logp"] <= 5.88 and p["tpsa"] <= 131.6, f"logP {p['logp']}, TPSA {p['tpsa']}"),
        "Ghose": (
            160 <= p["mw"] <= 480 and -0.4 <= p["logp"] <= 5.6 and 40 <= p["mr"] <= 130 and 20 <= p["heavy_atoms"] <= 70,
            f"MW {p['mw']}, logP {p['logp']}, MR {p['mr']}",
        ),
        "Muegge": (
            200 <= p["mw"] <= 600 and -2 <= p["logp"] <= 5 and p["tpsa"] <= 150 and p["rings"] <= 7
            and carbons > 4 and hetero > 1 and p["rotb"] <= 15 and p["hba"] <= 10 and p["hbd"] <= 5,
            f"MW {p['mw']}, logP {p['logp']}",
        ),
    }
    return {k: dict(passed=bool(v[0]), detail=v[1]) for k, v in rules.items()}


def alerts(mol) -> dict:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    out = {}
    for label, cat in (("PAINS", FilterCatalogParams.FilterCatalogs.PAINS), ("Brenk", FilterCatalogParams.FilterCatalogs.BRENK)):
        params = FilterCatalogParams()
        params.AddCatalog(cat)
        fc = FilterCatalog(params)
        out[label] = sorted({e.GetDescription() for e in fc.GetMatches(mol)})
    return out


def admet_ai_predict(smiles: str) -> dict | None:
    """Local ML ADMET predictions (Swanson et al., 2024) if admet_ai is installed."""
    global _ADMET_MODEL, _ADMET_TRIED
    if not _ADMET_TRIED:
        _ADMET_TRIED = True
        try:
            from admet_ai import ADMETModel  # type: ignore

            _ADMET_MODEL = ADMETModel()
        except Exception:  # noqa: BLE001
            _ADMET_MODEL = None
    if _ADMET_MODEL is None:
        return None
    try:
        preds = _ADMET_MODEL.predict(smiles=smiles)
        return {k: (round(float(v), 3) if isinstance(v, (int, float)) else v) for k, v in preds.items()}
    except Exception:  # noqa: BLE001
        return None


def link_outs(smiles: str) -> list[dict]:
    return [dict(l, smiles=smiles, copy=quote(smiles, safe="")) for l in LINK_OUTS]
