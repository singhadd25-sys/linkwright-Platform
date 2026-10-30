"""Structure quality checks and reproducible RDKit ligand preparation.

This module deliberately stops at an optimized, hydrogenated SDF. It does not
claim to create validated docking PDBQT or protein force-field topologies.
"""
from __future__ import annotations

import hashlib
import io
from datetime import datetime, timezone


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def prepare_ligand(filename: str, data: bytes | None = None, smiles: str = "", seed: int = 20261030) -> dict:
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from .chem_io import parse_ligand

    mol, name, warnings = parse_ligand(filename, data, smiles)
    mol = Chem.RemoveHs(mol)
    if mol.GetNumAtoms() == 0:
        raise ValueError("The ligand has no heavy atoms.")
    if mol.GetNumHeavyAtoms() > 500:
        raise ValueError("Ligands above 500 heavy atoms are not accepted by this preparation endpoint.")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchikey = Chem.MolToInchiKey(mol)
    work = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = False
    status = AllChem.EmbedMolecule(work, params)
    if status != 0:
        raise ValueError("RDKit could not generate a 3D conformer. Try another tautomer or stereochemical form.")
    method = "MMFF94s"
    if AllChem.MMFFHasAllMoleculeParams(work):
        props = AllChem.MMFFGetMoleculeProperties(work, mmffVariant="MMFF94s")
        ff = AllChem.MMFFGetMoleculeForceField(work, props)
    elif AllChem.UFFHasAllMoleculeParams(work):
        method = "UFF"
        ff = AllChem.UFFGetMoleculeForceField(work)
        warnings.append("MMFF94s parameters were incomplete; UFF was used instead.")
    else:
        raise ValueError("No complete MMFF94s or UFF parameterization is available for this molecule.")
    ff.Initialize()
    ff.Minimize(maxIts=1000)
    energy = float(ff.CalcEnergy())
    conf = work.GetConformer()
    for atom in work.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        if not all(abs(v) < 1e5 for v in (pos.x, pos.y, pos.z)):
            raise ValueError("Generated coordinates failed finite-coordinate validation.")
    work.SetProp("_Name", name)
    work.SetProp("Linkwright_Canonical_SMILES", canonical)
    work.SetProp("Linkwright_InChIKey", inchikey)
    work.SetProp("Linkwright_3D_Method", method)
    work.SetProp("Linkwright_Force_Field_Energy", f"{energy:.6f}")
    work.SetProp("Linkwright_Random_Seed", str(seed))
    work.SetProp("Linkwright_RDKit_Version", rdBase.rdkitVersion)
    work.SetProp("Linkwright_Prepared_UTC", datetime.now(timezone.utc).isoformat())
    buf = io.StringIO()
    writer = Chem.SDWriter(buf)
    writer.write(work)
    writer.close()
    sdf = buf.getvalue()
    sdf_bytes = sdf.encode("utf-8")
    return {
        "name": name,
        "canonical_smiles": canonical,
        "inchikey": inchikey,
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(float(Descriptors.MolWt(mol)), 4),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "stereochemistry_assigned": bool(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "optimization_method": method,
        "force_field_energy_kcal_mol": round(energy, 6),
        "random_seed": int(seed),
        "rdkit_version": rdBase.rdkitVersion,
        "input_sha256": _sha256(data if data is not None else smiles.encode("utf-8")),
        "prepared_sdf_sha256": _sha256(sdf_bytes),
        "prepared_sdf": sdf,
        "warnings": warnings + [
            "A single minimized conformer is not a conformational ensemble or a binding pose.",
            "This SDF is not PDBQT: docking atom types and partial charges still require a validated preparation toolchain.",
            "Review tautomer, protonation, stereochemistry and salt handling before scientific use.",
        ],
    }


def validate_protein_pdb(data: bytes, filename: str = "protein.pdb") -> dict:
    """Conservative PDB ATOM/HETATM inventory; not a repair or completeness validator."""
    if not data or len(data) > 100 * 1024 * 1024:
        raise ValueError("Protein PDB is empty or exceeds 100 MB.")
    text = data.decode("utf-8", errors="replace")
    atoms, residues, chains, altlocs, waters, hetero = 0, set(), set(), set(), 0, 0
    bad = []
    for line_no, line in enumerate(text.splitlines(), 1):
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) < 54:
                bad.append(line_no)
                continue
            try:
                xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
            except ValueError:
                bad.append(line_no)
                continue
            if not all(abs(v) < 1e5 for v in xyz):
                bad.append(line_no)
                continue
            atoms += 1
            chain = line[21:22].strip() or "_"
            resname = line[17:20].strip()
            resseq = line[22:27].strip()
            chains.add(chain)
            residues.add((chain, resname, resseq))
            if line[16:17].strip():
                altlocs.add(line[16:17].strip())
            if resname in {"HOH", "WAT", "H2O"}:
                waters += 1
            if line.startswith("HETATM"):
                hetero += 1
    if atoms == 0:
        raise ValueError("No valid ATOM/HETATM coordinate records were found in the PDB file.")
    warnings = [
        "This is an inventory check only; it does not detect or repair missing residues/atoms, validate protonation, or assign bond orders.",
        "Review chains, alternate conformations, waters, ions, cofactors and biological assembly before docking.",
    ]
    if bad:
        warnings.append(f"{len(bad)} malformed coordinate record(s) were skipped; first line numbers: {bad[:10]}.")
    if len(chains) > 1:
        warnings.append("Multiple chains are present; select the biologically relevant chain(s) explicitly.")
    if altlocs:
        warnings.append("Alternate-location identifiers are present; resolve them using a documented policy.")
    return {
        "filename": filename, "atom_records": atoms, "unique_residues": len(residues),
        "chains": sorted(chains), "alternate_location_ids": sorted(altlocs),
        "water_records": waters, "hetero_atom_records": hetero,
        "malformed_coordinate_lines": bad[:100], "input_sha256": _sha256(data),
        "warnings": warnings,
    }
