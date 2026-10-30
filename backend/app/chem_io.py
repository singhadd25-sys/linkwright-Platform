"""Ligand input parsing (SDF, MOL, MOL2, PDB, SMILES) and structure drawing."""
from __future__ import annotations

import os

MAX_BYTES = 5 * 1024 * 1024
SUPPORTED = (".sdf", ".sd", ".mol", ".mol2", ".pdb", ".smi", ".smiles", ".txt")


def parse_ligand(filename: str, data: bytes | None, smiles_text: str | None = None):
    """Return (mol, name, warnings). Raises ValueError with a user-readable message."""
    from rdkit import Chem

    warnings: list[str] = []
    mol = None
    name = ""

    if data:
        if len(data) > MAX_BYTES:
            raise ValueError(f"{filename}: file is larger than {MAX_BYTES // (1024 * 1024)} MB.")
        ext = os.path.splitext((filename or "").lower())[1]
        text = data.decode("utf-8", errors="replace")
        if ext in (".sdf", ".sd"):
            sup = Chem.SDMolSupplier()
            sup.SetData(text)
            mols = [m for m in sup if m is not None]
            if not mols:
                raise ValueError(f"{filename}: no readable molecule found in the SDF.")
            if len(mols) > 1:
                warnings.append(f"{filename}: {len(mols)} molecules found; only the first was used.")
            mol = mols[0]
        elif ext == ".mol":
            mol = Chem.MolFromMolBlock(text)
        elif ext == ".mol2":
            if text.count("@<TRIPOS>MOLECULE") > 1:
                warnings.append(f"{filename}: several molecules found; only the first was used.")
            mol = Chem.MolFromMol2Block(text)
        elif ext == ".pdb":
            mol = Chem.MolFromPDBBlock(text)
            warnings.append(
                f"{filename}: PDB files carry no bond orders, so the chemistry may be guessed wrongly. "
                "Prefer SDF, MOL2 or SMILES."
            )
        elif ext in (".smi", ".smiles", ".txt"):
            line = next((l.strip() for l in text.splitlines() if l.strip()), "")
            parts = line.split()
            if parts:
                mol = Chem.MolFromSmiles(parts[0])
                name = parts[1] if len(parts) > 1 else ""
        else:
            raise ValueError(f"{filename}: unsupported format. Use one of {', '.join(SUPPORTED)}.")
        if mol is None:
            raise ValueError(f"{filename}: RDKit could not read this structure.")
    elif smiles_text and smiles_text.strip():
        mol = Chem.MolFromSmiles(smiles_text.strip().split()[0])
        if mol is None:
            raise ValueError("The SMILES string could not be parsed.")
    else:
        raise ValueError("Provide a ligand file or a SMILES string for each parent.")

    mol = _standardise(mol, warnings)
    if not name and mol.HasProp("_Name"):
        name = mol.GetProp("_Name").strip()
    if not name:
        name = os.path.splitext(os.path.basename(filename or "ligand"))[0] or "ligand"
    return mol, name, warnings


def _standardise(mol, warnings):
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
        warnings.append("Salt or multi-fragment input: the largest fragment was kept.")
    try:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    except Exception:  # noqa: BLE001
        pass
    Chem.SanitizeMol(mol)
    return mol


def mol_svg(mol, highlight_atoms=None, size=(380, 280)) -> str:
    """Draw a molecule as SVG, optionally highlighting atoms (the reaction centre)."""
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    from rdkit import Chem

    m = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    opts = d.drawOptions()
    opts.clearBackground = False
    atoms = list(highlight_atoms or [])
    colours = {i: (0.94, 0.64, 0.01) for i in atoms}
    d.DrawMolecule(m, highlightAtoms=atoms, highlightAtomColors=colours, highlightBonds=[])
    d.FinishDrawing()
    return d.GetDrawingText()
