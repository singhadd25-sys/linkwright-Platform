"""Hybrid molecule generation.

Two strategies, both fully rule-based and deterministic (no black-box model), so every
hybrid can be explained bond by bond:

1. REACTION RULES - a named, literature-style coupling (amide, CuAAC triazole, ...) is applied
   when parent A carries one functional group and parent B carries the partner group.
2. LINKER JOINING - a hydrogen on parent A and a hydrogen on parent B are replaced by a short
   linker (direct bond, alkyl chain, piperazine, amide, triazole, ...).

Each product keeps a record of which atoms came from which parent and which atoms form the new
connection, so the explanation and the highlighted drawing are generated from real data.
"""
from __future__ import annotations

from dataclasses import dataclass, field

_AMINE = "[NX3;H2,H1;!$(NC=[O,S,N]);!$(NS=O)"  # nucleophilic N-H (not amide / sulfonamide)

# id, display name, reactant SMARTS pair >> product, role names, description, conditions, reference
RULES = [
    dict(
        id="amide",
        name="Amide coupling",
        smarts=f"[CX3:1](=[O:2])[OX2H1].{_AMINE}:3]>>[C:1](=[O:2])[N:3]",
        roles=("carboxylic acid", "amine"),
        bond="C(=O)-N",
        loss="H2O",
        conditions="Typical: HATU or EDC/HOBt, DIPEA, DMF, room temperature.",
        ref="Montalbetti & Falque, Tetrahedron 2005, 61, 10827.",
    ),
    dict(
        id="ester",
        name="Esterification",
        smarts="[CX3:1](=[O:2])[OX2H1].[OX2H1;!$(OC=[O,S]):3]>>[C:1](=[O:2])[O:3]",
        roles=("carboxylic acid", "alcohol/phenol"),
        bond="C(=O)-O",
        loss="H2O",
        conditions="Typical: EDC/DMAP (Steglich) or acid chloride route.",
        ref="Neises & Steglich, Angew. Chem. Int. Ed. 1978, 17, 522.",
    ),
    dict(
        id="ether",
        name="Williamson ether synthesis",
        smarts="[CX4:1][Cl,Br,I].[OX2H1;!$(OC=[O,S]):2]>>[C:1][O:2]",
        roles=("alkyl halide", "alcohol/phenol"),
        bond="C-O",
        loss="HX",
        conditions="Typical: K2CO3 or NaH, DMF or acetonitrile.",
        ref="",
    ),
    dict(
        id="n_alkyl",
        name="N-alkylation",
        smarts=f"[CX4:1][Cl,Br,I].{_AMINE}:2]>>[C:1][N:2]",
        roles=("alkyl halide", "amine"),
        bond="C-N",
        loss="HX",
        conditions="Typical: K2CO3 or DIPEA, DMF or acetonitrile.",
        ref="",
    ),
    dict(
        id="reductive_amination",
        name="Reductive amination",
        smarts=f"[CX3H1:1](=O)[*:3].{_AMINE}:2]>>[C:1]([*:3])[N:2]",
        roles=("aldehyde", "amine"),
        bond="C-N",
        loss="H2O (then reduction)",
        conditions="Typical: NaBH(OAc)3, DCE or THF, acetic acid.",
        ref="Abdel-Magid et al., J. Org. Chem. 1996, 61, 3849.",
    ),
    dict(
        id="sulfonamide",
        name="Sulfonamide formation",
        smarts=f"[SX4:1](=[O:2])(=[O:3])[Cl].{_AMINE}:4]>>[S:1](=[O:2])(=[O:3])[N:4]",
        roles=("sulfonyl chloride", "amine"),
        bond="S(=O)(=O)-N",
        loss="HCl",
        conditions="Typical: pyridine or Et3N, DCM, 0 C to room temperature.",
        ref="",
    ),
    dict(
        id="urea",
        name="Urea formation",
        smarts=f"[N:1]=[C:2]=[O:3].{_AMINE}:4]>>[N:1][C:2](=[O:3])[N:4]",
        roles=("isocyanate", "amine"),
        bond="N-C(=O)-N",
        loss="none (addition)",
        conditions="Typical: DCM or THF, room temperature.",
        ref="",
    ),
    dict(
        id="cuaac",
        name="CuAAC 'click' triazole",
        smarts="[#6:1][C:2]#[CH:3].[N;+0:5]=[N+:6]=[N-:7]>>[#6:1][c:2]1[cH:3][n:5][n:6][n:7]1",
        roles=("terminal alkyne", "organic azide"),
        bond="1,4-disubstituted 1,2,3-triazole ring",
        loss="none (cycloaddition)",
        conditions="Typical: CuSO4/sodium ascorbate, tBuOH-water, room temperature.",
        ref="Rostovtsev et al., Angew. Chem. Int. Ed. 2002, 41, 2596; Tornoe et al., J. Org. Chem. 2002, 67, 3057.",
    ),
    dict(
        id="suzuki",
        name="Suzuki-Miyaura coupling",
        smarts="[c:1][Cl,Br,I].[c:2][BX3]([OX2H1])[OX2H1]>>[c:1][c:2]",
        roles=("aryl halide", "aryl boronic acid"),
        bond="aryl-aryl C-C",
        loss="HX + boric acid",
        conditions="Typical: Pd(PPh3)4 or Pd(dppf)Cl2, K2CO3, dioxane-water, 80-100 C.",
        ref="Miyaura & Suzuki, Chem. Rev. 1995, 95, 2457.",
    ),
]

# '*' marks the two attachment points. Order matters for asymmetric linkers, so both
# orientations are tried.
LINKERS = [
    dict(id="direct", name="Direct bond", smiles=None),
    dict(id="ethylene", name="Ethylene", smiles="*CC*"),
    dict(id="propylene", name="Propylene", smiles="*CCC*"),
    dict(id="piperazine", name="Piperazine", smiles="*N1CCN(CC1)*"),
    dict(id="amide", name="Amide", smiles="*C(=O)N*"),
    dict(id="ether_glycol", name="Ethylene glycol ether", smiles="*OCCO*"),
    dict(id="urea", name="Urea", smiles="*NC(=O)N*"),
    dict(id="sulfonamide", name="Sulfonamide", smiles="*S(=O)(=O)N*"),
    dict(id="triazole", name="1,4-Triazole (via CH2)", smiles="*Cn1cc(*)nn1"),
]
DEFAULT_LINKERS = [l["id"] for l in LINKERS]


@dataclass
class Hybrid:
    smiles: str
    method: str  # "reaction" or "linker"
    rule_id: str
    rule_name: str
    explanation: str
    highlight_atoms: list[int] = field(default_factory=list)
    mol: object = None
    details: dict = field(default_factory=dict)


# --------------------------------------------------------------------------------------
# Reaction-rule mode
# --------------------------------------------------------------------------------------
def _fix_triazole_charges(mol):
    for a in mol.GetAtoms():
        if a.GetSymbol() == "N" and a.GetIsAromatic() and a.GetFormalCharge() != 0:
            a.SetFormalCharge(0)
            a.SetNoImplicit(False)


def _valid_single_product(mol):
    from rdkit import Chem

    try:
        _fix_triazole_charges(mol)
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001
        return False
    return len(Chem.GetMolFrags(mol)) == 1


def enumerate_reactions(a, b, name_a="A", name_b="B", max_per_rule=20):
    """Apply every rule in both parent orientations. Returns list[Hybrid]."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    found: dict[str, Hybrid] = {}
    parents = {"A": (a, name_a), "B": (b, name_b)}
    for rule in RULES:
        rxn = AllChem.ReactionFromSmarts(rule["smarts"])
        maps_r0 = {at.GetAtomMapNum() for at in rxn.GetReactantTemplate(0).GetAtoms() if at.GetAtomMapNum()}
        maps_r1 = {at.GetAtomMapNum() for at in rxn.GetReactantTemplate(1).GetAtoms() if at.GetAtomMapNum()}
        for first, second in (("A", "B"), ("B", "A")):
            m1, n1 = parents[first]
            m2, n2 = parents[second]
            try:
                outs = rxn.RunReactants((m1, m2), maxProducts=200)
            except Exception:  # noqa: BLE001
                continue
            kept = 0
            for tup in outs:
                p = tup[0]
                if not _valid_single_product(p):
                    continue
                smi = Chem.MolToSmiles(p)
                if smi in found:
                    continue
                site1, site2, centre = [], [], []
                for at in p.GetAtoms():
                    if not at.HasProp("old_mapno"):
                        continue
                    centre.append(at.GetIdx())
                    idx = at.GetIntProp("react_atom_idx") if at.HasProp("react_atom_idx") else None
                    mp = at.GetIntProp("old_mapno")
                    if mp in maps_r0 and idx is not None:
                        site1.append(idx)
                    elif mp in maps_r1 and idx is not None:
                        site2.append(idx)
                text = (
                    f"{rule['name']}. Parent {first} ({n1}) supplies the {rule['roles'][0]} "
                    f"(reactive atom index {site1 or 'n/a'} in the input structure); parent {second} ({n2}) "
                    f"supplies the {rule['roles'][1]} (reactive atom index {site2 or 'n/a'}). "
                    f"The new connection is a {rule['bond']} linkage; by-product: {rule['loss']}. "
                    f"{rule['conditions']}"
                )
                found[smi] = Hybrid(
                    smiles=smi,
                    method="reaction",
                    rule_id=rule["id"],
                    rule_name=rule["name"],
                    explanation=text,
                    highlight_atoms=centre,
                    mol=p,
                    details=dict(
                        first_parent=first,
                        second_parent=second,
                        sites_first=site1,
                        sites_second=site2,
                        reference=rule["ref"],
                        bond=rule["bond"],
                    ),
                )
                kept += 1
                if kept >= max_per_rule:
                    break
    return list(found.values())


# --------------------------------------------------------------------------------------
# Linker mode
# --------------------------------------------------------------------------------------
def candidate_sites(mol, requested=None, limit=6):
    """Atoms that carry at least one hydrogen. Symmetry-unique and peripheral first."""
    from rdkit import Chem

    n = mol.GetNumAtoms()
    if requested:
        sites = [i for i in requested if 0 <= i < n and mol.GetAtomWithIdx(i).GetTotalNumHs() > 0]
        return sites
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    dm = Chem.GetDistanceMatrix(mol)
    seen, scored = set(), []
    for at in mol.GetAtoms():
        i = at.GetIdx()
        if at.GetTotalNumHs() == 0 or at.GetFormalCharge() != 0 or ranks[i] in seen:
            continue
        seen.add(ranks[i])
        scored.append((float(dm[i].mean()), i))
    scored.sort(reverse=True)  # most peripheral first
    return [i for _, i in scored[:limit]]


def _join(a, ia, b, ib, linker_smiles, flip=False):
    from rdkit import Chem

    na, nb = a.GetNumAtoms(), b.GetNumAtoms()
    combo = Chem.CombineMols(a, b)
    for at in combo.GetAtoms():
        at.SetIntProp("lw_role", 0 if at.GetIdx() < na else 1)
    if linker_smiles is None:
        rw = Chem.RWMol(combo)
        rw.AddBond(ia, na + ib, Chem.BondType.SINGLE)
        junction = [ia, na + ib]
        linker_atoms: list[int] = []
        dummies: list[int] = []
    else:
        lm = Chem.MolFromSmiles(linker_smiles)
        dum = [x.GetIdx() for x in lm.GetAtoms() if x.GetAtomicNum() == 0]
        if len(dum) != 2:
            raise ValueError("linker needs exactly two attachment points")
        if flip:
            dum = dum[::-1]
        nbrs = [lm.GetAtomWithIdx(d).GetNeighbors()[0].GetIdx() for d in dum]
        combo = Chem.CombineMols(combo, lm)
        off = na + nb
        rw = Chem.RWMol(combo)
        rw.AddBond(ia, off + nbrs[0], Chem.BondType.SINGLE)
        rw.AddBond(na + ib, off + nbrs[1], Chem.BondType.SINGLE)
        junction = [ia, na + ib]
        linker_atoms = [off + i for i in range(lm.GetNumAtoms()) if i not in dum]
        dummies = [off + d for d in dum]
    keep = set(junction) | set(linker_atoms)
    for at in rw.GetAtoms():
        at.SetIntProp("lw_keep", 1 if at.GetIdx() in keep else 0)
    for d in sorted(dummies, reverse=True):
        rw.RemoveAtom(d)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    highlight = [at.GetIdx() for at in mol.GetAtoms() if at.GetIntProp("lw_keep") == 1]
    return mol, highlight


def enumerate_linkers(a, b, name_a="A", name_b="B", linker_ids=None, sites_a=None, sites_b=None,
                      max_sites=6, max_mw=800.0):
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    ids = linker_ids or DEFAULT_LINKERS
    linkers = [l for l in LINKERS if l["id"] in ids]
    sa = candidate_sites(a, sites_a, max_sites)
    sb = candidate_sites(b, sites_b, max_sites)
    found: dict[str, Hybrid] = {}
    for lk in linkers:
        flips = (False,) if lk["smiles"] is None else (False, True)
        for ia in sa:
            for ib in sb:
                for flip in flips:
                    try:
                        mol, hl = _join(a, ia, b, ib, lk["smiles"], flip)
                    except Exception:  # noqa: BLE001
                        continue
                    if len(Chem.GetMolFrags(mol)) != 1 or Descriptors.MolWt(mol) > max_mw:
                        continue
                    smi = Chem.MolToSmiles(mol)
                    if smi in found:
                        continue
                    sym_a = a.GetAtomWithIdx(ia).GetSymbol()
                    sym_b = b.GetAtomWithIdx(ib).GetSymbol()
                    what = "a direct bond" if lk["smiles"] is None else f"the '{lk['name']}' linker"
                    text = (
                        f"Linker joining. A hydrogen on atom {ia} ({sym_a}) of parent A ({name_a}) and a hydrogen "
                        f"on atom {ib} ({sym_b}) of parent B ({name_b}) are replaced by {what}"
                        f"{' (reversed orientation)' if flip else ''}. Both parent scaffolds are otherwise "
                        "unchanged. This is a design-level construction; check synthetic feasibility."
                    )
                    found[smi] = Hybrid(
                        smiles=smi,
                        method="linker",
                        rule_id=f"linker:{lk['id']}",
                        rule_name=f"Linker: {lk['name']}",
                        explanation=text,
                        highlight_atoms=hl,
                        mol=mol,
                        details=dict(site_a=ia, site_b=ib, linker=lk["id"], reversed=flip),
                    )
    return list(found.values())


def generate(a, b, mode="both", name_a="A", name_b="B", linker_ids=None, sites_a=None, sites_b=None):
    out: list[Hybrid] = []
    if mode in ("reactions", "both"):
        out += enumerate_reactions(a, b, name_a, name_b)
    if mode in ("linkers", "both"):
        seen = {h.smiles for h in out}
        out += [h for h in enumerate_linkers(a, b, name_a, name_b, linker_ids, sites_a, sites_b)
                if h.smiles not in seen]
    return out
