"""Explainable functional-group analysis and curated forward reaction execution.

RetroRules is a reaction-template/retrosynthesis resource. Its templates are
not blindly executed as two-parent synthesis rules. Instead, this module uses
RDKit SMARTS functional-group detection and the application's explicitly
curated forward reaction templates; RetroRules remains searchable provenance
and supporting reaction knowledge.
"""
from __future__ import annotations
from typing import Any

FUNCTIONAL_GROUPS = [
    ("carboxylic_acid", "Carboxylic acid", "[CX3](=O)[OX2H1]"),
    ("primary_amine", "Primary amine", "[NX3;H2;!$(NC=O);!$(NS=O)]"),
    ("secondary_amine", "Secondary amine", "[NX3;H1;!$(NC=O);!$(NS=O)]"),
    ("alcohol", "Alcohol", "[OX2H][CX4;!$(C=O)]"),
    ("phenol", "Phenol", "[OX2H]c"),
    ("aldehyde", "Aldehyde", "[CX3H1](=O)[#6]"),
    ("ketone", "Ketone", "[#6][CX3](=O)[#6]"),
    ("sulfonyl_chloride", "Sulfonyl chloride", "[SX4](=O)(=O)Cl"),
    ("isocyanate", "Isocyanate", "[NX2]=[CX2]=[OX1]"),
    ("alkyl_halide", "Alkyl halide", "[CX4][Cl,Br,I]"),
    ("aryl_halide", "Aryl halide", "[c][Cl,Br,I]"),
    ("boronic_acid", "Boronic acid", "[BX3]([OX2H])[OX2H]"),
    ("terminal_alkyne", "Terminal alkyne", "[C]#[CH]"),
    ("organic_azide", "Organic azide", "[N]=[N+]=[N-]"),
]


def detect_functional_groups(mol) -> list[dict[str, Any]]:
    from rdkit import Chem
    found = []
    for key, label, smarts in FUNCTIONAL_GROUPS:
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            continue
        matches = mol.GetSubstructMatches(query, uniquify=True, maxMatches=100)
        if matches:
            found.append({"id": key, "name": label, "smarts": smarts,
                          "matches": [list(m) for m in matches], "count": len(matches)})
    return found


def _rule_products(a, b, name_a: str, name_b: str, max_per_rule: int = 10):
    from . import hybrid
    from rdkit import Chem
    from rdkit.Chem import AllChem
    results = []
    for rule in hybrid.RULES:
        rxn = AllChem.ReactionFromSmarts(rule["smarts"])
        if rxn is None:
            continue
        for orientation, (left, right, left_name, right_name) in enumerate(
            ((a, b, name_a, name_b), (b, a, name_b, name_a))
        ):
            try:
                products = rxn.RunReactants((left, right), maxProducts=100)
            except Exception:
                continue
            count = 0
            for product_tuple in products:
                if not product_tuple:
                    continue
                product = Chem.Mol(product_tuple[0])
                try:
                    Chem.SanitizeMol(product)
                except Exception:
                    continue
                if len(Chem.GetMolFrags(product)) != 1:
                    continue
                smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                # Server-rendered 2D depiction for the Build Hybrids reaction preview.
                from rdkit.Chem.Draw import rdMolDraw2D
                drawer = rdMolDraw2D.MolDraw2DSVG(460, 280)
                drawer.DrawMolecule(product)
                drawer.FinishDrawing()
                product_svg = drawer.GetDrawingText()
                results.append({
                    "rule_id": rule["id"], "rule_name": rule["name"],
                    "reactant_roles": list(rule["roles"]),
                    "orientation": ["A", "B"] if orientation == 0 else ["B", "A"],
                    "reactant_names": [left_name, right_name],
                    "product_smiles": smiles,
                    "product_svg": product_svg,
                    "new_connection": rule["bond"],
                    "reported_byproduct": rule["loss"],
                    "byproduct_scope_note": "Idealized reaction metadata only; complete atom/electron balance, reagent-derived waste and work-up products are not simulated.",
                    "typical_conditions": rule["conditions"],
                    "activation_mechanism": {
                        "amide": "Carboxyl activation followed by amine nucleophilic acyl substitution",
                        "ester": "Acyl activation/substitution by alcohol or phenol oxygen",
                        "ether": "Base-assisted nucleophilic substitution at an alkyl electrophile",
                        "n_alkyl": "Nucleophilic substitution at an alkyl electrophile",
                        "reductive_amination": "Carbonyl-amine condensation followed by iminium reduction; this product graph omits the reduction step chemistry",
                        "sulfonamide": "Nucleophilic substitution at electrophilic sulfur(VI)",
                        "urea": "Nucleophilic addition of amine to isocyanate carbon",
                        "cuaac": "Cu(I)-catalyzed azide-alkyne cycloaddition (idealized product mapping)",
                        "suzuki": "Pd-catalyzed aryl cross-coupling (idealized product mapping)",
                    }.get(rule["id"], "Mechanism not annotated"),
                    "literature_reference": rule["ref"],
                    "validation": {"rdkit_sanitized": True, "single_connected_product": True, "canonical_smiles": smiles, "formula": __import__("rdkit.Chem", fromlist=["rdMolDescriptors"]).rdMolDescriptors.CalcMolFormula(product), "formal_charge": sum(at.GetFormalCharge() for at in product.GetAtoms())},
                    "evidence_class": "curated_forward_rule_computational_product",
                    "warning": "Computational product graph only; no experimental yield, selectivity, atom-balanced reagent/byproduct model, or synthetic feasibility is established.",
                })
                count += 1
                if count >= max_per_rule:
                    break
    unique = {}
    for item in results:
        unique[(item["rule_id"], item["product_smiles"])] = item
    return list(unique.values())


def analyze_pair(smiles_a: str, smiles_b: str, name_a: str = "Parent A", name_b: str = "Parent B") -> dict[str, Any]:
    from rdkit import Chem
    from .reaction_db import ReactionTemplateStore
    a = Chem.MolFromSmiles(smiles_a)
    b = Chem.MolFromSmiles(smiles_b)
    if a is None or b is None:
        raise ValueError("One or both input SMILES strings could not be parsed by RDKit.")
    groups_a = detect_functional_groups(a)
    groups_b = detect_functional_groups(b)
    products = _rule_products(a, b, name_a, name_b)
    try:
        db_status = ReactionTemplateStore().status()
    except Exception as exc:
        db_status = {"indexed": False, "message": f"Reaction database status unavailable: {exc}"}
    return {
        "parents": {
            "A": {"name": name_a, "canonical_smiles": Chem.MolToSmiles(a), "functional_groups": groups_a},
            "B": {"name": name_b, "canonical_smiles": Chem.MolToSmiles(b), "functional_groups": groups_b},
        },
        "compatible_forward_reactions": products,
        "reaction_count": len(products),
        "retrorules_database": {
            "indexed": bool(db_status.get("indexed")),
            "template_count": int(db_status.get("template_count", 0) or 0),
            "source_distribution": db_status.get("source_distribution", []),
            "role": "searchable reaction-template evidence; not automatically treated as a two-parent forward synthesis rule",
        },
        "interpretation": "Functional-group matches are substructure annotations. Candidate products are generated only by Linkwright's curated forward reaction SMARTS and pass RDKit sanitization/connectivity checks.",
        "limitations": [
            "RDKit sanitization does not prove experimental synthesis or stability.",
            "RetroRules templates may encode retrosynthetic or single-substrate transformations and are not blindly executed as forward two-parent reactions.",
            "Reaction conditions, selectivity, protecting groups, yield and substrate scope require literature review and experimental validation.",
        ],
    }
