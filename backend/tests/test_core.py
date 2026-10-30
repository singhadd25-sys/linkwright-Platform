"""Run with:  cd backend && python -m pytest -q
Chemistry tests are skipped automatically if RDKit is not installed."""
import io
import os
import subprocess
import sys
import tempfile
import zipfile

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import rdkit as _rd  # noqa: F401
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False
needs_rdkit = pytest.mark.skipif(not HAVE_RDKIT, reason="RDKit not installed")


# ---------------------------------------------------------------- chemistry (RDKit)
@needs_rdkit
def test_amide_rule_builds_expected_product():
    from rdkit import Chem
    from app import hybrid

    acid = Chem.MolFromSmiles("OC(=O)c1ccccc1")        # benzoic acid
    amine = Chem.MolFromSmiles("NCc1ccccc1")           # benzylamine
    out = hybrid.enumerate_reactions(acid, amine, "acid", "amine")
    smiles = {h.smiles for h in out if h.rule_id == "amide"}
    assert Chem.MolToSmiles(Chem.MolFromSmiles("O=C(NCc1ccccc1)c1ccccc1")) in smiles


@needs_rdkit
def test_cuaac_gives_14_triazole():
    from rdkit import Chem
    from app import hybrid

    alkyne = Chem.MolFromSmiles("C#Cc1ccccc1")
    azide = Chem.MolFromSmiles("[N-]=[N+]=NCc1ccccc1")
    out = [h for h in hybrid.enumerate_reactions(alkyne, azide) if h.rule_id == "cuaac"]
    assert out, "CuAAC produced nothing"
    expected = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccc(cc1)Cn1cc(nn1)-c1ccccc1"))
    assert expected in {h.smiles for h in out}


@needs_rdkit
def test_linker_join_direct_and_piperazine():
    from rdkit import Chem
    from app import hybrid

    a = Chem.MolFromSmiles("c1ccccc1")
    b = Chem.MolFromSmiles("c1ccncc1")
    out = hybrid.enumerate_linkers(a, b, linker_ids=["direct", "piperazine"])
    smis = {h.smiles for h in out}
    assert any(Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts("c-c")) for s in smis)
    assert any(Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmarts("N1CCNCC1")) for s in smis)
    assert all(len(Chem.GetMolFrags(Chem.MolFromSmiles(s))) == 1 for s in smis)


@needs_rdkit
def test_parser_formats(tmp_path):
    from app import chem_io

    mol, name, _ = chem_io.parse_ligand("x.smi", b"CCO ethanol\n")
    assert name == "ethanol" and mol.GetNumAtoms() == 3
    with pytest.raises(ValueError):
        chem_io.parse_ligand("x.xyz", b"junk")


@needs_rdkit
def test_properties_and_svg():
    from rdkit import Chem
    from app import chem_io, properties

    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    p = properties.physchem(m)
    assert abs(p["mw"] - 180.16) < 0.1
    assert properties.drug_likeness(p, m)["Lipinski (<=1 violation)"]["passed"]
    assert "<svg" in chem_io.mol_svg(m, [0, 1])


# ---------------------------------------------------------------- no RDKit needed
def test_mdp_100ns_step_count():
    from app import gromacs

    md = gromacs.mdp_files(100, 310.15)["md.mdp"]
    assert "nsteps = 50000000" in md and "ref_t = 310.15" in md


def test_gromacs_static_files_zip():
    from app import gromacs

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        gromacs.write_static_files(z, "H01", "CCO", 0, 100, 310.15)
    names = set(zipfile.ZipFile(io.BytesIO(buf.getvalue())).namelist())
    assert {"md.mdp", "em.mdp", "nvt.mdp", "npt.mdp", "ions.mdp", "run_md.sh", "build_complex.py",
            "parameterise_ligand.sh", "pdbqt_to_sdf.py", "analyze_md.py", "README.txt"} <= names


def test_build_complex_helper(tmp_path):
    """The merge helper must replace ligand coordinates and patch the topology."""
    from app import gromacs

    (tmp_path / "build_complex.py").write_text(gromacs.BUILD_TOP_PY)
    (tmp_path / "prot.gro").write_text(
        "p\n2\n    1ALA      N    1   0.100   0.100   0.100\n    1ALA     CA    2   0.200   0.100   0.100\n   3.0 3.0 3.0\n")
    (tmp_path / "lig.gro").write_text(
        "l\n2\n    1LIG      C    1   9.000   9.000   9.000\n    1LIG      O    2   9.100   9.000   9.000\n   3.0 3.0 3.0\n")
    sdf = ("x\n  t\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
           "   10.0000   20.0000   30.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
           "   11.0000   20.0000   30.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
           "  1  2  1  0\nM  END\n$$$$\n")
    (tmp_path / "pose.sdf").write_text(sdf)
    (tmp_path / "topol.top").write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n[ molecules ]\nProtein 1\n')
    r = subprocess.run([sys.executable, "build_complex.py", "prot.gro", "lig.gro", "pose.sdf", "topol.top", "LIG.itp"],
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    gro = (tmp_path / "complex.gro").read_text()
    assert "   1.000   2.000   3.000" in gro            # 10/20/30 angstrom -> nm
    top = (tmp_path / "topol.top").read_text()
    assert '#include "LIG.itp"' in top and "LIG" in top.split("[ molecules ]")[1]


def test_build_complex_rejects_atom_mismatch(tmp_path):
    from app import gromacs

    (tmp_path / "build_complex.py").write_text(gromacs.BUILD_TOP_PY)
    (tmp_path / "prot.gro").write_text("p\n1\n    1ALA      N    1   0.100   0.100   0.100\n   3.0 3.0 3.0\n")
    (tmp_path / "lig.gro").write_text("l\n2\n    1LIG      C    1   1.0 1.0 1.0\n    1LIG      O    2   1.1 1.0 1.0\n   3.0 3.0 3.0\n")
    (tmp_path / "pose.sdf").write_text("x\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                                      "    1.0000    1.0000    1.0000 C   0  0\nM  END\n$$$$\n")
    (tmp_path / "topol.top").write_text('#include "x.ff/forcefield.itp"\n')
    r = subprocess.run([sys.executable, "build_complex.py", "prot.gro", "lig.gro", "pose.sdf", "topol.top", "LIG.itp"],
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode != 0 and "atoms" in (r.stdout + r.stderr)


def _fake_result():
    props = dict(formula="C1", mw=100.0, logp=1.0, hbd=1, hba=2, tpsa=30.0, rotb=2, fsp3=0.3, qed=0.7)
    rules = {"Lipinski (<=1 violation)": dict(passed=True, detail="0")}
    return dict(
        job_id="abc", created="2026-10-30 00:00 UTC", params=dict(mode="both"),
        software=dict(linkwright="0.1.0", rdkit="test"), warnings=[],
        parents={k: dict(name=k, smiles="CCO", props=props, xref=dict(verdict="v", queried="t")) for k in "AB"},
        hybrids=[dict(rank=1, id="H01", rule_name="Amide coupling", smiles="CCO", inchikey="K", props=props,
                      explanation="why", details=dict(reference="Ref"), rules=rules,
                      alerts=dict(PAINS=[], Brenk=["x"]),
                      xref=dict(pubchem=dict(verdict="pv", queried="t"), chembl=dict(verdict="cv", queried="t")),
                      admet_ai=None)],
        gromacs=dict(selected=["H01"], ns=100, temp=310.15),
    )


def test_report_docx_and_markdown(tmp_path):
    from app import report

    res = _fake_result()
    md = report.build_markdown(res)
    assert "Amide coupling" in md and "Limitations" in md
    path = tmp_path / "r.docx"
    report.build_docx(res, str(path))
    assert path.stat().st_size > 5000


def test_admet_csv_matches_hybrid_id_and_smiles():
    from app.admet_csv import import_csv
    hybrids = [
        {"id": "H01", "smiles": "CCO"},
        {"id": "H02", "smiles": "CCN"},
    ]
    csv_data = b"hybrid_id,SMILES,hERG,AMES\nH01,CCO,0.12,0.8\n,CCN,0.3,0.4\nUNKNOWN,CCC,0.1,0.2\n"
    out = import_csv(csv_data, hybrids, "Test export")
    assert out["matched_count"] == 2
    assert out["matched"]["H01"]["herg"] == 0.12
    assert out["matched"]["H02"]["ames"] == 0.4
    assert out["unmatched_count"] == 1
    assert out["source"] == "Test export"


def test_admet_csv_requires_identity_column():
    from app.admet_csv import import_csv
    with pytest.raises(ValueError):
        import_csv(b"endpoint,value\\nhERG,0.2\\n", [{"id": "H01", "smiles": "CCO"}])


def test_md_analysis_reads_xvg_and_writes_svg(tmp_path):
    from app.md_analysis import read_xvg, svg_plot
    xvg = tmp_path / "rmsd.xvg"
    xvg.write_text('@ title "RMSD"\n# comment\n0.0 0.10\n1.0 0.20\n2.0 0.15\n')
    points = read_xvg(xvg)
    assert points == [(0.0, 0.1), (1.0, 0.2), (2.0, 0.15)]
    out = tmp_path / "rmsd.svg"
    svg_plot(points, "RMSD", "nm", out)
    assert "<svg" in out.read_text()


def test_surechembl_response_helpers():
    from app.surechembl import _find_key, _extract_records
    payload = {"data": {"searchHash": "abc123", "results": [{"id": "SCHEMBL1"}]}}
    assert _find_key(payload, {"searchhash"}) == "abc123"
    assert _extract_records(payload) == [{"id": "SCHEMBL1"}]


def test_api_key_gate(monkeypatch):
    from fastapi import HTTPException
    from app import main
    monkeypatch.setenv("LINKWRIGHT_API_KEY", "secret")
    with pytest.raises(HTTPException) as err:
        main._require_api_key("wrong")
    assert err.value.status_code == 401
    main._require_api_key("secret")


def test_docking_route_is_exposed_as_multipart_api():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from app.main import app
    schema = app.openapi()
    assert "/api/docking" in schema["paths"]
    assert "multipart/form-data" in schema["paths"]["/api/docking"]["post"]["requestBody"]["content"]
    assert "/docking" in schema["paths"]


def test_pdbqt_to_sdf_helper_refuses_atom_count_mismatch(tmp_path):
    from app import gromacs
    script=tmp_path/"pdbqt_to_sdf.py"
    script.write_text(gromacs.PDBQT_TO_SDF_PY)
    sdf=tmp_path/"ref.sdf"
    sdf.write_text("ethanol\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
                   "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                   "    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
                   "  1  2  1  0\nM  END\n$$$$\n")
    pdbqt=tmp_path/"pose.pdbqt"
    pdbqt.write_text("ATOM      1 C1   LIG A   1       1.000   2.000   3.000  1.00  0.00    0.000 C\\n")
    out=tmp_path/"pose.sdf"
    result=subprocess.run([sys.executable,str(script),str(sdf),str(pdbqt),str(out)],capture_output=True,text=True)
    assert result.returncode != 0
    assert "Atom-count mismatch" in (result.stdout+result.stderr)

@needs_rdkit
def test_ligand_preparation_records_identity_and_provenance():
    from app.structure_qc import prepare_ligand
    result = prepare_ligand("ethanol.smi", None, "CCO", seed=42)
    assert result["canonical_smiles"] == "CCO"
    assert result["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert result["optimization_method"] in {"MMFF94s", "UFF"}
    assert "V2000" in result["prepared_sdf"] or "V3000" in result["prepared_sdf"]
    assert len(result["input_sha256"]) == 64
    assert any("not PDBQT" in w for w in result["warnings"])


def test_protein_pdb_inventory_and_warnings():
    from app.structure_qc import validate_protein_pdb
    pdb = ("ATOM      1  N   ALA A   1      11.104  13.207   9.321  1.00 20.00           N\n"
           "ATOM      2  CA  ALA A   1      12.000  13.000   9.000  1.00 20.00           C\nEND\n").encode()
    result = validate_protein_pdb(pdb)
    assert result["atom_records"] == 2
    assert result["unique_residues"] == 1
    assert result["chains"] == ["A"]
    assert len(result["input_sha256"]) == 64


def test_v04_routes_are_exposed():
    from app.main import app
    paths = app.openapi()["paths"]
    assert "/api/ligand/prepare" in paths
    assert "/api/protein/validate" in paths
    assert "/structure-tools" in paths


def test_pdbqt_preflight_requires_ligand_root_and_valid_coordinates(tmp_path):
    from app.vina import validate_pdbqt
    atom = "ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00     0.000 C\n"
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text("ROOT\n" + atom + "ENDROOT\nTORSDOF 0\n")
    result = validate_pdbqt(str(ligand), "ligand")
    assert result["atom_records"] == 1 and len(result["sha256"]) == 64
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(atom)
    assert validate_pdbqt(str(receptor), "receptor")["atom_records"] == 1
    ligand.write_text(atom)
    with pytest.raises(ValueError, match="ROOT"):
        validate_pdbqt(str(ligand), "ligand")
