"""GROMACS input preparation for protein-ligand MD.

Linkwright does not run molecular dynamics itself. It writes a reproducible package per hybrid:
  ligand.sdf (3D, explicit H) + parameterisation script (AmberTools antechamber, GAFF2, AM1-BCC,
  then acpype) + .mdp files (EM, NVT, NPT, production) + topology/merge helpers + run script.

The user supplies protein.pdb and a docked pose of the ligand (same atom order as ligand.sdf,
hydrogens included). The merge helper stops with an error if atom counts or elements disagree.
"""
from __future__ import annotations

import io
import os
import zipfile

FORCE_FIELD = "amber99sb-ildn"
WATER = "tip3p"


def mdp_files(ns: int = 100, temp: float = 310.15, dt: float = 0.002) -> dict[str, str]:
    prod_steps = int(round(ns * 1000 / dt))
    common_nb = (
        "cutoff-scheme = Verlet\nnstlist = 20\ncoulombtype = PME\nrcoulomb = 1.0\n"
        "vdwtype = Cut-off\nrvdw = 1.0\nDispCorr = EnerPres\npbc = xyz\n"
    )
    groups = "tc-grps = System\ntau_t = 0.1\n" f"ref_t = {temp}\n"
    return {
        "ions.mdp": (
            "integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 50000\n" + common_nb +
            "constraints = none\n"
        ),
        "em.mdp": (
            "integrator = steep\nemtol = 1000.0\nemstep = 0.01\nnsteps = 50000\nnstenergy = 500\n"
            + common_nb + "constraints = h-bonds\nconstraint_algorithm = lincs\n"
        ),
        "nvt.mdp": (
            "define = -DPOSRES\nintegrator = md\n"
            f"nsteps = {int(round(100 / dt))}\ndt = {dt}\n"
            "nstxout-compressed = 5000\nnstenergy = 1000\nnstlog = 1000\n"
            "continuation = no\nconstraints = h-bonds\nconstraint_algorithm = lincs\n"
            "tcoupl = V-rescale\n" + groups + "pcoupl = no\ngen_vel = yes\n"
            f"gen_temp = {temp}\ngen_seed = -1\n" + common_nb
        ),
        "npt.mdp": (
            "define = -DPOSRES\nintegrator = md\n"
            f"nsteps = {int(round(100 / dt))}\ndt = {dt}\n"
            "nstxout-compressed = 5000\nnstenergy = 1000\nnstlog = 1000\n"
            "continuation = yes\nconstraints = h-bonds\nconstraint_algorithm = lincs\n"
            "tcoupl = V-rescale\n" + groups +
            "pcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\n"
            "compressibility = 4.5e-5\nrefcoord_scaling = com\ngen_vel = no\n" + common_nb
        ),
        "md.mdp": (
            "integrator = md\n"
            f"nsteps = {prod_steps}   ; {ns} ns at dt = {dt} ps\ndt = {dt}\n"
            "nstxout-compressed = 5000   ; one frame per 10 ps\nnstenergy = 5000\nnstlog = 5000\n"
            "continuation = yes\nconstraints = h-bonds\nconstraint_algorithm = lincs\n"
            "tcoupl = V-rescale\n" + groups +
            "pcoupl = C-rescale\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\n"
            "compressibility = 4.5e-5\ngen_vel = no\n" + common_nb
        ),
    }


def ligand_sdf(mol, name: str) -> str:
    """3D conformer with explicit hydrogens (ETKDGv3 + MMFF94 clean-up)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(m, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(m, params) != 0:
            raise ValueError("3D embedding failed for this hybrid")
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    except Exception:  # noqa: BLE001
        pass
    m.SetProp("_Name", name)
    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    w.write(m)
    w.close()
    return buf.getvalue()


PARAMETERISE_SH = """#!/usr/bin/env bash
# Ligand parameterisation: GAFF2 + AM1-BCC charges (AmberTools), converted to GROMACS by acpype.
# Requires: antechamber, parmchk2, acpype on PATH.
set -euo pipefail
CHARGE={charge}
antechamber -i ligand.sdf -fi mdl -o LIG.mol2 -fo mol2 -c bcc -nc "$CHARGE" -at gaff2 -rn LIG -s 2
acpype -i LIG.mol2 -b LIG -c user -a gaff2 -n "$CHARGE"
echo "Done. Check LIG.acpype/LIG_GMX.itp and LIG_GMX.gro."
"""

PDBQT_TO_SDF_PY = '#!/usr/bin/env python3\n"""Transfer Vina PDBQT coordinates onto a reference ligand SDF in strict atom order."""\nimport sys\nfrom rdkit import Chem\n\ndef pdbqt_atoms(path):\n    out=[]\n    for line in open(path,encoding="utf-8",errors="replace"):\n        if line.startswith(("ATOM  ","HETATM")):\n            name=line[12:16].strip()\n            x,y,z=map(float,(line[30:38],line[38:46],line[46:54]))\n            typ=line[77:79].strip().upper() if len(line)>=79 else ""\n            if typ in {"C","A"}: el="C"\n            elif typ.startswith("N"): el="N"\n            elif typ.startswith("O"): el="O"\n            elif typ.startswith("S"): el="S"\n            elif typ=="P": el="P"\n            elif typ.startswith("CL"): el="Cl"\n            elif typ.startswith("BR"): el="Br"\n            elif typ in {"F","I","H","HD"}: el="H" if typ=="HD" else typ\n            else:\n                raw="".join(c for c in name if c.isalpha())\n                el=raw[:1].upper()+raw[1:2].lower() if raw else ""\n            out.append((el,(x,y,z),name))\n    return out\n\ndef main():\n    if len(sys.argv)!=4:\n        raise SystemExit("Usage: python pdbqt_to_sdf.py reference_ligand.sdf vina_pose.pdbqt docked_pose.sdf")\n    ref=Chem.SDMolSupplier(sys.argv[1],removeHs=False)[0]\n    if ref is None or ref.GetNumConformers()==0:\n        raise SystemExit("Reference SDF is invalid or has no 3D coordinates.")\n    atoms=pdbqt_atoms(sys.argv[2])\n    if len(atoms)!=ref.GetNumAtoms():\n        raise SystemExit(f"Atom-count mismatch: reference SDF has {ref.GetNumAtoms()} atoms; PDBQT has {len(atoms)}. "\n                         "The PDBQT must retain all atoms, including hydrogens, in matching order.")\n    conf=ref.GetConformer()\n    for i,(el,(x,y,z),name) in enumerate(atoms):\n        expected=ref.GetAtomWithIdx(i).GetSymbol()\n        if el.capitalize()!=expected:\n            raise SystemExit(f"Element/order mismatch at atom {i+1}: SDF={expected}, PDBQT={el} ({name}).")\n        conf.SetAtomPosition(i,(x,y,z))\n    w=Chem.SDWriter(sys.argv[3]); w.write(ref); w.close()\n    print("Wrote docked_pose.sdf with reference bond orders and Vina coordinates.")\nif __name__=="__main__": main()\n'

BUILD_TOP_PY = '''#!/usr/bin/env python3
"""Merge protein + docked ligand and patch topol.top. Usage:
    python build_complex.py protein_processed.gro LIG.acpype/LIG_GMX.gro docked_pose.sdf topol.top LIG.acpype/LIG_GMX.itp
Writes complex.gro and edits topol.top in place (backup: topol.top.bak)."""
import re, shutil, sys

prot_gro, lig_gro, pose_sdf, top, itp = sys.argv[1:6]


def read_sdf_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[3][:3])
    out = []
    for l in lines[4:4 + n]:
        out.append((float(l[0:10]), float(l[10:20]), float(l[20:30]), l[31:34].strip()))
    return out


def gro_atoms(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    return L[2:2 + n], L[2 + n]


def element(name):
    name = re.sub(r"[^A-Za-z]", "", name)
    if name[:2] in ("Cl", "CL", "Br", "BR"):
        return name[:2].capitalize()
    return name[:1].upper()


pose = read_sdf_xyz(pose_sdf)
lig_lines, _ = gro_atoms(lig_gro)
if len(pose) != len(lig_lines):
    sys.exit(f"ERROR: docked pose has {len(pose)} atoms but the ligand topology has {len(lig_lines)}. "
             "The pose must include hydrogens in the same atom order as ligand.sdf.")
new_lig = []
for (x, y, z, el), line in zip(pose, lig_lines):
    nm = line[10:15].strip()
    if element(nm) != el.capitalize():
        sys.exit(f"ERROR: element mismatch ({nm} vs {el}). Atom order differs from ligand.sdf.")
    new_lig.append(line[:20] + f"{x/10:8.3f}{y/10:8.3f}{z/10:8.3f}")
p_lines, box = gro_atoms(prot_gro)
with open("complex.gro", "w") as f:
    f.write("Protein-ligand complex\\n")
    f.write(f"{len(p_lines) + len(new_lig)}\\n")
    f.write("\\n".join(p_lines + new_lig) + "\\n" + box + "\\n")

shutil.copy(top, top + ".bak")
t = open(top).read()
inc = f'#include "{itp}"\\n'
m = re.search(r'#include ".*forcefield\\.itp"\\n', t)
if not m:
    sys.exit("ERROR: could not find the forcefield include in topol.top")
if inc not in t:
    t = t[:m.end()] + inc + t[m.end():]
if not re.search(r"^LIG\\s+1", t, re.M):
    t = t.rstrip() + "\\nLIG                 1\\n"
open(top, "w").write(t)
print("Wrote complex.gro and updated topol.top")
'''


def run_md_sh(ns: int) -> str:
    return f"""#!/usr/bin/env bash
# Protein-ligand MD workflow ({ns} ns production). Needs GROMACS >= 2021, a GPU build for speed.
# Inputs in this folder: protein.pdb, docked_pose.sdf (hydrogens, same atom order as ligand.sdf)
# Run first: bash parameterise_ligand.sh
set -euo pipefail
GMX=${{GMX:-gmx}}

$GMX pdb2gmx -f protein.pdb -ff {FORCE_FIELD} -water {WATER} -ignh -o protein_processed.gro -p topol.top
python3 build_complex.py protein_processed.gro LIG.acpype/LIG_GMX.gro docked_pose.sdf topol.top LIG.acpype/LIG_GMX.itp

$GMX editconf -f complex.gro -o box.gro -c -d 1.2 -bt dodecahedron
$GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
$GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

$GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -v -deffnm em

# Single System thermostat group avoids non-portable index group numbers.
# Protein position restraints are supplied by pdb2gmx; ligand restraints still require review.
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$GMX mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$GMX mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
$GMX mdrun -deffnm md -nb gpu -pme gpu -bonded gpu

# Basic analysis (select groups when prompted)
printf "Protein\\nSystem\\n" | $GMX trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center
printf "Protein\\nProtein\\n" | $GMX rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg
printf "Protein\\n" | $GMX rmsf -s md.tpr -f md_noPBC.xtc -o rmsf.xvg -res
printf "Protein\\n" | $GMX gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg
echo "Hydrogen-bond analysis needs system-specific donor/acceptor selections; run manually."
python3 analyze_md.py --dir .
"""


README = """Linkwright GROMACS package
==========================
Hybrid: {name}
SMILES: {smiles}
Net charge used for parameterisation: {charge}

Steps
1. bash parameterise_ligand.sh      (AmberTools + acpype required)
2. Place protein.pdb and docked_pose.sdf in this folder.
   The docked pose must contain hydrogens and keep the atom order of ligand.sdf.
3. If starting from a Vina output, run:
   python3 pdbqt_to_sdf.py ligand.sdf vina_pose.pdbqt docked_pose.sdf
   This requires every atom, including hydrogens, in matching order.
4. bash run_md.sh                   (GROMACS 2021 or newer; GPU recommended)

Notes
- Protein force field: {ff} with {water}. Ligand: GAFF2 with AM1-BCC charges.
  If your study uses AMBER ff14SB, install a GROMACS port of it and change the pdb2gmx line.
- The thermostat uses one System group for portability. Assess whether this is suitable for
  your system. Protein restraints are generated by pdb2gmx; ligand restraints require manual
  generation and topology review before production.
- pdbqt_to_sdf.py transfers coordinates only and preserves reference SDF bond orders. It refuses
  atom-count/order mismatches. Many PDBQT workflows merge non-polar hydrogens, so their output
  will not pass this strict conversion unless a validated full-atom workflow is used.
- analyze_md.py creates descriptive CSV/JSON summaries and SVG plots from GROMACS XVG files.
  These outputs do not establish convergence; assess equilibration, replicate runs and uncertainty.
- Hydrogen-bond analysis is left manual because donor/acceptor selections are system-specific.
- Linkwright does not execute MD on the web server; you run the package on your own system.
"""


def build_package(mol, name: str, smiles: str, ns: int = 100, temp: float = 310.15) -> bytes:
    from rdkit import Chem

    charge = Chem.GetFormalCharge(mol)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("ligand.sdf", ligand_sdf(mol, name))
        write_static_files(z, name, smiles, charge, ns, temp)
    return buf.getvalue()


def write_static_files(z: zipfile.ZipFile, name: str, smiles: str, charge: int, ns: int, temp: float):
    for fn, txt in mdp_files(ns, temp).items():
        z.writestr(fn, txt)
    z.writestr("parameterise_ligand.sh", PARAMETERISE_SH.format(charge=charge))
    z.writestr("build_complex.py", BUILD_TOP_PY)
    z.writestr("pdbqt_to_sdf.py", PDBQT_TO_SDF_PY)
    z.writestr("run_md.sh", run_md_sh(ns))
    from pathlib import Path
    z.writestr("analyze_md.py", Path(__file__).with_name("md_analysis.py").read_text(encoding="utf-8"))
    z.writestr("README.txt", README.format(name=name, smiles=smiles, charge=charge, ff=FORCE_FIELD, water=WATER))


def safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in s)[:60] or "hybrid"
