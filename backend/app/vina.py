"""Optional AutoDock Vina CLI docking adapter.

Requires `vina` and prepared receptor/ligand PDBQT files. This module deliberately does not
pretend to prepare proteins or convert RDKit molecules: those steps need validated toolchains
(Meeko/AutoDockTools) and user review.
"""
from __future__ import annotations
import os
import hashlib
import shutil
import subprocess
from pathlib import Path

def validate_pdbqt(path: str, role: str) -> dict:
    """Basic PDBQT syntax preflight; not a chemical correctness guarantee."""
    if role not in {"receptor", "ligand"}:
        raise ValueError("PDBQT role must be receptor or ligand.")
    raw = Path(path).read_bytes()
    if not raw:
        raise ValueError(f"{role} PDBQT is empty.")
    text = raw.decode("utf-8", errors="replace")
    atom_count = 0
    malformed = []
    has_root = any(line.strip() == "ROOT" for line in text.splitlines())
    has_endroot = any(line.strip() == "ENDROOT" for line in text.splitlines())
    for line_no, line in enumerate(text.splitlines(), 1):
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line) < 70:
                malformed.append(line_no)
                continue
            try:
                xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
                charge = float(line[70:76]) if len(line) >= 76 and line[70:76].strip() else 0.0
            except ValueError:
                malformed.append(line_no)
                continue
            atom_type = line.split()[-1] if line.split() else ""
            if not atom_type or not all(abs(v) < 1e5 for v in xyz) or not abs(charge) < 1e4:
                malformed.append(line_no)
                continue
            atom_count += 1
    if atom_count == 0:
        raise ValueError(f"{role} PDBQT contains no valid ATOM/HETATM coordinate records.")
    if malformed:
        raise ValueError(f"{role} PDBQT has malformed atom records at line(s): {malformed[:10]}.")
    if role == "ligand" and not (has_root and has_endroot):
        raise ValueError("Ligand PDBQT must contain ROOT and ENDROOT records.")
    return {"atom_records": atom_count, "sha256": hashlib.sha256(raw).hexdigest(),
            "has_root": has_root, "has_endroot": has_endroot}


def run_vina(receptor: str, ligand: str, output: str, center: tuple[float,float,float],
             size: tuple[float,float,float], exhaustiveness: int = 8, num_modes: int = 9,
             seed: int = 2026, timeout: int = 3600) -> dict:
    exe = shutil.which("vina")
    if not exe:
        raise RuntimeError("AutoDock Vina executable not found on PATH. Install Vina 1.2+.")
    for label, path in (("receptor", receptor), ("ligand", ligand)):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise ValueError(f"{label} PDBQT file is missing or empty.")
    receptor_qc = validate_pdbqt(receptor, "receptor")
    ligand_qc = validate_pdbqt(ligand, "ligand")
    if len(center) != 3 or len(size) != 3 or any(float(v) <= 0 for v in size):
        raise ValueError("Docking box requires three center coordinates and three positive sizes.")
    if not 1 <= exhaustiveness <= 128 or not 1 <= num_modes <= 20:
        raise ValueError("exhaustiveness must be 1–128 and num_modes 1–20.")
    out = str(Path(output).resolve())
    cmd = [exe, "--receptor", str(Path(receptor).resolve()), "--ligand", str(Path(ligand).resolve()),
           "--out", out, "--center_x", str(center[0]), "--center_y", str(center[1]),
           "--center_z", str(center[2]), "--size_x", str(size[0]), "--size_y", str(size[1]),
           "--size_z", str(size[2]), "--exhaustiveness", str(exhaustiveness),
           "--num_modes", str(num_modes), "--seed", str(seed)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(f"Vina exceeded the {timeout}s time limit.") from e
    if proc.returncode != 0 or not os.path.isfile(out) or os.path.getsize(out) == 0:
        raise RuntimeError("Vina failed: " + (proc.stderr or proc.stdout)[-1500:])
    scores = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                scores.append({"mode": int(parts[0]), "affinity_kcal_mol": float(parts[1])})
            except ValueError:
                pass
    try:
        version_proc = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=15, check=False)
        version = (version_proc.stdout or version_proc.stderr).strip()[:200]
    except (OSError, subprocess.TimeoutExpired):
        version = "version unavailable"
    return {"engine": "AutoDock Vina", "engine_version": version, "command": cmd, "scores": scores,
            "receptor_qc": receptor_qc, "ligand_qc": ligand_qc,
            "output_sha256": hashlib.sha256(Path(out).read_bytes()).hexdigest(),
            "output": out, "stdout": proc.stdout[-5000:], "stderr": proc.stderr[-2000:],
            "interpretation": "Docking scores are approximate ranking estimates, not measured binding affinities."}
