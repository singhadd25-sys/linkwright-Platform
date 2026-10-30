# Changelog

## 0.4.2
- Added a streaming importer for the RetroRules compressed CSV template dataset.
- Added persistent SQLite indexing with source fingerprinting, SHA-256 provenance, template counts and per-source counts.
- Added API endpoints for dataset status, indexing, template search and single-template retrieval.
- Added optional RDKit reaction-SMARTS parsing/validation for retrieved templates.
- Explicitly separated database indexing/search from forward reaction execution; no claim is made that all RetroRules entries are direct two-parent synthesis rules.
- Added automated tests using a small synthetic compressed dataset.

## 0.4.0
- Added reproducible RDKit ligand preparation: hydrogen addition, ETKDGv3 conformer generation, MMFF94s/UFF minimization, canonical SMILES, InChIKey, formula, input/output SHA-256 and recorded seed/version.
- Added conservative protein PDB coordinate inventory endpoint for chains, residues, alternate locations, waters, hetero atoms and malformed coordinate records.
- Added a browser-based Structure Quality Tools page with prepared SDF download.
- Added basic PDBQT syntax preflight and SHA-256 provenance for receptor, ligand and docked pose outputs.
- Added explicit warnings that minimized SDF is not PDBQT and that protein inventory is not protein repair or completeness validation.
- Fixed a job-creation timestamp variable bug in the hybrid-generation API.


## 0.3.0
- Added a dedicated browser interface and background batch endpoint for AutoDock Vina PDBQT docking.
- Added SureChEMBL 2.0 identity/similarity structure-search adapter with pending/unavailable states.
- Added GROMACS XVG summarisation to CSV/JSON and SVG plots.
- Added strict PDBQT-to-reference-SDF coordinate transfer helper for compatible full-atom atom-order-preserving poses.
- Removed permissive `-maxwarn 1` flags from generated GROMACS workflow.
- Switched to the portable `System` thermostat group; documented manual ligand restraint review.
- Added optional shared API-key protection, upload size limits and shared per-IP hourly rate limiting.
- Expanded reports, documentation and automated tests.

## Validation status
- Python syntax compilation: passed.
- Automated test suite: 17 passed.
- FastAPI health endpoint and docking page smoke checks: passed.
- Live SureChEMBL request, actual AutoDock Vina docking and real GROMACS simulation: not executed in this build environment.

## 0.4.0 validation status
- RDKit ligand preparation smoke test: passed (ethanol, MMFF94s, deterministic seed).
- Protein PDB inventory smoke test: passed (synthetic two-atom fixture).
- Full automated suite: 21 passed.
- FastAPI health, structure-tools page, ligand preparation API and protein validation API smoke tests: passed.
- No biological benchmark, redocking, Vina run or MD run is claimed by this release.
