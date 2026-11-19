"""PCR Nexus job adapter for the original qpcr_auto analysis engine."""
from __future__ import annotations

import json
import os
import shutil
import zipfile
from pathlib import Path
from typing import Any

from . import qpcr_auto

SUPPORTED_SUFFIXES = {".csv", ".tsv", ".txt", ".xlsx", ".xlsm"}


def analyse_pcr_file(source_path: str, job_dir: str, plate_path: str | None = None,
                     reference_gene: str | None = None, control_group: str | None = None,
                     dpcr_partition_volume_ul: float = 0.00085) -> dict[str, Any]:
    """Run the uploaded, unchanged qpcr_auto engine and package its generated outputs."""
    suffix = Path(source_path).suffix.lower()
    if suffix not in SUPPORTED_SUFFIXES:
        raise ValueError("Unsupported file type. Upload CSV, TSV, TXT, XLSX or XLSM.")
    out_dir = os.path.join(job_dir, "analysis")
    detection, analysis, paths = qpcr_auto.run_pipeline(
        source_path, plate=plate_path, ref=reference_gene, control=control_group,
        dpcr_vol=dpcr_partition_volume_ul, outdir=out_dir
    )
    if analysis is None or paths is None:
        raise ValueError("PCR file format was not recognised. Check that it contains numeric Ct/Cq, amplification, melt/HRM or digital-PCR data.")
    # Copy stable download names to the job root; keep the complete engine output in analysis/.
    outputs = {
        "report": os.path.join(job_dir, "pcr_report.html"),
        "tables": os.path.join(job_dir, "pcr_tables.xlsx"),
        "calculations": os.path.join(job_dir, "pcr_calculations.txt"),
    }
    shutil.copy2(paths["report"], outputs["report"])
    shutil.copy2(paths["xlsx"], outputs["tables"])
    calc_src = os.path.join(out_dir, "calculations.txt")
    if os.path.isfile(calc_src):
        shutil.copy2(calc_src, outputs["calculations"])
    package = os.path.join(job_dir, "pcr_package.zip")
    with zipfile.ZipFile(package, "w", zipfile.ZIP_DEFLATED) as archive:
        for root, _, files in os.walk(out_dir):
            for name in files:
                full = os.path.join(root, name)
                archive.write(full, os.path.relpath(full, out_dir))
    summary = {
        "pcr_type": detection["pcr_type"],
        "instrument": detection["instrument"],
        "confidence": detection["confidence"],
        "evidence": detection["evidence"],
        "qc": analysis.qc,
        "tables": list(analysis.tables.keys()),
        "figures": [f["key"] for f in analysis.figs],
        "well_count": analysis.meta.get("n_wells"),
        "target_count": len(analysis.meta.get("targets", [])),
        "warnings": analysis.qc,
        "files": {
            "report": "/api/pcr/jobs/{job_id}/files/pcr_report.html",
            "tables": "/api/pcr/jobs/{job_id}/files/pcr_tables.xlsx",
            "calculations": "/api/pcr/jobs/{job_id}/files/pcr_calculations.txt",
            "package": "/api/pcr/jobs/{job_id}/files/pcr_package.zip",
        },
    }
    with open(os.path.join(job_dir, "pcr_results.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary
