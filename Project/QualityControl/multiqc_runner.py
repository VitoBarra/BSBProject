from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from Log.log_util import log

from . import QualityControlConfig
from .fastqc_runner import _windows_to_wsl_path

LOG_PREFIX = "multiqc"


def _log(message: str) -> None:
    log(message, LOG_PREFIX)


def run_multiqc(
    config: QualityControlConfig,
    executable: str = "multiqc",
    report_name: str = "multiqc_report",
) -> Path:
    qc_input_dir = config.resolved_qc_root()
    out_dir = config.resolved_multiqc_report_out()
    if not qc_input_dir.exists():
        raise FileNotFoundError(f"Missing QC directory: {qc_input_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    wsl_executable = shutil.which("wsl") or r"C:\Windows\System32\wsl.exe"
    command = [
        wsl_executable,
        "bash",
        "-lc",
        f"command -v {executable} >/dev/null 2>&1 || {{ echo 'MultiQC not found in WSL PATH' >&2; exit 127; }}; "
        f"PYTHONNOUSERSITE=1 {executable} --force --outdir '{_windows_to_wsl_path(out_dir)}' --filename '{report_name}.html' '{_windows_to_wsl_path(qc_input_dir)}'",
    ]

    _log(f"QC input directory: {qc_input_dir}")
    _log(f"MultiQC output directory: {out_dir}")
    _log(f"Running command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "Unable to start the WSL MultiQC command. Verify that WSL is installed and reachable from Python."
        ) from exc
    _log("Done")
    return out_dir / f"{report_name}.html"
