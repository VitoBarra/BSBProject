from __future__ import annotations

from pathlib import Path

from ExternalTools import ExternalToolRunner
from Log.log_util import log

from . import QualityControlConfig

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
    _log(f"QC input directory: {qc_input_dir}")
    _log(f"MultiQC output directory: {out_dir}")
    runner = ExternalToolRunner(executable=executable, display_name="MultiQC", log=_log)
    runner.run(
        [
            "--force",
            "--dirs",
            "--dirs-depth",
            "1",
            "--outdir",
            runner.path_arg(out_dir),
            "--filename",
            f"{report_name}.html",
            runner.path_arg(qc_input_dir),
        ],
        env={"PYTHONNOUSERSITE": "1"},
        missing_message="MultiQC not found in PATH",
    )
    _log("Done")
    return out_dir / f"{report_name}.html"
