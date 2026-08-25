from __future__ import annotations

import logging
from pathlib import Path

from ExternalTools import ExternalToolRunner

LOGGER = logging.getLogger("multiqc")


def run_multiqc(
    input_dirs: list[Path],
    output_dir: Path,
    report_name: str,
) -> None:
    """Aggregate reports from ``input_dirs`` into ``output_dir``.

    The generated HTML report is named ``report_name``.
    """

    missing = [path for path in input_dirs if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing MultiQC input directories: {', '.join(str(path) for path in missing)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Input directories: %s", ", ".join(str(path) for path in input_dirs))
    LOGGER.info("Output directory: %s", output_dir)

    runner = ExternalToolRunner(executable="multiqc", display_name="MultiQC", logger=LOGGER)
    runner.run(
        [
            "--force",  # Overwrite an existing report with the same name.
            "--dirs",  # Prefix sample names with their source report directory.
            # Keep only the final directory name in that prefix.
            "--dirs-depth",
            "1",
            # Write the HTML report and MultiQC data directory here.
            "--outdir",
            runner.path_arg(output_dir),
            # Set the output HTML report filename.
            "--filename",
            report_name,
            # Search each supplied FastQC or FastP report directory.
            *(runner.path_arg(path) for path in input_dirs),
        ],
        missing_message="MultiQC not found in PATH",
    )
    LOGGER.info("Done")
