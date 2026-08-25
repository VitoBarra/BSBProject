from __future__ import annotations

import os
import logging
import shlex
import subprocess
from pathlib import Path
from typing import Mapping, Sequence


class ExternalToolRunner:
    """Run a command-line app from the Docker runtime environment."""

    def __init__(self, executable: str, display_name: str, logger: logging.Logger) -> None:
        self.executable = executable
        self.display_name = display_name
        self.logger = logger

    def path_arg(self, path: Path) -> str:
        return str(path.resolve())

    def run(
        self,
        args: Sequence[str],
        *,
        env: Mapping[str, str] | None = None,
        missing_message: str | None = None,
    ) -> None:
        command = [self.executable, *args]
        error_context = missing_message or (
            f"Unable to start {self.display_name}. "
            f"Verify that '{self.executable}' is installed and available in PATH."
        )

        self.logger.info("Running command: %s", " ".join(shlex.quote(part) for part in command))
        try:
            subprocess.run(command, check=True, env=self._merged_env(env))
        except FileNotFoundError as exc:
            raise FileNotFoundError(error_context) from exc

    def _merged_env(self, env: Mapping[str, str] | None) -> Mapping[str, str]:
        merged = os.environ.copy()
        if env is not None:
            merged.update(env)
        return merged
