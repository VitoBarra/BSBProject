from __future__ import annotations

from pathlib import Path

from . import PROJECT_ROOT


def portable_path(value: str | Path) -> str:
    """Return project paths as relative strings and external paths as absolute strings."""
    resolved = Path(value).resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(resolved)
