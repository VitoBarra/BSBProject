"""Quantification configuration and workflow entry points."""

from .config import QuantificationConfig, available_memory_bytes
from .salmon_runner import build_salmon_index, run_salmon
