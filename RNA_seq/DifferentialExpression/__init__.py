"""Differential-expression configuration and workflow entry points."""

from .config import (
    DESEQ2_NUMERIC_RESULT_COLUMNS,
    DESEQ2_RESULT_COLUMNS,
    DifferentialExpressionConfig,
)
from .prepare_deseq2_inputs import prepare_deseq2_inputs
from .deseq2_runner import generate_deseq2_plots, run_deseq2_analysis
