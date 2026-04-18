"""Quality-control configuration and workflow entry points."""

from .config import QualityControlConfig
from .fastp_runner import run_fastp
from .fastqc_runner import run_fastqc
from .multiqc_runner import run_multiqc
