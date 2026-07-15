PROJECT_DIR := $(CURDIR)/Project
IMAGE := bsbproject-r-analysis:latest
# Linux bind mounts benefit from matching the host UID/GID. Native Windows
# does not provide `id`, and Docker Desktop does not require this mapping.
ifeq ($(OS),Windows_NT)
DOCKER_USER :=
else
HOST_UID := $(shell id -u)
HOST_GID := $(shell id -g)
DOCKER_USER := --user "$(HOST_UID):$(HOST_GID)"
endif

DOCKER_RUN := docker run --rm $(DOCKER_USER) -e HOME=/tmp -v "$(PROJECT_DIR):/project" -w /project $(IMAGE)
DOCKER_SHELL := docker run --rm -it $(DOCKER_USER) -e HOME=/tmp -v "$(PROJECT_DIR):/project" -w /project $(IMAGE)
PYTHON := python3

.PHONY: help docker-build shell run prepare-dea dea plots enrichment analysis multiqc-raw multiqc-clean multiqc-compare

help:
	@echo "Targets:"
	@echo "  docker-build   Build the Docker image with Python, R, DESeq2, and enrichment packages"
	@echo "  run ARGS=...   Run Project/main.py inside the Docker container"
	@echo "  prepare-dea    Aggregate Salmon transcript estimates into gene-level DEA inputs"
	@echo "  dea            Run paired DESeq2 differential expression analysis"
	@echo "  plots          Regenerate DEA plots from existing result tables"
	@echo "  multiqc-raw    Generate a MultiQC report from raw FastQC results"
	@echo "  multiqc-clean  Generate a MultiQC report from trimmed FastQC and fastp results"
	@echo "  enrichment     Run GO enrichment from DESeq2 results"
	@echo "  analysis       Run DEA followed by enrichment"
	@echo "  shell          Open a shell inside the analysis container"

docker-build:
	docker build -t $(IMAGE) -f Project/Dockerfile Project

run:
	$(DOCKER_RUN) $(PYTHON) main.py $(ARGS)

prepare-dea:
	$(MAKE) run ARGS="--prepare-dea-inputs"

dea:
	$(MAKE) run ARGS="--run-dea"

plots:
	$(DOCKER_RUN) Rscript DifferentialExpression/scripts/plot_de_results.R data/GSE103001/de/results/deseq2_all_genes.csv data/GSE103001/de/results/normalized_counts.csv data/GSE103001/de/sample_table.tsv data/GSE103001/de/results

multiqc-raw:
	$(DOCKER_RUN) multiqc --force --dirs --dirs-depth 1 --outdir /project/data/GSE103001/qc/multiqc_raw --filename multiqc_raw.html /project/data/GSE103001/qc/fastqc

multiqc-clean:
	$(DOCKER_RUN) multiqc --force --dirs --dirs-depth 1 --outdir /project/data/GSE103001/qc/multiqc_clean --filename multiqc_clean.html /project/data/GSE103001/qc/fastqc_trimmed /project/data/GSE103001/qc/fastp

enrichment:
	$(MAKE) run ARGS="--run-enrichment"

analysis:
	$(MAKE) run ARGS="--prepare-dea-inputs --run-dea --run-enrichment"

shell:
	$(DOCKER_SHELL) bash
