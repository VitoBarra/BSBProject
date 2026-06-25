PROJECT_DIR := $(CURDIR)/Project
IMAGE := bsbproject-r-analysis:latest
HOST_UID := $(shell id -u)
HOST_GID := $(shell id -g)
DOCKER_RUN := docker run --rm --user "$(HOST_UID):$(HOST_GID)" -e HOME=/tmp -v "$(PROJECT_DIR):/project" -w /project $(IMAGE)
PYTHON := python3

.PHONY: help docker-build shell run prepare-dea dea plots enrichment analysis

help:
	@echo "Targets:"
	@echo "  docker-build   Build the Docker image with Python, R, DESeq2, and enrichment packages"
	@echo "  run ARGS=...   Run Project/main.py inside the Docker container"
	@echo "  prepare-dea    Aggregate Salmon transcript estimates into gene-level DEA inputs"
	@echo "  dea            Run paired DESeq2 differential expression analysis"
	@echo "  plots          Regenerate DEA plots from existing result tables"
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
	Rscript Project/DifferentialExpression/scripts/plot_de_results.R \
		Project/data/GSE103001/de/results/deseq2_all_genes.csv \
		Project/data/GSE103001/de/results/normalized_counts.csv \
		Project/data/GSE103001/de/sample_table.tsv \
		Project/data/GSE103001/de/results

enrichment:
	$(MAKE) run ARGS="--run-enrichment"

analysis:
	$(MAKE) run ARGS="--prepare-dea-inputs --run-dea --run-enrichment"

shell:
	$(DOCKER_RUN) bash
