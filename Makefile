PROJECT_DIR := $(CURDIR)
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

DOCKER_RUN := docker run --rm $(DOCKER_USER) -e HOME=/tmp -e PYTHONHASHSEED=0 -v "$(PROJECT_DIR):/workspace" -w /workspace $(IMAGE)
DOCKER_SHELL := docker run --rm -it $(DOCKER_USER) -e HOME=/tmp -e PYTHONHASHSEED=0 -v "$(PROJECT_DIR):/workspace" -w /workspace $(IMAGE)
PYTHON := python3

.PHONY: help init-env docker-build shell run download download-data salmon-index qc-raw trim qc-trimmed quantify dea enrichment analysis network-bc-download network-bc-analyze network-bc network-disgenet-download network-separation-download network-separation-analyze network-separation

help: ## Show this help message
	@echo ""
	@echo "================================================================================"
	@echo "  BSBProject — Makefile targets"
	@echo "================================================================================"
	@echo ""
	@echo "  ENVIRONMENT:"
	@grep -E '^(help|init-env|docker-build|run|shell):.*##' Makefile | \
		awk -F ':.*## ' '{printf "    \033[36m%-30s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "  PART 1 — RNA-seq:"
	@grep -E '^(download-data|qc-raw|trim|qc-trimmed|quantify|dea|enrichment|analysis):.*##' Makefile | \
		awk -F ':.*## ' '{printf "    \033[36m%-30s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "  PART 2 — NETWORK MEDICINE:"
	@grep -E '^(network-bc-download|network-bc-analyze|network-bc|network-disgenet-download|network-separation-download|network-separation-analyze|network-separation):.*##' Makefile | \
		awk -F ':.*## ' '{printf "    \033[36m%-30s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "================================================================================"
	@echo ""

##@ Environment

init-env: ## Create .env.local with the required DisGeNET API-key variable
	@if [ -e .env.local ]; then \
		echo ".env.local already exists; leaving it unchanged"; \
	else \
		printf 'DISGENET_API_KEY=\n' > .env.local; \
		chmod 600 .env.local; \
		echo "Created .env.local; add your DisGeNET API key after DISGENET_API_KEY="; \
	fi

docker-build: ## Build the Docker analysis image
	docker build -t $(IMAGE) -f Dockerfile .

run: ## Run main.py in Docker; pass options with ARGS="..."
	$(DOCKER_RUN) $(PYTHON) main.py $(ARGS)

shell: ## Open a shell in the analysis container
	$(DOCKER_SHELL) bash

##@ Part 1 - Data acquisition

download-data: ## Prepare metadata, paired FASTQ files, and decoy-aware-index references
	$(MAKE) run ARGS="--download-RnaSeq-data"

##@ Part 1 - Raw-read quality control

qc-raw: ## Run FastQC and MultiQC on raw reads
	$(MAKE) run ARGS="--run-raw-fastqc --run-raw-multiqc"

##@ Part 1 - Read trimming

trim: ## Clean paired reads with fastp
	$(MAKE) run ARGS="--run-fastp"

##@ Part 1 - Trimmed-read quality control

qc-trimmed: ## Run FastQC and MultiQC on trimmed reads
	$(MAKE) run ARGS="--run-trimmed-fastqc --run-trimmed-multiqc"

##@ Part 1 - Quantification

quantify: ## Build the Salmon index if needed and quantify reads
	$(MAKE) run ARGS="--run-salmon --build-salmon-index"

##@ Part 1 - Differential expression

dea: ## Prepare inputs, run paired DESeq2, and generate DEA plots
	$(MAKE) run ARGS="--run-dea"

##@ Part 1 - Enrichment

enrichment: ## Run GO enrichment when needed and regenerate its outputs
	$(MAKE) run ARGS="--run-enrichment"

##@ Part 1 - Complete workflow

analysis: ## Run Part 1 from data download through GO enrichment
	$(MAKE) download-data
	$(MAKE) qc-raw
	$(MAKE) trim
	$(MAKE) qc-trimmed
	$(MAKE) quantify
	$(MAKE) dea
	$(MAKE) enrichment

##@ Part 2 - Network medicine

network-bc-download: ## Download the filtered physical STRING network for BC genes
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.bc_network download --score-threshold 0.70 --network-type physical

network-bc-analyze: ## Calculate BC centralities, communities, and the selected module
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.bc_network analyze --target-size 10 --seed 42 --null-iterations 500

network-bc: network-bc-download network-bc-analyze ## Download and analyze the BC network

network-disgenet-download: ## Download curated RA and DM associations from DisGeNET
	@test -f .env.local || (echo "Missing .env.local with DISGENET_API_KEY" && exit 1)
	docker run --rm $(DOCKER_USER) --env-file .env.local -e HOME=/tmp -e PYTHONHASHSEED=0 -v "$(PROJECT_DIR):/workspace" -w /workspace $(IMAGE) $(PYTHON) -m NetworkMedicine.disgenet_download

network-separation-download: ## Download the joint BC, RA, and DM STRING network
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.disease_separation download --score-threshold 0.70 --network-type physical --additional-interactors 650

network-separation-analyze: ## Calculate BC-RA and BC-DM network separation
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.disease_separation analyze

network-separation: network-separation-download network-separation-analyze ## Download and analyze the disease-separation network
