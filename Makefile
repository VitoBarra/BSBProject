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

DOCKER_RUN := docker run --rm $(DOCKER_USER) -e HOME=/tmp -e PYTHONHASHSEED=0 -v "$(PROJECT_DIR):/project" -w /project $(IMAGE)
DOCKER_SHELL := docker run --rm -it $(DOCKER_USER) -e HOME=/tmp -e PYTHONHASHSEED=0 -v "$(PROJECT_DIR):/project" -w /project $(IMAGE)
PYTHON := python3

.PHONY: help docker-build shell run prepare-dea dea-plots enrichment analysis multiqc-raw multiqc-clean network-bc-download network-bc-analyze network-bc network-disgenet-download

help:
	@echo "Targets:"
	@echo "  docker-build   Build the Docker image with Python, R, DESeq2, and enrichment packages"
	@echo "  run ARGS=...   Run Project/main.py inside the Docker container"
	@echo "  download-data  Build metadata and download the selected paired FASTQ files"
	@echo "  download-reference  Download the Ensembl GRCh38 transcriptome"
	@echo "  trim-clean     Run raw QC, fastp cleaning, post-cleaning QC, and MultiQC"
	@echo "  prepare-reads  Download, quality-control, trim, and clean the reads"
	@echo "  quantify       Build the Salmon index if needed and quantify all cleaned reads"
	@echo "  prepare-quantification  Download the transcriptome and run Salmon"
	@echo "  prepare-dea    Aggregate Salmon transcript estimates into gene-level DEA inputs"
	@echo "  dea            Run paired DESeq2 and write analysis tables"
	@echo "  dea-plots      Generate DEA plots from existing tables with Python"
	@echo "  multiqc-raw    Generate a MultiQC report from raw FastQC results"
	@echo "  multiqc-clean  Generate a MultiQC report from trimmed FastQC and fastp results"
	@echo "  enrichment     Run GO enrichment and write the raw result table"
	@echo "  enrichment-plots  Generate GO tables, summary, and plot with Python"
	@echo "  analysis       Run the complete Module 1 workflow from download through enrichment"
	@echo "  network-bc-download  Download class-filtered physical STRING network (experiments/databases >= 0.70)"
	@echo "  network-bc-analyze   Calculate BC centralities, community significance, and final module"
	@echo "  network-bc           Download and analyze the BC network"
	@echo "  network-disgenet-download  Download all curated RA and DM associations"
	@echo "  network-separation  Build class-filtered joint STRING network and compare BC with RA and DM"
	@echo "  shell          Open a shell inside the analysis container"

docker-build:
	docker build -t $(IMAGE) -f Project/Dockerfile Project

run:
	$(DOCKER_RUN) $(PYTHON) main.py $(ARGS)

download-data:
	$(MAKE) run ARGS="--build-metadata-table --download-fastq"

download-reference:
	$(MAKE) run ARGS="--download-reference-transcriptome"

trim-clean:
	$(MAKE) run ARGS="--run-raw-fastqc --run-fastp --run-trimmed-fastqc"
	$(MAKE) multiqc-raw
	$(MAKE) multiqc-clean

prepare-reads:
	$(MAKE) download-data
	$(MAKE) trim-clean

quantify:
	$(MAKE) run ARGS="--run-salmon --build-salmon-index"

prepare-quantification:
	$(MAKE) download-reference
	$(MAKE) quantify

prepare-dea:
	$(MAKE) run ARGS="--prepare-dea-inputs"

dea:
	$(MAKE) run ARGS="--run-dea"

dea-plots:
	$(MAKE) run ARGS="--plot-dea-results"

multiqc-raw:
	$(DOCKER_RUN) multiqc --force --dirs --dirs-depth 1 --outdir /project/data/GSE103001/qc/multiqc_raw --filename multiqc_raw.html /project/data/GSE103001/qc/fastqc

multiqc-clean:
	$(DOCKER_RUN) multiqc --force --dirs --dirs-depth 1 --outdir /project/data/GSE103001/qc/multiqc_clean --filename multiqc_clean.html /project/data/GSE103001/qc/fastqc_trimmed /project/data/GSE103001/qc/fastp

enrichment:
	$(MAKE) run ARGS="--run-enrichment"

enrichment-plots:
	$(MAKE) run ARGS="--plot-enrichment-results"

analysis:
	$(MAKE) prepare-reads
	$(MAKE) prepare-quantification
	$(MAKE) run ARGS="--prepare-dea-inputs --run-dea --plot-dea-results --run-enrichment --plot-enrichment-results"

network-bc-download:
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.bc_network download --score-threshold 0.70 --network-type physical

network-bc-analyze:
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.bc_network analyze --target-size 10 --seed 42 --null-iterations 500

network-bc: network-bc-download network-bc-analyze

network-disgenet-download:
	@test -f .env.local || (echo "Missing .env.local with DISGENET_API_KEY" && exit 1)
	docker run --rm $(DOCKER_USER) --env-file .env.local -e HOME=/tmp -e PYTHONHASHSEED=0 -v "$(PROJECT_DIR):/project" -w /project $(IMAGE) $(PYTHON) -m NetworkMedicine.disgenet_download

network-separation-download:
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.disease_separation download --score-threshold 0.70 --network-type physical --additional-interactors 650

network-separation-analyze:
	$(DOCKER_RUN) $(PYTHON) -m NetworkMedicine.disease_separation analyze

network-separation: network-separation-download network-separation-analyze

shell:
	$(DOCKER_SHELL) bash
