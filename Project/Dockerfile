FROM bioconductor/bioconductor_docker:RELEASE_3_21

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        fastp \
        fastqc \
        procps \
        python3 \
        python3-pip \
        python3-venv \
        salmon \
        unzip \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/bsbproject-python \
    && /opt/bsbproject-python/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/bsbproject-python/bin/pip install --no-cache-dir \
        matplotlib~=3.10.8 \
        multiqc \
        networkx~=3.6.1 \
        numpy~=2.4.1 \
    && Rscript -e "BiocManager::install(c('DESeq2', 'tximport', 'SummarizedExperiment', 'clusterProfiler', 'org.Hs.eg.db', 'AnnotationDbi'), ask = FALSE, update = FALSE)" \
    && Rscript -e "install.packages('ggplot2', repos = 'https://cloud.r-project.org')"

ENV PATH="/opt/bsbproject-python/bin:${PATH}"
ENV MPLCONFIGDIR="/tmp/matplotlib"

WORKDIR /project
