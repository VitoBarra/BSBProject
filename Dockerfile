FROM bioconductor/bioconductor_docker:RELEASE_3_21

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        libcurl4-openssl-dev \
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
        biopython==1.87 \
        matplotlib==3.10.8 \
        multiqc==1.35 \
        networkx==3.6.1 \
        numpy==2.4.1 \
        pycurl==7.46.0 \
    && Rscript -e "BiocManager::install(c('DESeq2', 'tximport', 'SummarizedExperiment', 'clusterProfiler', 'org.Hs.eg.db', 'AnnotationDbi'), ask = FALSE, update = FALSE)" \
    && fastp --version 2>&1 | grep -F "0.23.4" \
    && fastqc --version 2>&1 | grep -F "0.12.1" \
    && salmon --version 2>&1 | grep -F "1.10.2" \
    && Rscript -e "expected <- c(DESeq2='1.48.2', tximport='1.36.1', clusterProfiler='4.16.0', org.Hs.eg.db='3.21.0'); actual <- vapply(names(expected), function(package) as.character(packageVersion(package)), character(1)); if (!identical(actual, expected)) stop(paste('Unexpected R package versions:', paste(names(actual), actual, collapse=', ')))"

RUN Rscript -e "install.packages('getopt', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.7.5.tar.gz', repos = NULL, type = 'source')" \
    && Rscript -e "expected <- '1.7.5'; actual <- as.character(packageVersion('optparse')); if (!identical(actual, expected)) stop(paste('Unexpected optparse version:', actual))"

ENV PATH="/opt/bsbproject-python/bin:${PATH}"
ENV MPLCONFIGDIR="/tmp/matplotlib"

WORKDIR /workspace
