FROM bioconductor/bioconductor_docker:RELEASE_3_21

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 \
    && rm -rf /var/lib/apt/lists/* \
    && Rscript -e "BiocManager::install(c('DESeq2', 'SummarizedExperiment', 'clusterProfiler', 'org.Hs.eg.db', 'AnnotationDbi'), ask = FALSE, update = FALSE)" \
    && Rscript -e "install.packages('ggplot2', repos = 'https://cloud.r-project.org')"

WORKDIR /project
