# Data Sources and Accessions

This note summarizes the repositories and accession types used to recover the raw RNA-seq data for `GSE103001`.

## Repository Roles

### GEO
GEO (Gene Expression Omnibus) is the repository used to describe the biological dataset and its samples.

Useful here because it provides:
- the series accession `GSE103001`
- the sample accessions `GSM...`
- the sample titles such as `Pat_12-02_normal`
- links from GEO samples to SRA experiments

### SRA
SRA (Sequence Read Archive) is the archive that organizes sequencing metadata and raw run-level objects.

Useful here because it provides the accession chain from experiment to run:
- `SRP...` for the study
- `SRS...` for the sample
- `SRX...` for the experiment
- `SRR...` for the run

### ENA
ENA (European Nucleotide Archive) is the European archive that mirrors the same sequencing records and often exposes direct FASTQ download links in a simpler way.

Useful here because it provides:
- direct `fastq.gz` URLs
- file sizes
- simple API endpoints that are easy to script

## Accession Types

### GEO accessions
- `GSE...`: GEO series
- `GSM...`: GEO sample

Example:
- `GSE103001` = the whole dataset
- `GSM2752350` = one specific sample in the dataset

### SRA accessions
- `SRP...`: study
- `SRS...`: sample
- `SRX...`: experiment
- `SRR...`: run

Example interpretation:
- `SRP116023` = the whole sequencing study
- `SRX3120292` = one sequencing experiment
- `SRR5962198` = one run that produces FASTQ files

## Conceptual Hierarchy

The SRA hierarchy is:

```text
Study -> Sample -> Experiment -> Run
SRP      SRS       SRX           SRR
```

A practical interpretation is:
- **Study (`SRP`)**: the project container
- **Sample (`SRS`)**: the biological specimen registered in SRA/BioSample
- **Experiment (`SRX`)**: how that sample was sequenced
- **Run (`SRR`)**: the concrete sequencing output associated with FASTQ files

## Practical Resolution Chain Used in This Project

To go from the biological sample listed in GEO to the actual FASTQ files, we use this chain:

```text
GSM -> SRX -> SRR -> FASTQ URL
```

Example:

```text
GSM2752350
-> Pat_12-02_normal
-> SRX3120292
-> SRR5962198
-> SRR5962198_1.fastq.gz + SRR5962198_2.fastq.gz
```

## Why Multiple Steps Are Needed

The raw data are not stored in one single flat table.
Different repositories expose different parts of the information:
- GEO explains what the biological samples are
- SRA explains the sequencing structure
- ENA gives convenient FASTQ download links

That is why the scripts first parse GEO metadata, then resolve SRA accessions, then ask ENA for downloadable files.

## Useful Links

### GEO series
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001

### GEO SOFT metadata
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103001/soft/GSE103001_family.soft.gz

### NCBI E-utilities: SRX -> SRR
Template:
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=<SRX>&rettype=runinfo&retmode=text`

Example:
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRX3120292&rettype=runinfo&retmode=text

### ENA: SRR -> FASTQ URLs
Template:
- `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<SRR>&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes`

Example:
- https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR5962198&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes
