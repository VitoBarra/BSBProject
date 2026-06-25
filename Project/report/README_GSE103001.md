# GSE103001 Data Acquisition

This folder contains a reproducible path to rebuild the metadata table for the RNA-seq assignment dataset and download the selected FASTQ files.

## Files
- `GSE103001_family.soft.gz`: GEO series metadata file for GSE103001.
- `build_metadata_table.py`: rebuilds a TSV with matched patient pairs, SRX/SRR accessions, FASTQ URLs, and file sizes.
- `download_GSE103001_4pairs.py`: downloads FASTQ files listed in a TSV.
- `GSE103001_selected_4pairs.tsv`: example TSV already built for 4 matched pairs.

## Manual Reproducible Procedure

### 1. GEO series page
Open the GEO series page:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103001

This page describes the study and lists the 44 samples.

### 2. GEO SOFT file
Open or download the GEO SOFT metadata file:
- https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103001/soft/GSE103001_family.soft.gz

This file contains, for each sample:
- the GSM accession
- the sample title such as `Pat_12-02_normal`
- the SRA experiment accession `SRX...`

The logic used in the script is:
- parse every `^SAMPLE` block
- read `!Sample_title`
- extract `patient_id` and `condition` from titles like `Pat_12-02_normal`
- read the `SRX...` accession from `!Sample_relation = SRA: ...`

### 3. From SRX to SRR using NCBI E-utilities
For each experiment accession `SRX...`, resolve the corresponding run accession `SRR...` with:
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=<SRX>&rettype=runinfo&retmode=text`

Example:
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=SRX3120292&rettype=runinfo&retmode=text

In the returned CSV, read the `Run` column.

### 4. From SRR to FASTQ URLs using ENA
For each run accession `SRR...`, resolve the direct FASTQ URLs with:
- `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<SRR>&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes`

Example:
- https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR5962198&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes

The output includes:
- `fastq_ftp`: semicolon-separated FASTQ paths
- `fastq_bytes`: semicolon-separated file sizes

To download manually, prepend `https://` to each `fastq_ftp` entry.

### 5. Matched-pair selection rule used here
The script selects patients in lexicographic order and keeps only complete matched pairs:
- a pair is valid only if both `normal` and `tumor` are present
- `--num-pairs N` keeps the first `N` complete pairs

For `--num-pairs 4`, the selected patients are:
- `12-02`
- `12-03`
- `13-02`
- `13-03`

## Automatic Rebuild

Build a metadata table for 4 pairs:

```powershell
python .\DataSourcer\build_metadata_table.py --num-pairs 4
```

Build a metadata table for 6 pairs and choose the output name explicitly:

```powershell
python .\DataSourcer\build_metadata_table.py --num-pairs 6 --output .\data\GSE1030001\GSE103001_selected_6pairs.tsv
```

## Download FASTQ Files

Download the FASTQ files listed in the default 4-pair TSV:

```powershell
python .\DataSourcer\download_fastq_from_tsv.py
```

Download FASTQ files from a different TSV and destination:

```powershell
python .\DataSourcer\download_fastq_from_tsv.py --metadata .\GSE103001_selected_6pairs.tsv --dest .\raw_fastq\GSE103001_6pairs
```

## Notes
- FASTQ files are raw sequencing reads, not trimmed data.
- You should still run QC and decide whether trimming is needed.
- One run in this dataset may expose an additional singleton FASTQ file besides `_1` and `_2`; this comes from the archive export, not from preprocessing done by the authors.

