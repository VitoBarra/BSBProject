BIOINFORMATICS AND SYSTEMS BIOLOGY (BSB)

ASSIGNMENT TO:

Vito Barra

INSPIRING PAPER:

Transcriptome-wide analysis of natural antisense transcripts shows their potential role in breast

cancer

https://doi.org/10.1038/s41598-017-17811-2

PART 1 (Prof. Galfrè):

Dataset: GSE103001

Breast cancer is a heterogeneous disease in which transcriptomic alterations affect both protein-

coding and non-coding genes. Estrogen receptor-positive (ER+) breast tumors represent a major

clinical subtype, and the comparison between tumor tissue and matched adjacent non-malignant

tissue provides a useful framework for identifying disease-associated transcriptional programs

while controlling, at least in part, for inter-patient variability. In this assignment, you will reproduce

a simplified gene-level RNA-seq analysis on ER+ breast cancer samples, focusing on differential

expression and pathway enrichment.

Goal

1.

Reproduce a scaled-down version of the paper’s tumor-versus-adjacent-tissue RNA-seq

analysis, but using a more recent computational setup:

2.

download FASTQ

3.

quantify with Salmon

4.

perform DEA (tumor vs adjacent non-malignant tissue)

5.

perform gene/pathway enrichment

6.

compare your results with the published conclusions, and discuss how the use of Salmon

+ GRCh38 may affect the outcome relative to the original STAR + HTSeq + GRCh37

analysis.

Dataset

Use GSE103001.

Comparison: ER+ breast tumor vs matched adjacent non-malignant tissue.

Raw sequencing data are available through SRA SRP116023 / BioProject PRJNA399721.

To keep the analysis manageable, use 8 samples total, corresponding to 4 matched pairs.

Important note

This is a paired design: each tumor sample must be compared with the adjacent tissue from the

same patient. Your statistical analysis must therefore account for patient pairing. The original

study also used a paired design.

Tasks

1. Download and basic QC

2. Quantification with Salmon

Quantify transcript abundance with Salmon, using a current GRCh38-based human

transcriptome annotation such as a recent GENCODE or Ensembl release.

3. Differential expression analysis

Perform differential expression analysis using DESeq2 or edgeR. Because the dataset is paired,

use a design that accounts for the patient effect, for example:

4. Enrichment analysis

Perform gene or pathway enrichment analysis on the significant DEGs, or on a ranked gene list.

5. Comparison and discussion

Briefly compare your results with the original paper.

You are not required to reproduce the full antisense-transcript analysis of the paper. The focus of

this assignment is the gene-level tumor-versus-normal comparison and its biological

interpretation. The original paper’s broader aim was to study natural antisense transcripts in ER+

breast cancer.

Deliverables

Submit:

•

a short report

•

the scripts or notebook used for the analysis

•

the DEG table

•

the enrichment results and discussion

PART 2 (Prof. Milazzo):

Analyze the pathology considered in the paper (Breast Cancer - BC) by using a network medicine

approach (i.e., investigation of PPI network through network theory methods).

Objective: Investigate the relationship between BC and another pathology: Rheumatoid Arthritis

(RA).

What to do:

•

Create a PPI network (from STRING, by filtering interaction types as done in class) by

including the differently expressed genes identified in the paper (the 72 genes listed in table

S11 in Supplemental File 1). Rank the 72 genes on the basis of network-based node

importance measures seen in class (centralities, etc.). In addition, apply community

detection methods to identify significant cluster in the PPI network. Use the obtained

rankings and clusters to select a small subset of (10-20) “important” genes out of the given

72 in a reasonable way. Such a subset of genes will constitute the BC disease module

•

Create a RA disease module by including RA relevant genes obtained through the

DISGENET database (RA identifier in DISGENET = C0003873). Consider Diabetes Mellitus

(DM) as a “control” disease, creating a disease module taking relevant genes from

DISGENET (DM identifier in DISGENET = C0011849). Construct a PPI network involving all

diseases (or two PPI networks with BC+RA and BC+DM) and compute network separation

of diseases on it (or them).

Notes:

•

the analysis can be done in Python and/or by using any other useful tool/language/library.

•

support by AI assistants (ChatGPT, etc…) is welcome, but be ready to answer questions on

the AI-generated code, text, images, …
