# Guide to Supplemental Table S11

## Source

Table S11 belongs to Supplemental File 1 of:

> Wenric, S. et al. (2017), *Transcriptome-wide analysis of natural
> antisense transcripts shows their potential role in breast cancer*,
> Scientific Reports 7, 17452.

- DOI: <https://doi.org/10.1038/s41598-017-17811-2>
- Local paper: [Transcriptome-wide analysis of natural antisense transcripts shows their potential role in breast cancer.pdf](./Transcriptome-wide%20analysis%20of%20natural%20antisense%20transcripts%20shows%20their%20potential%20role%20in%20breast%20cancer.pdf)
- Local Supplemental File 1: [41598_2017_17811_MOESM1_ESM.pdf](./41598_2017_17811_MOESM1_ESM.pdf)
- Extracted gene list: [bc_genes_s11.tsv](../data/network_medicine/processed/bc_genes_s11.tsv)

Table S11 occupies PDF pages 54–62 of Supplemental File 1. The gene
symbols and Ensembl identifiers are shown on PDF pages 54, 57, and 60.

## What Table S11 represents

Table S11 reports 72 protein-coding cancer genes selected through the
paper's analysis of natural antisense transcription in breast cancer.
It is an annotated candidate-gene table, not a standard differential
expression results table.

The authors studied pairs formed by:

- a protein-coding transcript (PCT);
- an overlapping non-coding natural antisense transcript (ncNAT).

They generated candidate lists using three analyses:

1. **DiffCorr**: the correlation between the PCT and its ncNAT changes
   between tumour and matched adjacent non-malignant tissue.
2. **ncNATDiffExp**: the ncNAT is differentially expressed between
   tumour and adjacent tissue.
3. **VarRatio**: the expression ratio between the PCT and ncNAT changes
   substantially between the two conditions.

The protein-coding genes obtained from these analyses were compared with
the COSMIC Cancer Gene Census. Their intersection produced the 72 genes
reported in Table S11:

```text
Genes selected by DiffCorr, ncNATDiffExp, or VarRatio
                            ∩
                COSMIC Cancer Gene Census
                            =
              72 genes reported in Table S11
```

Each S11 gene therefore met two conditions:

1. it exhibited an interesting coding/antisense deregulation pattern in
   the study;
2. it was recognized as a cancer-associated gene by the version of the
   COSMIC Cancer Gene Census used by the authors.

This does not mean that all 72 genes are specific causes of breast
cancer. Some are better known for their roles in other tumour types.

## Experimental results versus COSMIC prior knowledge

Table S11 combines two distinct sources of information:

```text
The paper's RNA-seq analyses
→ identify genes with altered PCT–ncNAT behaviour

COSMIC Cancer Gene Census
→ supplies pre-existing knowledge about cancer involvement

Intersection and annotation
→ produce the 72 genes and the descriptive columns in Table S11
```

The paper contributes the experimental selection of candidate genes.
COSMIC contributes general knowledge accumulated from cancer research,
including whether alterations of a gene have been observed as somatic or
germline, associated tumour types, roles in cancer, mutation categories,
and translocation partners.

The authors did not first establish that these were the 72 "most
important" genes. They identified protein-coding genes with interesting
antisense-expression behaviour and retained those already present in the
Cancer Gene Census. The result is therefore a biologically motivated
candidate list rather than a definitive importance ranking.

Most importantly, the COSMIC columns are not patient-level mutation
results from the 22 patients analysed in the paper. An entry such as
`BRCA1 — Germline: yes` means that pathogenic germline variants of
`BRCA1` were already known to confer cancer susceptibility. It does not
mean that the authors detected a germline BRCA1 variant in one or more
patients in this cohort.

## Why the table is divided across nine pages

Table S11 contains many columns and could not fit on one page. Its
columns are consequently divided into three blocks, each repeated for
the same groups of genes:

- pages 1, 4, and 7 of Table S11: gene identity and tumour types;
- pages 2, 5, and 8: cancer role and mutation information;
- pages 3, 6, and 9: aliases and additional syndrome information.

These pages describe the same 72 genes; they are not separate gene
lists.

## Explanation of the columns

### Gene Symbol

The conventional symbol used to identify the gene, for example `BRCA1`,
`EGFR`, `PIK3CA`, `RB1`, or `RUNX1`.

These symbols will be mapped to human proteins in STRING. Because the
table was published in 2017, all symbols and aliases must be checked
against the current STRING version.

### Ensembl PC

The Ensembl identifier of the protein-coding gene. `PC` means
protein-coding.

Examples:

```text
BRCA1 → ENSG00000012048
EGFR  → ENSG00000146648
```

Stable identifiers are valuable when gene symbols have aliases or have
changed over time.

### Somatic

Indicates whether somatic alterations in the gene were documented in
cancer. Somatic alterations are acquired during a person's lifetime,
occur in tumour cells, and are not normally inherited by offspring.

### Germline

Indicates whether inherited, disease-associated alterations were
documented. Germline alterations can be transmitted to offspring and
may confer a hereditary predisposition to cancer.

For example, `BRCA1` has both somatic and germline involvement in the
table.

Everyone normally inherits genes such as `BRCA1`, `RB1`, and `MSH2`.
The `Germline` column does not indicate whether the gene itself was
inherited: it indicates that pathogenic variants of that gene are known
to occur in the germline.

A germline variant is present in the egg or sperm from which the person
developed, or arises very early in embryonic development. It is therefore
normally present in almost every cell. It may have been inherited from a
parent or may be de novo. If the person has children, the variant may be
transmitted to them.

Such a variant usually represents cancer predisposition rather than a
certainty that cancer will develop.

```text
Normal gene inherited by everyone
≠
Pathogenic germline variant inherited by a particular person
```

In tumour-suppressor genes, an inherited pathogenic variant may act as a
first hit. A later somatic alteration affecting the remaining functional
copy in a particular cell can provide a second hit and contribute to
tumour development.

### How somatic and germline variants are distinguished

A reference genome provides common coordinates and a baseline for
detecting sequence differences, but a difference from the reference is
not automatically pathogenic. Every individual carries many benign
variants.

To identify candidate somatic alterations, tumour DNA is ideally
compared with matched normal DNA from the same patient:

```text
Reference genome
       ↓
Patient's normal DNA → inherited and individual variants
       ↓
Patient's tumour DNA compared with matched normal DNA
       ↓
Alterations found only in the tumour → candidate somatic variants
```

A germline variant should also be detectable in non-tumour DNA, commonly
obtained from blood or saliva. A somatic variant is normally restricted
to the tumour or to a subset of cells.

Detected variants require interpretation. Depending on the available
population, clinical, functional, and segregation evidence, a variant
may be classified as benign, likely benign, of uncertain significance,
likely pathogenic, or pathogenic.

These principles explain how mutations can be identified in genomic
studies, but Table S11 does not report this kind of matched
tumour–normal variant calling for the study participants. Its mutation
columns summarize previously known COSMIC annotations.

### Tumour Types (Somatic)

Lists tumour types in which somatic alterations of the gene were known.
These are not limited to breast cancer because the annotation comes from
the general COSMIC Cancer Gene Census.

Examples include:

```text
BRCA1  → ovarian
PIK3CA → colorectal, gastric, glioblastoma, breast
RB1    → retinoblastoma, sarcoma, breast, small-cell lung carcinoma
```

### Tumour Types (Germline)

Lists tumour types associated with inherited alterations of the gene.
For example, germline `BRCA1` alterations are associated with breast and
ovarian cancer.

### Cancer Syndrome

Names a hereditary cancer syndrome associated with the gene.

Examples:

```text
BRCA1  → hereditary breast/ovarian cancer
CDKN2A → familial malignant melanoma
MSH2   → hereditary non-polyposis colorectal cancer
RB1    → familial retinoblastoma
```

A blank cell does not mean that the gene has no cancer relevance. It
means that no corresponding syndrome was reported in that column.

### Molecular Genetics

Describes the inheritance or action pattern used in the Cancer Gene
Census:

- `Dom`: dominant;
- `Rec`: recessive;
- `Rec/X`: X-linked recessive.

This field does not indicate whether gene expression is upregulated or
downregulated.

### Role in Cancer

Describes the established functional role of the gene:

- `oncogene`: activation can promote cancer;
- `TSG`: tumour-suppressor gene;
- `oncogene/TSG`: the gene can have either role depending on context.

Examples:

```text
BRCA1  → TSG
CDKN2A → TSG
PIK3CA → oncogene
RUNX1  → oncogene/TSG
```

### Mutation Types

Reports the types of cancer-associated alterations annotated by COSMIC.
Common abbreviations include:

- `A`: amplification;
- `D`: large deletion;
- `F`: frameshift mutation;
- `Mis`: missense mutation;
- `N`: nonsense mutation;
- `S`: splice-site mutation;
- `T`: translocation;
- `O`: other mutation type.

Examples:

```text
BRCA1  → D, Mis, N, F, S
MDM2   → A
PIK3CA → Mis
```

### Translocation Partner

For genes involved in chromosomal rearrangements, this column identifies
known fusion or translocation partners.

Examples:

```text
AFF1 → KMT2A
FLI1 → EWSR1
HLF  → TCF3
LIFR → PLAG1
```

### Other Germline Mutation and Other Syndrome

These fields describe inherited alterations and syndromes that are not
necessarily classified as cancer syndromes. For example, `HNF1A` is
associated with maturity-onset diabetes of the young.

This demonstrates that a cancer-associated gene may also contribute to
other diseases.

### Synonyms

Contains alternative symbols and database identifiers, including:

- older gene symbols;
- protein names;
- Entrez Gene identifiers;
- Ensembl identifiers;
- UniProt identifiers.

For example, `ERBB1` is an alias of `EGFR`. This field is useful when
mapping the historical gene list to current databases.

## What Table S11 does not show

Table S11 does not provide:

- log2 fold changes;
- expression counts;
- p-values or adjusted p-values;
- whether a coding gene is upregulated or downregulated;
- STRING interactions;
- centrality measurements;
- community assignments;
- a ranking of the 72 genes;
- proof that every gene causes breast cancer.

It must therefore not be interpreted as the output of the Part 1 DESeq2
analysis. It is a curated candidate list produced by the original
paper's ncNAT analyses and its comparison with COSMIC.

## Role in Part 2 of the assignment

Table S11 provides the starting genes for the Breast Cancer network:

```text
72 genes from Table S11
            ↓
Map identifiers and retrieve interactions from STRING
            ↓
Construct the Breast Cancer PPI network
            ↓
Calculate node centralities
            ↓
Detect network communities
            ↓
Select 10–20 important genes
            ↓
Final Breast Cancer disease module
```

The S11 cancer annotations can support the biological interpretation of
the final genes. They must not replace the requested network analysis.
For example, a familiar cancer gene should not be selected merely
because it is well known: selection should also be justified through its
centrality, community membership, bridge role, and relevance to the
network.

The final Breast Cancer module will subsequently be compared with:

- a Rheumatoid Arthritis module obtained from DisGeNET identifier
  `C0003873`;
- a Diabetes Mellitus control module obtained from DisGeNET identifier
  `C0011849`.

The comparison will be performed by mapping the disease genes onto a
common PPI network and calculating network separation.

## Data-quality considerations

Before querying STRING:

1. preserve the original gene symbol and Ensembl identifier;
2. map both identifiers against the current STRING database;
3. record updated symbols and aliases separately;
4. retain unmapped and ambiguous genes in a mapping report;
5. do not silently discard isolated proteins;
6. document the STRING version, organism, evidence filters, network
   type, and confidence threshold.

The cleaned project table contains exactly:

```text
72 data rows
72 unique gene symbols
72 unique Ensembl identifiers
```
