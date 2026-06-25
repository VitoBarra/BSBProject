# Part 2 — Network Medicine Analysis of Breast Cancer, Rheumatoid Arthritis and Diabetes Mellitus

## Objective

This analysis addresses Part 2 of the assignment by prioritizing the 72 Breast
Cancer (BC) genes reported in Supplemental Table S11 and comparing the final BC
disease module with Rheumatoid Arthritis (RA). Diabetes Mellitus (DM) is the
control disease.

The analysis has two deliberately different network scopes:

1. an induced 72-gene network for ranking and clustering the BC candidates;
2. a broader physical PPI background for BC–RA and BC–DM separation.

## Data and methods

### Disease genes

- **BC:** 72 protein-coding genes from Table S11 of Wenric et al. (2017),
  extracted without dropping any candidate.
- **RA:** 225 unique genes from all curated DisGeNET associations for
  `C0003873`.
- **DM:** 88 unique genes from the same curated DisGeNET policy for
  `C0011849`.

The DisGeNET associations were downloaded through the Academic API on
2026-08-09. The current database release at the time of finalization was
[DisGeNET 25.4](https://www.disgenet.com/docs?version=25.4). No minimum
association-score cutoff was imposed, so RA and DM were treated consistently
without selecting an arbitrary threshold.

### STRING policy

The network policy follows Lectures 29 and 30:

```text
Organism:                 Homo sapiens (NCBI taxonomy 9606)
STRING version:           12.0
Network type:             physical
Evidence retained:        experiments and curated databases
Evidence excluded:        text mining and other predicted/indirect channels
Minimum confidence:       0.70
BC additional proteins:   0
Separation background:    up to 650 additional interaction partners
```

The stable STRING 12.0 endpoint is used for reproducibility. STRING's API
returns individual channel scores but does not expose the web interface's
evidence checkboxes. Therefore, the experimental (`escore`) and database
(`dscore`) scores were combined using STRING's probabilistic procedure,
including removal and restoration of its prior probability (p=0.041). Edges
with the resulting score below 0.70 were discarded. The original all-channel
score is retained in the edge table but is not used for ranking. See the
[official STRING scoring documentation](https://string-db.org/help/scores/).

### Centralities and communities

The 72-gene undirected network was analyzed using:

- degree and degree centrality;
- unweighted betweenness centrality;
- normalized harmonic centrality, which remains defined for disconnected
  graphs;
- PageRank weighted by the class-filtered STRING score.

Each measure was converted to a within-network percentile. The composite score
was the unweighted mean of the four percentiles. Weighted Louvain community
detection used seed 42.

Community structure was tested against 500 degree-preserving randomized
networks. Each null network underwent edge swaps, observed edge weights were
shuffled onto the randomized topology, and Louvain communities were detected
again. The empirical upper-tail p-value was calculated with a plus-one
correction.

### BC module selection

Strict class filtering left only 19 non-isolated BC genes. Selecting all of
them would not constitute meaningful prioritization, so the final target was
set to 10 genes. The rule was:

1. take the highest composite-ranked representative of every nontrivial
   Louvain community;
2. fill the remaining positions by composite centrality;
3. do not select isolates solely because of tied zero centralities;
4. document a topological and biological rationale for every selected gene.

### Disease separation

The final BC, RA and DM genes were mapped onto a shared physical STRING
background. Additional proteins provide paths but are not disease-module
members.

Following Lecture 29, module distance is the all-pairs average shortest-path
distance:

$$
\langle d_{XY}\rangle =
\frac{1}{|M_X||M_Y|}
\sum_{i\in M_X}\sum_{j\in M_Y}d_{ij}.
$$

Disease separation is:

$$
S_{AB}=\langle d_{AB}\rangle-
\frac{\langle d_{AA}\rangle+\langle d_{BB}\rangle}{2}.
$$

Self-pairs are included in intra-module Cartesian products, consistently with
the denominator in the lecture formula. Unreachable pairs are excluded from
the mean and reported explicitly.

## BC network results

### Mapping and topology

All 72 genes mapped to STRING. At the class-aligned threshold the induced
network contained:

```text
Nodes:                         72
Edges:                         14
Density:                       0.00548
Connected components:          61
Largest connected component:    4 genes
Isolated genes:                53
Weighted average clustering:    0.0753
```

The strong fragmentation is an expected consequence of requiring direct,
high-confidence physical evidence among only the 72 submitted proteins. An
isolate is not biologically unimportant; it simply has no accepted edge to
another Table S11 protein under this policy.

### Ranking

The leading genes were:

| Rank | Gene | Degree | Composite percentile | Community |
|---:|---|---:|---:|---:|
| 1 | MSH2 | 3 | 0.9930 | 1 |
| 2 | MSH6 | 3 | 0.9930 | 1 |
| 3 | RUNX1 | 2 | 0.8327 | 2 |
| 4 | TAL1 | 2 | 0.8292 | 2 |
| 5 | XPC | 2 | 0.7887 | 1 |
| 6 | BRCA1 | 2 | 0.7852 | 1 |
| 7 | FLI1 | 2 | 0.7835 | 2 |
| 8 | CDKN2A | 1 | 0.7500 | 3 |
| 9 | CHD4 | 1 | 0.7500 | 4 |
| 10 | EGFR | 1 | 0.7500 | 5 |

Complete rankings are in
[bc_gene_centralities.tsv](../data/network_medicine/results/bc/bc_gene_centralities.tsv).

### Communities and significance

The 19 connected genes formed eight nontrivial communities:

| Community | Members |
|---:|---|
| 1 | MSH2, MSH6, XPC, BRCA1 |
| 2 | RUNX1, TAL1, FLI1 |
| 3 | CDKN2A, MDM2 |
| 4 | CHD4, H3F3A |
| 5 | EGFR, PIK3CA |
| 6 | GNAS, TSHR |
| 7 | IL6ST, LIFR |
| 8 | MTCP1, TCL1A |

Observed weighted modularity was 0.7996. Across 500 degree-preserving null
networks, mean modularity was 0.7165 (SD 0.0321), giving \(z=2.58\) and an
empirical upper-tail \(p=0.0040\). The observed partition is therefore more
modular than expected under this null model, although inference is limited by
the network's very small number of edges.

### Final 10-gene BC module

| Gene | Community | Selection and biological rationale |
|---|---:|---|
| MSH2 | 1 | Highest-ranked repair-community hub; mismatch-repair partner of MSH6. |
| MSH6 | 1 | Tied highest composite score; forms MutSα with MSH2. |
| RUNX1 | 2 | Top transcription-factor-community representative and cancer gene. |
| TAL1 | 2 | Second high-ranked member of the RUNX1–TAL1–FLI1 community. |
| CDKN2A | 3 | Cell-cycle tumour suppressor representing the CDKN2A–MDM2 component. |
| CHD4 | 4 | NuRD chromatin-remodelling factor representing CHD4–H3F3A. |
| EGFR | 5 | Receptor tyrosine kinase relevant to breast-cancer signalling; represents EGFR–PIK3CA. |
| GNAS | 6 | Signal-transduction representative of GNAS–TSHR. |
| IL6ST | 7 | Shared cytokine-receptor subunit representing IL6ST–LIFR. |
| MTCP1 | 8 | Cancer-gene representative of MTCP1–TCL1A. |

The repair and cell-cycle interpretations are consistent with descriptions in
the [NCBI Gene database](https://www.ncbi.nlm.nih.gov/gene/) and with published
discussion of MSH2 and CDKN2A in hereditary breast-cancer gene panels
([reviewed study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7723566/)). The
annotations aid interpretation but did not replace the reproducible network
selection rule.

The final machine-readable module is
[bc_module_selected.tsv](../data/network_medicine/results/bc/bc_module_selected.tsv).

## Disease-module network

### Mapping and background

| Module | Input genes | STRING-mapped genes |
|---|---:|---:|
| BC | 10 | 10 |
| RA | 225 | 222 |
| DM | 88 | 87 |

There were 309 unique submitted genes and 306 mapped. The unmapped genes were
`AP4B1-AS1`, `GSTT1`, and `LOC128462409`; none was silently discarded.

The final background contained:

```text
Nodes:                         836
Edges:                       2,369
Background connector nodes:   530
Connected components:          151
Largest component:             638 nodes
Isolated nodes:                117
```

BC overlaps RA directly at `RUNX1` and `IL6ST`; BC and DM have no direct gene
overlap.

### Network separation

| Comparison | \(d_{AB}\) | \(d_{AA}\) | \(d_{BB}\) | \(S_{AB}\) | Finite cross pairs |
|---|---:|---:|---:|---:|---:|
| BC–RA | 4.6473 | 3.3171 | 4.5099 | 0.7339 | 1,035/2,220 (46.6%) |
| BC–DM | 4.4317 | 3.3171 | 4.1538 | 0.6963 | 315/870 (36.2%) |

Both separations are positive, so BC occupies a distinct topological region
from both comparison modules under this network policy. BC–DM has the smaller
separation by 0.0376; therefore, this analysis does **not** support the
hypothesis that BC is topologically closer to RA than to the DM control. The
difference is small and no randomized hypothesis test was performed for the
difference between the two separation values.

Expanding the background increased finite cross-path coverage relative to the
150-interactor analysis, but 53.4% of BC–RA and 63.8% of BC–DM pairs remain
unreachable. This is the principal limitation and reflects both incomplete
interactome coverage and the strict physical/evidence filter. Excluding
unreachable pairs avoids imposing an arbitrary finite penalty, but it means
the reported distances describe reachable pairs rather than the entire
Cartesian product.

## Figures and result files

- [BC network coloured by community](../data/network_medicine/results/bc/bc_network.png)
- [Largest BC component and selected genes](../data/network_medicine/results/bc/bc_largest_component.png)
- [BC composite ranking](../data/network_medicine/results/bc/bc_top20_composite_rank.png)
- [Community modularity null model](../data/network_medicine/results/bc/bc_modularity_null.png)
- [Representative BC-to-RA/DM shortest-path subnetwork](../data/network_medicine/results/separation/bc_ra_dm_joint_network.png)
- [Separation comparison](../data/network_medicine/results/separation/disease_separation.png)
- [Finite-path coverage](../data/network_medicine/results/separation/finite_path_coverage.png)
- [Separation table](../data/network_medicine/results/separation/disease_separation.tsv)

GraphML versions of both networks are provided for inspection in Cytoscape.

## Reproducibility

The final run used Python 3.12.3, NetworkX 3.6.1, Matplotlib 3.10.9, STRING
12.0, Louvain seed 42, and `PYTHONHASHSEED=0`. Raw responses, clean mappings, filtered edges,
metadata JSON files, null-model values, figures and GraphML files are retained.

Run the analysis inside the project Docker environment:

```bash
make docker-build
make network-bc
make network-separation
```

DisGeNET data can be refreshed separately, using the ignored API-key file:

```bash
make network-disgenet-download
```

## Conclusion

The strict, class-aligned physical-interaction policy produced a sparse but
significantly modular BC candidate network. The reproducible community-aware
selection rule prioritized a 10-gene module spanning DNA repair,
transcriptional regulation, cell-cycle control, chromatin remodelling,
signalling and cytokine-receptor components. In the broader interactome, BC
was separated from both RA and DM. BC–DM was marginally closer than BC–RA, so
the results do not support the proposed greater network proximity between BC
and RA. This negative result is biologically interpretable but should be read
in light of incomplete finite-path coverage and the restricted evidence
policy.

## Limitations

1. STRING is incomplete and biased toward well-studied proteins.
2. Physical, experiments/databases-only filtering improves interpretability of
   path length but produces a sparse BC network.
3. The 650 additional proteins approximate the broader class network; they are
   connector nodes, not disease genes.
4. RA and DM module sizes differ because the same evidence policy was applied
   without an arbitrary score cutoff.
5. Unreachable pairs remain substantial and are excluded from distance means.
6. The module is a network-based prioritization, not proof that its genes are
   causal drivers or clinical biomarkers.

## References

1. Wenric S, ElGuendi S, Caberg JH, et al. *Transcriptome-wide analysis of
   natural antisense transcripts shows their potential role in breast cancer*.
   Scientific Reports. 2017;7:17452.
   [doi:10.1038/s41598-017-17811-2](https://doi.org/10.1038/s41598-017-17811-2).
2. Szklarczyk D, Kirsch R, Koutrouli M, et al. *The STRING database in 2023:
   protein–protein association networks and functional enrichment analyses for
   any sequenced genome of interest*. Nucleic Acids Research.
   2023;51(D1):D638–D646.
   [doi:10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000).
3. Piñero J, Ramírez-Anguita JM, Saüch-Pitarch J, et al. *The DisGeNET
   knowledge platform for disease genomics: 2019 update*. Nucleic Acids
   Research. 2020;48(D1):D845–D855.
   [doi:10.1093/nar/gkz1021](https://doi.org/10.1093/nar/gkz1021).
4. Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E. *Fast unfolding of
   communities in large networks*. Journal of Statistical Mechanics: Theory
   and Experiment. 2008;2008(10):P10008.
   [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008).
5. Menche J, Sharma A, Kitsak M, et al. *Uncovering disease-disease
   relationships through the incomplete interactome*. Science.
   2015;347(6224):1257601.
   [doi:10.1126/science.1257601](https://doi.org/10.1126/science.1257601).
