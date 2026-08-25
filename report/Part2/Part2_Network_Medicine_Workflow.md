# Part 2 — Network Medicine Workflow

## Objective

Investigate the relationship between Breast Cancer (BC) and Rheumatoid
Arthritis (RA) using protein–protein interaction networks. Diabetes
Mellitus (DM) is used as a control disease.

The work is divided into two analyses:

1. prioritize the 72 BC genes and select a final BC module of 10–20
   genes;
2. compare the selected BC module with RA and DM through network
   separation.

These analyses use different networks and must not be confused.

## Input data

### Breast Cancer

The BC input consists of the 72 protein-coding genes reported in Table
S11 of Supplemental File 1.

- Source explanation: [Table_S11_Guide.md](./Table_S11_Guide.md)
- Extracted genes:
  [bc_genes_s11.tsv](../data/network_medicine/processed/bc_genes_s11.tsv)

The genes were selected by intersecting:

```text
Genes with altered PCT–ncNAT behaviour in the breast-cancer study
                              ∩
              COSMIC Cancer Gene Census genes
                              =
                         72 BC genes
```

### Rheumatoid Arthritis

- Source: DisGeNET
- Disease identifier: `C0003873`

### Diabetes Mellitus

- Source: DisGeNET
- Disease identifier: `C0011849`
- Role: control disease

## Phase 1 — Prioritize the 72 BC genes

### Purpose

Select 10–20 topologically important genes from the 72 candidates.
This selected subset will constitute the final BC disease module.

### Network scope

The complete human interactome is not required for this phase.

```text
Nodes = the 72 Table S11 genes
Edges = accepted STRING interactions among those genes
Additional STRING interactors = 0
```

This is the **BC prioritization network**.

### Step 1. Validate and map identifiers

Use both the published gene symbols and Ensembl identifiers to map the
72 genes to current human STRING proteins.

Record:

- original gene symbol;
- original Ensembl identifier;
- current mapped symbol;
- STRING protein identifier;
- mapping status;
- aliases or corrections;
- ambiguous or unmapped identifiers.

No gene should be silently removed.

Expected output:

```text
data/network_medicine/processed/bc_string_mapping.tsv
```

### Step 2. Retrieve STRING interactions

Query STRING using:

```text
Organism: Homo sapiens
Taxonomy ID: 9606
Input: mapped Table S11 genes
Additional interactors: 0
```

The course method was confirmed from Lectures 29 and 30:

- physical STRING subnetwork;
- experiments and databases evidence channels only;
- text mining and other indirect/predicted channels excluded;
- high-confidence threshold of 0.70;
- no additional interactors for prioritizing the original 72 BC genes.

Because the STRING API does not expose the web interface's evidence-channel
checkboxes, the implementation retrieves physical associations and calculates
an experiments-plus-databases score using STRING's probabilistic combination
procedure, including removal and restoration of its prior probability of
0.041. It retains edges with this score at least 0.70. The original channel
scores and API all-channel combined score remain in the output.

Save the original response and the filtered edge list separately.

Expected outputs:

```text
data/network_medicine/raw/string_bc_interactions.tsv
data/network_medicine/processed/bc_ppi_edges.tsv
```

### Step 3. Construct the induced BC network

Build an undirected PPI network in which:

- nodes are the 72 BC genes;
- edges are accepted STRING interactions;
- isolated mapped genes remain in the node table;
- optional edge weights contain STRING confidence scores.

Report:

- number of input genes;
- successfully mapped genes;
- unmapped or ambiguous genes;
- number of nodes and edges;
- isolated nodes;
- number and size of connected components;
- largest connected component;
- network density;
- average clustering coefficient;
- degree distribution.

### Step 4. Calculate centralities

For every BC gene, calculate:

- degree centrality;
- betweenness centrality;
- closeness or harmonic centrality;
- eigenvector centrality or PageRank.

Interpretation:

| Measure | Network role |
|---|---|
| Degree | Local hub with many direct interactions |
| Betweenness | Bridge or bottleneck between network regions |
| Closeness/harmonic | Proximity to the rest of the network |
| Eigenvector/PageRank | Connection to other influential nodes |

Disconnected components must be handled explicitly. Harmonic centrality
or component-aware closeness is preferable to a formula that assumes a
fully connected graph. PageRank can be used if eigenvector centrality is
unstable on a disconnected network.

Expected output:

```text
data/network_medicine/results/bc_gene_centralities.tsv
```

### Step 5. Detect communities

Apply a community-detection method used in class, such as Louvain.
Infomap may be used as a secondary comparison.

Record:

- community membership of every connected gene;
- community sizes;
- internal and external edges;
- modularity;
- representative high-centrality genes;
- isolated genes or trivial components.

Expected output:

```text
data/network_medicine/results/bc_communities.tsv
```

### Step 6. Select 10–20 genes

Do not select genes from a single centrality alone.

A reproducible selection strategy is:

1. convert each centrality to a rank or percentile;
2. calculate a transparent composite centrality score;
3. identify the most central genes in each nontrivial community;
4. include important high-betweenness bridges;
5. preserve representation of the main communities;
6. integrate the Table S11 cancer annotations for interpretation;
7. select 10 genes when strict filtering leaves a very sparse network, avoiding
   the uninformative selection of nearly every connected gene;
8. record a reason for every selection.

Possible composite score:

```text
mean(
    degree percentile,
    betweenness percentile,
    harmonic/closeness percentile,
    PageRank percentile
)
```

The COSMIC annotations support biological interpretation but must not
replace the network-based selection.

Expected output:

```text
data/network_medicine/results/bc_module_selected.tsv
```

Suggested columns:

```text
gene
composite_rank
degree_rank
betweenness_rank
closeness_or_harmonic_rank
pagerank
community
selected
selection_reason
```

## Phase 2 — Obtain RA and DM disease genes

### Purpose

Construct RA and DM gene sets using documented gene–disease
associations.

DisGeNET supplies disease-associated genes; it does not supply PPI
centrality rankings.

### Step 1. Retrieve associations

Retrieve:

```text
RA → C0003873
DM → C0011849
```

Save the unmodified responses.

Expected outputs:

```text
data/network_medicine/raw/disgenet_ra_associations.tsv
data/network_medicine/raw/disgenet_dm_associations.tsv
```

### Step 2. Apply a consistent evidence filter

The same policy must be used for RA and DM:

- accepted evidence sources;
- curated versus inferred associations;
- minimum association score, if used;
- gene identifier requirements;
- handling of obsolete or ambiguous symbols.

The primary analysis should use the same evidence threshold for both
diseases. If their module sizes differ substantially, an equal-size
top-ranked sensitivity analysis may be added without replacing the
primary result.

Expected outputs:

```text
data/network_medicine/processed/ra_genes_disgenet.tsv
data/network_medicine/processed/dm_genes_disgenet.tsv
```

## Phase 3 — Compare BC with RA and DM

### Disease modules

The three module definitions are:

```text
BC module = 10–20 genes selected from the 72-gene PPI analysis
RA module = filtered DisGeNET genes associated with C0003873
DM module = filtered DisGeNET genes associated with C0011849
```

RA and DM are not automatically PPI-prioritized. They are selected using
gene–disease association evidence and then mapped onto the PPI network.

### Comparison-network options

The assignment permits either:

```text
One network:
BC + RA + DM
```

or:

```text
Network 1: BC + RA
Network 2: BC + DM
```

Use the same STRING version and interaction-filtering policy for all
diseases.

### Initial induced-network approach

First construct networks containing the disease genes and the accepted
STRING edges among them.

Measure:

- mapping coverage;
- number of connected components;
- size of the largest connected component;
- proportion of module genes with finite paths to the other module.

If most relevant nodes are connected, these networks may be sufficient
for the requested comparison.

### Common-interactome approach

If the induced networks are too fragmented for meaningful shortest
paths, map the disease modules onto a common filtered human STRING
interactome:

```text
Common human PPI background
├── selected BC genes
├── RA genes
└── DM genes
```

The additional proteins provide paths but are not considered disease
module members.

Only this separation phase may require the broader interactome. It is
not needed for prioritizing the original 72 BC genes.

The reason for choosing the induced or background-interactome approach
must be documented using connectivity statistics.

## Phase 4 — Calculate network separation

Calculate:

```text
separation(BC, RA)
separation(BC, DM)
```

Use one explicitly documented separation definition consistently for
both comparisons.

The standard form is:

$$
s_{AB}=d_{AB}-\frac{d_{AA}+d_{BB}}{2}
$$

where:

- $d_{AB}$ summarizes shortest-path distances between modules A and B;
- $d_{AA}$ summarizes internal distances in A;
- $d_{BB}$ summarizes internal distances in B.

Following Lecture 29, these quantities use the all-pairs average shortest-path
distance. Closest-node distance belongs to the separate drug–disease proximity
method and is not used for disease–disease separation here.

Interpretation:

- $s_{AB}<0$: modules overlap or are topologically close;
- $s_{AB}\approx0$: weak separation;
- $s_{AB}>0$: modules occupy more distinct network regions.

Also report:

- direct gene overlap;
- Jaccard similarity;
- $d_{AA}$, $d_{BB}$, and $d_{AB}$;
- mapped and unmapped genes;
- disconnected genes;
- finite-path coverage.

A stronger optional analysis can compare the observed separation with
degree-matched randomized modules and report a z-score or empirical
p-value.

## Biological question

The final comparison asks:

> Is the selected Breast Cancer module topologically closer to
> Rheumatoid Arthritis than it is to the Diabetes Mellitus control
> module?

A smaller BC–RA separation than BC–DM would support greater topological
proximity between BC and RA in the selected PPI network. This suggests
potentially shared molecular mechanisms but does not, by itself, prove
causality or clinical comorbidity.

## Required figures

Recommended figures:

1. induced PPI network of the 72 BC genes;
2. BC network with node size representing centrality;
3. BC network coloured by community;
4. selected 10–20-gene BC module;
5. BC–RA and BC–DM comparison networks;
6. bar plot comparing the two separation values;
7. optional randomized null distributions.

## Required result tables

At minimum:

```text
bc_string_mapping.tsv
bc_ppi_edges.tsv
bc_gene_centralities.tsv
bc_communities.tsv
bc_module_selected.tsv
ra_genes_disgenet.tsv
dm_genes_disgenet.tsv
disease_module_mapping.tsv
network_separation_results.tsv
```

## Reproducibility requirements

Record:

- source and version of every database;
- access dates;
- gene identifiers and mappings;
- STRING organism and network type;
- STRING evidence channels and threshold;
- DisGeNET filtering policy;
- software and package versions;
- random seeds for community detection and null models;
- exact separation equation;
- rules for disconnected nodes.

Following `agentRule.md`, long-running analysis and environment-dependent
operations must run through Docker-backed Makefile targets.

Suggested targets:

```text
make network-bc
make network-disease-genes
make network-separation
make network-analysis
```

## Course-method decisions

The lecture material confirms physical STRING interactions, experiments and
databases evidence, confidence at least 0.70, Louvain or MCL community
detection, and all-pairs average shortest paths for disease separation.

For the separation network, the class demonstration added intermediate
interactors so that module proteins could be connected through the PPI
background. This project requests up to 650 additional STRING interactors,
approximating the 150 first-neighbour plus 500 second-neighbour expansion used
in the class demonstration. The API performs one background expansion rather
than exposing separate shells. These proteins provide paths but are never
treated as disease-module members.
