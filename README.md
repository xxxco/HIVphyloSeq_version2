# HIVphyloSeq V2

All analyses used publicly available sequences downloaded on Feb. 4, 2025, from the [LANL HIV Sequence Database](https://www.hiv.lanl.gov/): 496 full-length HIV-1 subtype B genomes (one per patient), sampled 2005–2024 in the United States, Canada, and 19 European countries. See [`sequence_selection.md`](sequence_selection.md) for the full filtering workflow, including how to download the raw data from LANL.

---

## Pipeline Overview

```
LANL sequences (all fields, no filters, 14 batches by Se ID range)
    │
    ▼
Sequence filtering (→ 496 full-length subtype B genomes, one per patient)
    │
    ▼
Extract 7 subregions (protease, RT, integrase, vif, env, nef, gag-pol)
    │
    ▼
Multiple sequence alignment (DECIPHER)
    │
    ▼
Maximum-likelihood phylogeny (RAxML-NG)
    │
    ├─► ClusterPicker ──┐
    │                   ├─► Multi-metric concordance (compare_cluster.R)
    └─► HIV-TRACE ──────┘               │
                                        ▼
                                vmeasure_results.tsv
```

---

## Software Requirements

| Tool | Version | Notes |
|------|---------|-------|
| R | ≥ 4.5 | |
| DECIPHER | Bioconductor | alignment |
| Biostrings | Bioconductor | FASTA I/O |
| clevr | CRAN | V-measure, Adjusted Rand Index |
| dplyr, readr, stringr, tibble, jsonlite, purrr | CRAN | data handling |
| RAxML-NG | 1.2.2 | phylogeny |
| ClusterPicker | 1.2.5 | clustering |
| HIV-TRACE | 0.9.2 | clustering |

All tools above (DECIPHER, RAxML-NG, ClusterPicker, HIV-TRACE, and their R/Python dependencies) were run inside a single Singularity/Apptainer container image, built to bundle every tool the pipeline needs at the pinned versions listed above, so each analysis step below runs from the same reproducible environment instead of a separately installed tool:
```bash
singularity shell /path/to/multiple-alignment-container.sif
```

---

## Step 1 — Sequence Filtering

Filtering was run against a folder of raw LANL downloads (see [`sequence_selection.md`](sequence_selection.md) for the download procedure and full filtering criteria).

Output: one row per patient, tagged with HXB2 coordinates, ready for subregion extraction.

---

## Step 2 — Subregion Extraction

Seven subregions plus the full-length genome are extracted by HXB2 coordinate (see the table in [`sequence_selection.md`](sequence_selection.md)) into one FASTA file per region.

---

## Step 3 — Multiple Sequence Alignment (DECIPHER)

Alignments were generated using `DECIPHER::AlignSeqs()` with iterative refinement (5 iterations, 3 refinements). Run once per region FASTA file:

```r
library(DECIPHER)
library(Biostrings)

sequences <- readDNAStringSet("path/to/region_sequences.fasta")
aligned   <- AlignSeqs(sequences, iterations = 5, refinements = 3, processors = NULL)
writeXStringSet(aligned, "DECIPHER_aligned_region_sequences.fasta")
```

---

## Step 4 — Maximum-Likelihood Phylogeny (RAxML-NG)

Run once per region. The example below is for `Env`; repeat for each of the 8 regions (7 subregions + full-length), adjusting `--msa` and `--prefix` accordingly.

```bash
raxml-ng --all \
    --msa    path/to/phy_file/DECIPHER_aligned_Env_sequences.phy \
    --model  GTR+G \
    --prefix path/to/tree/Env/bootstrap_1000/Env_mytree_bs \
    --seed   123 \
    --bs-trees 1000 \
    --threads auto{N}
```

**Model**: GTR + Gamma-distributed rate heterogeneity
**Bootstrap**: 1,000 nonparametric replicates; support values mapped onto the best-scoring ML tree (`.raxml.support` file used downstream)

---

## Step 5 — Transmission Cluster Identification

### ClusterPicker

Run once per region × bootstrap threshold × genetic distance threshold combination. Bootstrap thresholds evaluated: 70, 90, 99. Genetic distance thresholds evaluated: 0.5%–4.5% in 0.5% steps.

| Argument | Value | Description |
|----------|-------|-------------|
| 1 | `DECIPHER_aligned_Env_sequences.fasta` | DECIPHER-aligned FASTA file |
| 2 | `Env_mytree_bs.raxml.support` | RAxML-NG bootstrap support tree |
| 3 | `90` | Initial bootstrap support threshold (%) for cluster detection |
| 4 | `90` | Final bootstrap support threshold (%) — set equal to argument 3 |
| 5 | `0.03` | Maximum pairwise genetic distance within a cluster (3%) |
| 6 | `2` | Minimum cluster size (at least 2 sequences) |
| 7 | `gap` | Distance method: gap positions are ignored in distance calculation |

```bash
java -jar ClusterPicker_1.2.5.jar \
    DECIPHER_aligned_Env_sequences.fasta \
    Env_mytree_bs.raxml.support \
    90 90 \
    0.03 \
    2 \
    gap
```

Repeat for each region, adjusting the FASTA file, tree file, bootstrap threshold (90 → 70 or 99), and genetic distance (0.03 → 0.005–0.045) accordingly.

### HIV-TRACE

Run once per region × genetic distance threshold combination. Genetic distance thresholds evaluated: 0.5%–4.5% in 0.5% steps. Use the appropriate HXB2 reference per region (`HXB2_pol` for gag-pol, protease, RT, integrase, and full-length; `HXB2_env` for env; `HXB2_vif` for vif; `HXB2_nef` for nef).

| Flag | Value | Description |
|------|-------|-------------|
| `-i` | `DECIPHER_aligned_Env_sequences.fasta` | DECIPHER-aligned FASTA file |
| `-a` | `0.015` | Ambiguity threshold: sequences with >1.5% ambiguous bases are excluded |
| `-r` | `HXB2_env` | HXB2 reference name for the genomic region |
| `-t` | `0.03` | Genetic distance threshold for cluster membership (3%) |
| `-m` | `200` | Minimum nucleotide overlap required between two sequences |
| `-g` | `0.9` | Minimum overlap fraction (90% of the shorter sequence must overlap) |
| `-o` | `Env_hivtrace_GD030.json` | Output JSON file |

```bash
hivtrace \
    -i DECIPHER_aligned_Env_sequences.fasta \
    -a 0.015 \
    -r HXB2_env \
    -t 0.03 \
    -m 200 \
    -g 0.9 \
    -o Env_hivtrace_GD030.json
```

Repeat for each region, adjusting `-i`, `-r`, `-t`, and `-o` accordingly.

---

## Step 6 — Cluster Concordance Analysis

`compare_cluster.R` compares subregion clusters (ClusterPicker and HIV-TRACE, swept across bootstrap and genetic-distance thresholds) against the full-length genome reference clustering (ClusterPicker, GD = 3.0%, BS = 90%). It reports:

| Metric | Meaning |
|--------|---------|
| `V_measure_refclustered` | **Headline metric.** V-measure scored only over sequences whose full-length reference cluster has ≥ 2 members, so the large mutually-singleton background typical of HIV clustering doesn't saturate the score near 1. |
| `pw_recall` | Pairwise recall: fraction of full-length-linked pairs the subregion also links (structure recovered). |
| `pw_precision` | Pairwise precision: fraction of subregion-linked pairs also linked in the full-length reference (spurious-link rate). |
| `pw_f1` | Harmonic mean of `pw_precision` and `pw_recall`. |

Edit the path variables at the top of `compare_cluster.R`, then:

```bash
Rscript compare_cluster.R regions.txt
```

Output: `vmeasure_results.tsv` — one row per region × method × GD threshold × bootstrap combination.

---

## File Structure

```
HIVphyloSeq_V2/
├── README.md                     # this file
├── sequence_selection.md         # LANL download instructions, filtering criteria, cohort description
├── compare_cluster.R             # multi-metric cluster concordance analysis (clevr R package)
└── accession_numbers.txt         # GenBank accessions for the 496 sequences in the final analysis set
```
