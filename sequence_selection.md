# Sequence Selection

Sequences were obtained from the [Los Alamos National Laboratory (LANL) HIV Sequence Database](https://www.hiv.lanl.gov/), downloaded Feb. 4, 2025 via the [Advanced Search Interface](https://www.hiv.lanl.gov/components/sequence/HIV/asearch/map_db.comp). The filtering pipeline below reduced the full LANL dataset to the final analysis cohort of 496 full-length HIV-1 subtype B genomes (one per patient) sampled 2005–2024 in the US, Canada, and Europe.

## Data Download

LANL's search interface caps the number of records returned per query, so the full database was retrieved in several sequential batches by Se ID range, each saved as a tab-delimited `.txt` file in one folder.

### Step-by-step

1. Go to the [HIV Sequence Database](https://www.hiv.lanl.gov/), click **SEARCH** in the top menu, then **Advanced Search** (or go directly to the [Advanced Search Interface](https://www.hiv.lanl.gov/components/sequence/HIV/asearch/map_db.comp)). This is a different tool from the default "Search DB" page — it lets you choose exactly which database fields to return, which the plain search does not.
2. The page shows one box per LANL database table (sequence data, genome coordinates, publication linkage, cluster annotation, patient metadata, sample metadata, accession numbers, etc.), each with individual field checkboxes. Select the fields needed to reproduce the filtering criteria below: sequence and length, HXB2 genome-coordinate start/stop, subtype, country, sampling year, QC flag, patient-linkage ID, days-from-sample fields (for tie-breaking), and GenBank accession. Free-text fields not used by the filtering criteria (comments, journal titles, clinical fields like viral load or CD4/CD8 counts, etc.) aren't needed.
3. Leave all search-criteria fields blank so the query returns every record in the database — apply no subtype/country/year restriction at download time; all filtering happens afterward, per the criteria below.
4. Submit the search, then export/download the results as tab-delimited text. Because LANL caps records per query, repeat the search restricting by `Se ID` range (or another field) to pull the full database in batches, saving each batch as its own `.txt` file in one folder for the filtering scripts to read.

Each downloaded file has two header lines (`Number of records retrieved: N`, then a blank line) before the tab-delimited data, which the filtering scripts skip.

## Filtering Steps

| Step | Criterion | Remaining (n) | Excluded (n) |
|------|-----------|---------------|--------------|
| Full LANL dataset | — | 2,727,128 | — |
| Subtype filter | HIV-1 subtype B only | 1,361,074 | 1,366,054 |
| Geographic filter | US, Canada, or one of 19 European countries | 952,093 | 408,981 |
| Temporal filter | Sampling year 2005–2024 | 450,151 | 501,942 |
| Quality control | Passed LANL QC (`Problematic Sequence == "0"`) | 409,263 | 40,888 |
| Length filter | Sequence length > 7,000 bp | 23,670 | 385,593 |
| Sequence-ID dedup | One row per `Se ID` (removes duplicate publication rows) | 11,520 | 12,150 |
| Patient linkage | Sequence linked to a patient ID (`PAT id(SSAM)` non-missing) | 11,334 | 186 |
| One per patient | Earliest sample retained; ties broken by lowest `Days from first Sample`, then lowest `Accession` | 610 | 10,724 |
| Complete gag coverage | Sequences fully spanning the gag region (HXB2 790–2292) | 496 | 114 |

**Final analysis set: n = 496**

### Geographic filter — 21 countries

United States, Canada, Portugal, Spain, France, Belgium, Luxembourg, Switzerland, Italy, Germany, Austria, Czech Republic, Slovenia, Denmark, Poland, Slovakia, Hungary, Croatia, Sweden, Norway, United Kingdom.

## Rationale for Criteria

- **Subtype B**: Restricts analysis to a single subtype to reduce confounding in phylogenetic inference; subtype B predominates in North America and Western/Central Europe.
- **US, Canada & Europe**: A high-income-country sampling frame restricted to countries where subtype B circulates as the dominant clade.
- **2005–2024**: Two decades of sampling, giving dense temporal coverage while still relying on molecular-clock-informative sampling.
- **LANL QC**: Removes sequences flagged for excessive ambiguous nucleotides, frameshifts, hypermutation, or poor HXB2 alignment — artifacts that distort phylogenetic inference and genetic distance estimation.
- **Length > 7,000 bp**: Retains near-complete genomes suitable for full-genome analysis.
- **Sequence-ID dedup**: A single sequence can appear on multiple rows because LANL emits one row per linked publication; this collapses those to one row per `Se ID`.
- **Patient linkage**: Sequences that cannot be tied to a patient ID cannot be deduplicated at the patient level and are excluded.
- **One per patient**: Prevents over-representation of individuals with multiple longitudinal sequences; the earliest available sample is kept to approximate incident/baseline virus.
- **Complete gag coverage**: Required for consistent subregion extraction across all downstream analyses.

## Subregions Analyzed

Seven commonly deposited HIV-1 genomic subregions were extracted from the final 496 sequences using their mapped HXB2 coordinates:

| Subregion | HXB2 Start | HXB2 End | Length (bp) | Final N |
|-----------|-----------|---------|-------------|---------|
| Protease | 2,253 | 2,549 | 297 | 496 |
| Reverse Transcriptase (p51/RT) | 2,550 | 3,869 | 1,320 | 496 |
| Integrase (p31) | 4,230 | 5,096 | 867 | 496 |
| Vif | 5,041 | 5,619 | 579 | 496 |
| Nef | 8,797 | 9,417 | 621 | 463 |
| Env | 6,225 | 8,795 | 2,571 | 491 |
| Gag-Pol | 790 | 5,096 | 4,307 | 496 |

These regions correspond to genomic portions frequently sequenced for drug resistance testing, vaccine studies, and immune-escape analyses, and represent the most commonly deposited subregions in the LANL HIV Sequence Database. A sequence is retained per subregion only if it fully spans that subregion's HXB2 coordinates; regions with incomplete coverage in a given genome are dropped from that subregion's alignment, which is why Env and Nef fall below the full N = 496.
