# Ulcerative Colitis Meta-Analysis

Code for the paper: **"Large-scale transcriptomic meta-analysis identifies constitutive and inflammation-dependent layers in ulcerative colitis"**

This repository implements the complete analytical pipeline: from downloading and processing GEO microarray datasets, through differential expression and pathway enrichment, to random-effects meta-analysis, pathway clustering, cross-comparison synthesis, and generation of all paper figures and tables.

## Overview

The discovery analysis covers **17 GEO microarray datasets** (1,677 samples, 10 array
platforms) and runs **five parallel meta-analyses**. An independent **9-dataset RNA-seq
corpus** (1,250 samples) validates them through the same pipeline (see Step 6).

| Comparison | Datasets | Description |
|---|---|---|
| UC inflamed vs healthy controls | 15 | Primary analysis |
| UC uninflamed vs healthy controls | 6 | Constitutive dysregulation |
| UC inflamed vs UC uninflamed | 6 | Direct inflammation contrast |
| CD inflamed vs healthy controls | 9 | Cross-disease comparison |
| UC vs CD (direct) | 7 | Disease-specificity analysis |

Each comparison goes through the same pipeline: per-dataset limma differential expression, Hedges' g effect sizes, GSEA (Reactome 2022), REML random-effects meta-analysis at both gene and pathway levels, knee-point filtering, and Louvain community detection on pathway overlap networks.

## Setup

Install the package:
```bash
pip install .
```

### R dependency

The experiment pipeline uses **limma** (via `rpy2`) for moderated t-statistics. You need R installed with the limma package:
```r
install.packages("BiocManager")
BiocManager::install("limma")
```

## Usage

### Step 1: Download and process data

The data processing pipeline downloads datasets from GEO, maps probes or transcript
identifiers to genes, and harmonizes metadata.

Two corpora are processed, selected with `--source`. The microarray corpus is the
discovery set; the RNA-seq corpus is the independent validation set (see Step 6).

| `--source` | Dataset definitions | Written to |
|---|---|---|
| `microarray` (default) | `data/DAQ.csv` | `db/` |
| `rnaseq` | `data/DAQ_rnaseq.csv` | `db_rnaseq/` |

```bash
# Process all datasets of a corpus
python -m ibd run-pipeline
python -m ibd run-pipeline --source rnaseq

# Process specific datasets
python -m ibd run-pipeline -d GSE11223,GSE75214,GSE6731

# Force recomputation
python -m ibd run-pipeline -r   # recompute expression data
python -m ibd run-pipeline -m   # recompute metadata
```

Both corpora are saved as parquet (expression, samples x genes on a log2 scale) and
CSV (metadata) files. Datasets listed in an acquisition table but without a
processor are reported and skipped, so the tables can also record the candidates
that were screened out and why.

### Step 2: Run the meta-analysis experiment

Each comparison is run separately by editing two things in
`ibd/experiments/Experiment.py` and then running the command below with the
matching `--source` and output directory:

1. the group constants `GROUP_1`, `INFLAMMATION_1`, `GROUP_2`, `INFLAMMATION_2`
2. `min_datasets` in `ExperimentConfig`, set to 50% of the contributing cohorts
   (rounded up), which is the presence threshold used throughout the paper

`--source` selects the corpus: it points the run at `db/` or `db_rnaseq/` and, for
RNA-seq, skips the Entrez-to-symbol mapping, since those matrices are already
keyed by gene symbol.

| Output directory | `--source` | GROUP_1 / INFLAMMATION_1 | GROUP_2 / INFLAMMATION_2 | Cohorts | `min_datasets` |
|---|---|---|---|---|---|
| `results/revision_v2/uc_vs_hc` | microarray | UC / Inflamed | Ctrl / Uninflamed | 15 | 8 |
| `results/revision_v2/uninf_vs_hc` | microarray | UC / Uninflamed | Ctrl / Uninflamed | 6 | 3 |
| `results/revision_v2/cd_vs_hc` | microarray | CD / Inflamed | Ctrl / Uninflamed | 9 | 5 |
| `results/revision_v2/uc_vs_cd` | microarray | UC / Inflamed | CD / Inflamed | 7 | 4 |
| `results/revision_v2/infl_vs_uninfl` | microarray | UC / Inflamed | UC / Uninflamed | 6 | 3 |
| `results/revision_v2/rnaseq_uc_vs_hc` | rnaseq | UC / Inflamed | Ctrl / Uninflamed | 8 | 4 |
| `results/revision_v2/rnaseq_uninf_vs_hc` | rnaseq | UC / Uninflamed | Ctrl / Uninflamed | 3 | 2 |
| `results/revision_v2/rnaseq_infl_vs_uninfl` | rnaseq | UC / Inflamed | UC / Uninflamed | 2 | 1 |

When both arms come from the same disease group, as in the last two rows, the
arms are labelled by disease and inflammation status together so that they can be
told apart.

```bash
python -m ibd run-experiment \
    --output-dir results/revision_v2/uc_vs_hc

python -m ibd run-experiment --source rnaseq \
    --output-dir results/revision_v2/rnaseq_uc_vs_hc
```

The experiment produces:
- `per_dataset_gene_stats.csv` -- per-dataset gene-level statistics (Hedges' g, limma t-stats)
- `dataset_summary.csv` -- sample counts per dataset
- `per_dataset_gsea.csv` -- per-dataset GSEA results
- `gene_meta_all.csv` / `gene_meta_significant.csv` -- gene-level meta-analysis
- `gene_rra.csv` -- Robust Rank Aggregation (validation)
- `pathway_meta_all.csv` / `pathway_meta_significant.csv` -- pathway-level meta-analysis
- `pathway_leading_edge_detail.csv` -- leading-edge gene detail

Key options:
```
--alpha             FDR threshold (default: 0.05)
--min-samples       Minimum samples in a dataset (default: 6)
--min-per-group     Minimum samples per group (default: 3)
--gsea-permutations Number of GSEA permutations (default: 1000)
--geneset-libraries Comma-separated GSEA libraries (default: Reactome_2022)
--skip-gene-mapping Skip Entrez-to-symbol mapping
```

### Step 3: Pathway filtering and clustering

After running all four comparisons, apply knee-point filtering and hierarchical community detection:

```bash
# Knee-point filtering on sorted |NES| to select top pathways
python ibd/scripts/_cut_pathways_at_knee.py

# Hierarchical Louvain community detection on leading-edge gene overlap networks
# Edit the DIR variable for each comparison directory
python ibd/scripts/_hierarchical_community_detection.py
```

These scripts process results from all five microarray comparisons (in
`results/revision_v2/{uc_vs_hc,uninf_vs_hc,cd_vs_hc,uc_vs_cd,infl_vs_uninfl}/`) and produce:
- `pathway_meta_significant_knee.csv` -- pathways passing the knee-point filter
- `hierarchical_pathway_clusters.csv` -- pathway community assignments
- `pathway_meta_significant_knee_with_clusters.csv` -- merged result

### Step 4: Cross-comparison analysis

```bash
# Inflamed vs uninflamed effect size comparison (Wald test, constitutive/inflammation-dependent classification)
python ibd/scripts/S_effect_size_comparison.py \
    --inflamed results/revision_v2/uc_vs_hc/gene_meta_all.csv \
    --uninflamed results/revision_v2/uninf_vs_hc/gene_meta_all.csv \
    --output-dir results/revision_v2/

# Compare significant pathways between inflamed and uninflamed analyses
python ibd/scripts/S_compare_pathways.py

# Classify genes into activity layers by triage between the three UC contrasts
python ibd/scripts/S_gene_layers.py
```

### Step 5: Generate figures and tables

Scripts are numbered in order. Prefixes indicate output type: `fig` = figure, `tab` = table, `S_` = supplementary.

```bash
# Figures
python ibd/scripts/00_fig_generate_methods_figure.py          # Methods flowchart (SVG)
python ibd/scripts/03_fig_volcano_plot.py                      # Volcano plot
python ibd/scripts/06_fig_scatter_comparison.py                # Activity-layer scatter (constitutive vs amplified)
python ibd/scripts/07_fig_pathway_heatmap.py                   # Three-way pathway heatmap
python ibd/scripts/08_fig_uc_cd_concordance.py                 # UC-CD concordance + pathway clusters

# Tables (LaTeX)
python ibd/scripts/01_tab_dataset_summary.py                   # Dataset summary
python ibd/scripts/02_tab_pipeline_summary.py                  # Pipeline summary across 4 analyses
python ibd/scripts/04_tab_top_genes.py                         # Top up/downregulated genes
python ibd/scripts/05_tab_pathway_clusters.py                  # Pathway cluster tables
python ibd/scripts/09_tab_cross_comparison_table.py            # Cross-comparison gene table

# Supplementary
python ibd/scripts/S_generate_supplementary.py --data-dir results/revision_v2/  # Full supplementary PDF + CSV tables
```

### Step 6: RNA-seq validation

The microarray results are validated against independent RNA-seq cohorts, using
the same pipeline throughout: only `--source rnaseq` differs.

The validation set is every RNA-seq dataset in GEO satisfying the same inclusion
criteria as the discovery analysis, so that it is exhaustive rather than a
selection. `_screen_rnaseq.py` reproduces that screen from the GEO search onwards
and writes `results/revision_v2/rnaseq_screening.csv` (258 candidates, 9
included). Its outcome is already recorded in `data/DAQ_rnaseq.csv`, so re-running
it is only necessary to re-derive the screen.

```bash
# Optional: re-run the GEO screen that defines the validation set
python ibd/scripts/_screen_rnaseq.py

# 1. Build the nine cohorts into db_rnaseq/
python -m ibd run-pipeline --source rnaseq

# 2. Run the two comparisons (see the table in Step 2 for the constants each needs)
python -m ibd run-experiment --source rnaseq --output-dir results/revision_v2/rnaseq_uc_vs_hc
python -m ibd run-experiment --source rnaseq --output-dir results/revision_v2/rnaseq_uninf_vs_hc

# 3. Compare the platforms and report effective cohort counts
python ibd/scripts/S_compare_platforms.py
```

Step 1 downloads roughly 0.5 GB from GEO into `db_rnaseq_raw/` and caches it, so
re-running is cheap. Step 2 is the slow part; GSEA over eight cohorts takes a
while. Step 3 writes gene- and pathway-level concordance tables to
`results/revision_v2/comparison/`, together with the effective number of cohorts
contributing to each meta-analysis, reported for the discovery analyses as well.

### Utility scripts

```bash
# Build pathway overlap network graph (JSON for visualization)
python ibd/scripts/_build_pathway_graph.py --pathways <path> --genes <path>

# Build Reactome hierarchy tree with enrichment overlay (JSON)
python ibd/scripts/_build_reactome_hierarchy.py --pathways <path>
```

## Project structure

The two corpora are processed by parallel classes: `Dataset` reads microarray
series through a per-platform probe-to-gene mapper and a per-dataset metadata
processor, while `RnaSeqDataset` reads RNA-seq series through a per-dataset
processor that knows how that study deposited its data. Both write the same
parquet + CSV pair, so everything downstream is shared.

```
ibd/
├── __main__.py                     # CLI entry point (run-pipeline, run-experiment)
├── core/
│   ├── data/
│   │   ├── processing/
│   │   │   ├── metadata_processors.py   # Per-dataset metadata parsers (microarray)
│   │   │   └── rnaseq_processors.py     # Per-dataset expression + metadata parsers (RNA-seq)
│   │   ├── structures/
│   │   │   ├── Dataset.py               # Microarray series (download, process, store)
│   │   │   └── RnaSeqDataset.py         # RNA-seq series (download, process, store)
│   │   └── utils.py                # Data loading utilities
│   ├── platforms/                   # Platform-specific probe-to-gene mapping (GPL570, GPL6244, etc.)
│   ├── ensembl/                     # Ensembl gene mapping
│   └── utils/genes.py              # Entrez ID to gene symbol conversion
├── experiments/
│   ├── Experiment.py               # Full experiment pipeline (DE, GSEA, meta-analysis)
│   └── Config.py                   # Configuration dataclass
└── scripts/                        # Post-processing, figures, and tables (see above)

data/
├── DAQ.csv                         # Microarray acquisition table (accessions, platforms)
└── DAQ_rnaseq.csv                  # RNA-seq acquisition table, with screening outcomes

db/                                 # Processed microarray expression + metadata
db_rnaseq/                          # Processed RNA-seq expression + metadata

results/revision_v2/                # Every result reported in the paper
├── uc_vs_hc/                       # UC inflamed vs healthy controls
├── uninf_vs_hc/                    # UC uninflamed vs healthy controls
├── cd_vs_hc/                       # CD inflamed vs healthy controls
├── uc_vs_cd/                       # UC vs CD direct comparison
├── infl_vs_uninfl/                 # UC inflamed vs UC uninflamed
├── rnaseq_uc_vs_hc/                # UC inflamed vs healthy controls, RNA-seq
├── rnaseq_uninf_vs_hc/             # UC uninflamed vs healthy controls, RNA-seq
├── rnaseq_infl_vs_uninfl/          # UC inflamed vs uninflamed, RNA-seq
└── comparison/                      # Cross-platform concordance

results/paper/                      # Superseded first-submission run, kept for reference
```
