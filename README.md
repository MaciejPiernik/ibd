# IBD Meta-Analysis

Code for the paper: **"A Multi-Level Meta-Analysis of Transcriptomic Data Reveals Constitutive and Inflammation-Dependent Dysregulation in Inflammatory Bowel Disease"**

This repository implements the complete analytical pipeline: from downloading and processing GEO microarray datasets, through differential expression and pathway enrichment, to random-effects meta-analysis, pathway clustering, cross-comparison synthesis, and generation of all paper figures and tables.

## Overview

The analysis covers **14 GEO datasets** (1,470 samples, 9 microarray platforms) and runs **four parallel meta-analyses**:

| Comparison | Datasets | Description |
|---|---|---|
| UC inflamed vs healthy controls | 12 | Primary analysis |
| UC uninflamed vs healthy controls | 5 | Constitutive dysregulation |
| CD inflamed vs healthy controls | 8 | Cross-disease comparison |
| UC vs CD (direct) | 6 | Disease-specificity analysis |

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

The data processing pipeline downloads datasets from GEO, maps probes to genes, and harmonizes metadata. Dataset definitions are in `data/DAQ.csv`.

```bash
# Process all datasets
python -m ibd run-pipeline

# Process specific datasets
python -m ibd run-pipeline -d GSE11223,GSE75214,GSE6731

# Force recomputation
python -m ibd run-pipeline -r   # recompute expression data
python -m ibd run-pipeline -m   # recompute metadata
```

Processed data is saved to the `db/` directory as parquet (expression) and CSV (metadata) files.

### Step 2: Run the meta-analysis experiment

Each of the four comparisons is run separately by modifying the group constants in `ibd/experiments/experiment.py` (`GROUP_1`, `GROUP_2`, `INFLAMMATION_1`, `INFLAMMATION_2`) and specifying an output directory.

```bash
python -m ibd run-experiment \
    --data-dir ./db \
    --output-dir results/paper/uc_vs_hc
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

These scripts process results from all four comparisons (in `results/paper/{uc_vs_hc,uninf_vs_hc,cd_vs_hc,uc_vs_cd}/`) and produce:
- `pathway_meta_significant_knee.csv` -- pathways passing the knee-point filter
- `hierarchical_pathway_clusters.csv` -- pathway community assignments
- `pathway_meta_significant_knee_with_clusters.csv` -- merged result

### Step 4: Cross-comparison analysis

```bash
# Inflamed vs uninflamed effect size comparison (Wald test, constitutive/inflammation-dependent classification)
python ibd/scripts/S_effect_size_comparison.py \
    --inflamed results/paper/uc_vs_hc/gene_meta_all.csv \
    --uninflamed results/paper/uninf_vs_hc/gene_meta_all.csv \
    --output-dir results/paper/

# Compare significant pathways between inflamed and uninflamed analyses
python ibd/scripts/S_compare_pathways.py
```

### Step 5: Generate figures and tables

Scripts are numbered in order. Prefixes indicate output type: `fig` = figure, `tab` = table, `S_` = supplementary.

```bash
# Figures
python ibd/scripts/00_fig_generate_methods_figure.py          # Methods flowchart (SVG)
python ibd/scripts/03_fig_volcano_plot.py                      # Volcano plot
python ibd/scripts/06_fig_scatter_comparison.py --input results/paper/comparison_stats.csv  # Constitutive gene scatter
python ibd/scripts/07_fig_pathway_heatmap.py                   # Three-way pathway heatmap
python ibd/scripts/08_fig_uc_cd_concordance.py                 # UC-CD concordance + pathway clusters

# Tables (LaTeX)
python ibd/scripts/01_tab_dataset_summary.py                   # Dataset summary
python ibd/scripts/02_tab_pipeline_summary.py                  # Pipeline summary across 4 analyses
python ibd/scripts/04_tab_top_genes.py                         # Top up/downregulated genes
python ibd/scripts/05_tab_pathway_clusters.py                  # Pathway cluster tables
python ibd/scripts/09_tab_cross_comparison_table.py            # Cross-comparison gene table

# Supplementary
python ibd/scripts/S_generate_figures.py                       # Forest plots, I² diagnostics, Stouffer's check
python ibd/scripts/S_generate_supplementary.py --data-dir results/paper/  # Full supplementary PDF + CSV tables
```

### Utility scripts

```bash
# Build pathway overlap network graph (JSON for visualization)
python ibd/scripts/_build_pathway_graph.py --pathways <path> --genes <path>

# Build Reactome hierarchy tree with enrichment overlay (JSON)
python ibd/scripts/_build_reactome_hierarchy.py --pathways <path>
```

## Project structure

```
ibd/
├── __main__.py                     # CLI entry point (run-pipeline, run-experiment)
├── core/
│   ├── data/
│   │   ├── processing/             # Metadata parsers, normalization
│   │   ├── structures/Dataset.py   # Dataset class (download, process, store)
│   │   └── utils.py                # Data loading utilities
│   ├── platforms/                   # Platform-specific probe-to-gene mapping (GPL570, GPL6244, etc.)
│   ├── ensembl/                     # Ensembl gene mapping
│   └── utils/genes.py              # Entrez ID to gene symbol conversion
├── experiments/
│   ├── experiment.py               # Full experiment pipeline (DE, GSEA, meta-analysis)
│   └── Config.py                   # Configuration dataclass
└── scripts/                        # Post-processing, figures, and tables (see above)

data/
├── DAQ.csv                         # Dataset acquisition table (accession numbers, platforms)
└── geo_cache/                      # Cached GEO downloads

db/                                 # Processed expression + metadata files

results/paper/                      # Experiment outputs
├── uc_vs_hc/                       # UC inflamed vs healthy controls
├── uninf_vs_hc/                    # UC uninflamed vs healthy controls
├── cd_vs_hc/                       # CD inflamed vs healthy controls
└── uc_vs_cd/                       # UC vs CD direct comparison
```
