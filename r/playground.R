# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary packages
BiocManager::install(c("GEOquery", "limma", "biomaRt", "AnnotationDbi", "org.Hs.eg.db"))

# Load necessary libraries
library(GEOquery)
library(limma)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)

dataset_dir <- "data/datasets"
geo_datasets <- list.files(dataset_dir, pattern = "GSE", full.names = TRUE)

gse <- getGEO(filename = geo_datasets[3])
platform_id <- annotation(gse)

expr_matrix <- exprs(gse)

expr_matrix_t <- t(expr_matrix)

expr_df <- as.data.frame(expr_matrix_t)

gene_symbols <- fData(gse)$GENE_SYMBOL
colnames(expr_df) <- gene_symbols
empty_columns <- which(colnames(expr_df) == "")
expr_df <- expr_df[ , -empty_columns]

expr_df
