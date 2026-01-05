########################################
# WGCNA Analysis for Crohn's Disease
# 01: Data Preprocessing
########################################

source("00_config.R")

# Load required libraries
required_packages <- c(
  "DESeq2", "WGCNA", "matrixStats", "dplyr", "AnnotationDbi", "org.Hs.eg.db"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed"))
  }
  library(pkg, character.only = TRUE)
}

cat("\n=== Step 1: Loading data ===\n")

# Load count matrix
cts <- read.table(
  count_matrix_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
cat("Count matrix dimensions:", dim(cts), "\n")

# Load metadata
coldata <- read.csv(metadata_path, check.names = FALSE)
rownames(coldata) <- coldata$sample

# Filter to CD samples only
coldata <- coldata[coldata$disease == "CD", , drop = FALSE]
cat("CD samples:", nrow(coldata), "\n")

# Align count matrix and metadata
common_samples <- intersect(colnames(cts), rownames(coldata))
cts <- cts[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

# Filter low-expression genes (keep genes with >=5 counts in >=25% of samples)
keep_genes <- rowSums(cts >= 5) >= (0.25 * ncol(cts))
cts <- cts[keep_genes, , drop = FALSE]
cat("Genes after filtering:", nrow(cts), "\n")

# Clean gene IDs (remove version numbers)
rownames(cts) <- gsub("\\..*", "", rownames(cts))
rownames(cts) <- substr(rownames(cts), 1, 18)

cat("\n=== Step 2: VST normalization ===\n")

# Prepare data for DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(cts),
  colData = coldata,
  design = ~1
)

# Additional filter: remove genes with <5 total counts
dds <- dds[rowSums(counts(dds)) > 5, ]

# Variance stabilizing transformation (blind = TRUE for unsupervised)
vst_dds <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_dds)

cat("VST matrix dimensions:", dim(vst_mat), "\n")

cat("\n=== Step 3: Covariate adjustment ===\n")

# Prepare covariates
coldata$batch_detail <- as.factor(coldata$batch_detail)
coldata$sex <- as.factor(coldata$sex)
coldata$TIN_median <- as.numeric(coldata$TIN_median)

covariates <- data.frame(
  TIN = coldata$TIN_median,
  batch = as.numeric(coldata$batch_detail),
  sex = as.numeric(coldata$sex),
  row.names = rownames(coldata)
)

# Ensure samples match
covariates <- covariates[colnames(vst_mat), , drop = FALSE]

# Remove covariate effects using empiricalBayesLM
vst_mat_adjusted <- empiricalBayesLM(
  data = t(vst_mat),
  removedCovariates = as.matrix(covariates)
)$adjustedData

# Transpose back to genes x samples
vst_mat_adjusted <- t(vst_mat_adjusted)

cat("Adjusted matrix dimensions:", dim(vst_mat_adjusted), "\n")

cat("\n=== Step 4: Variance-based gene filtering ===\n")

# Calculate variance for each gene
gene_variance <- matrixStats::rowVars(vst_mat_adjusted)

# Select top N most variable genes
top_var_genes_idx <- names(sort(gene_variance, decreasing = TRUE))[1:min(top_var_genes, length(gene_variance))]
vst_mat_adjusted <- vst_mat_adjusted[top_var_genes_idx, , drop = FALSE]

cat("Final gene count:", nrow(vst_mat_adjusted), "\n")
cat("Final sample count:", ncol(vst_mat_adjusted), "\n")

cat("\n=== Data preprocessing complete ===\n")
cat("Output: vst_mat_adjusted (genes x samples)\n")
cat("Output: coldata (metadata)\n")

