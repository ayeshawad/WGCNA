########################################
# WGCNA Analysis for Crohn's Disease
# 02: WGCNA Network Construction
########################################

source("00_config.R")
source("01_data_preprocessing.R")

# Load required libraries
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cat("\n=== Step 1: Prepare expression matrix ===\n")

# Transpose: genes x samples -> samples x genes (WGCNA requirement)
expr_cd <- t(vst_mat_adjusted)

# Remove zero-variance genes
gene_var <- apply(expr_cd, 2, var, na.rm = TRUE)
expr_cd <- expr_cd[, gene_var > 0, drop = FALSE]

cat("Expression matrix dimensions (samples x genes):", dim(expr_cd), "\n")

cat("\n=== Step 2: Quality checks ===\n")

# WGCNA quality check
gsg <- goodSamplesGenes(expr_cd, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing", sum(!gsg$goodGenes), "bad genes\n")
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing", sum(!gsg$goodSamples), "bad samples\n")
  }
  expr_cd <- expr_cd[gsg$goodSamples, gsg$goodGenes]
}

# Align metadata
meta_cd <- coldata[rownames(coldata) %in% rownames(expr_cd), , drop = FALSE]
common_samples <- intersect(rownames(expr_cd), rownames(meta_cd))
expr_cd <- expr_cd[common_samples, , drop = FALSE]
meta_cd <- meta_cd[common_samples, , drop = FALSE]

cat("\n=== Step 3: Sample outlier detection ===\n")

# Create sample dendrogram
sampleTree <- hclust(dist(expr_cd), method = "average")
pdf(file.path(results_dir, "sample_tree_raw.pdf"),
    width = 18, height = 5)
plot(sampleTree, main = "CD samples clustering", xlab = "", sub = "")
dev.off()

# Calculate sample connectivity
sampleAdj <- adjacency(t(expr_cd), type = "signed")
sampleConn <- rowSums(sampleAdj) - 1
Z.k <- scale(sampleConn)
Z.k <- as.numeric(Z.k)
names(Z.k) <- rownames(expr_cd)

cat("Sample connectivity Z-scores:\n")
print(summary(Z.k))

# Visualize connectivity on dendrogram
zColors <- numbers2colors(Z.k, signed = TRUE)
pdf(file.path(results_dir, "sample_tree_Zk.pdf"),
    width = 9, height = 6)
plotDendroAndColors(
  sampleTree,
  colors = zColors,
  groupLabels = "Z.k (standardized connectivity)",
  main = "CD sample dendrogram with connectivity"
)
dev.off()

# Identify outliers
outlierSamples <- Z.k < outlier_zcut
if (sum(outlierSamples) > 0) {
  cat("Removing", sum(outlierSamples), "outlier samples (Z.k <", outlier_zcut, ")\n")
  outlier_names <- names(Z.k)[outlierSamples]
  cat("Outliers:", paste(outlier_names, collapse = ", "), "\n")
  
  # Visualize outliers
  sampleColor_out <- ifelse(outlierSamples, "red", "grey")
  pdf(file.path(results_dir, "sample_tree_outlier_flag.pdf"),
      width = 15, height = 6)
  plotDendroAndColors(
    sampleTree,
    colors = sampleColor_out,
    groupLabels = "Outlier by Z.k",
    main = paste0("CD sample dendrogram outliers (Z.k <", outlier_zcut, ")")
  )
  dev.off()
  
  expr_cd <- expr_cd[!outlierSamples, , drop = FALSE]
  meta_cd <- meta_cd[rownames(expr_cd), , drop = FALSE]
  
  # Final sample dendrogram after outlier removal
  sampleTree2 <- hclust(dist(expr_cd), method = "average")
  pdf(file.path(results_dir, "sample_tree_after_outliers.pdf"),
      width = 8, height = 5)
  plot(sampleTree2, main = "CD samples after outlier removal", xlab = "", sub = "")
  dev.off()
} else {
  cat("No outliers detected\n")
}

cat("Final sample count:", nrow(expr_cd), "\n")
cat("Final gene count:", ncol(expr_cd), "\n")

cat("\n=== Step 4: Soft threshold selection ===\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(
  expr_cd,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# Visualize soft threshold selection
pdf(file.path(results_dir, "soft_threshold_selection.pdf"),
    width = 12, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")
dev.off()

cat("Selected soft power:", soft_power, "\n")
cat("R^2 at selected power:", 
    sft$fitIndices[sft$fitIndices$Power == soft_power, "SFT.R.sq"], "\n")

cat("\n=== Step 5: Module detection ===\n")

net <- blockwiseModules(
  expr_cd,
  power              = soft_power,
  networkType        = "signed",
  TOMType            = "signed",
  corType            = "bicor",
  maxPOutliers       = 0.1,
  minModuleSize      = min_module_size,
  deepSplit          = 4,
  reassignThreshold  = 0,
  mergeCutHeight     = 0,
  numericLabels      = TRUE,
  pamRespectsDendro  = FALSE,
  saveTOMs           = FALSE,
  verbose            = 3
)

# Extract module assignments
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)

cat("\nModule summary:\n")
print(table(moduleColors))

# Save gene dendrogram
pdf(file.path(results_dir, "gene_tree_with_modules.pdf"),
    width = 10, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

cat("\n=== WGCNA network construction complete ===\n")
cat("Total modules (excluding grey):", sum(moduleColors != "grey"), "\n")
cat("Output: net (network object)\n")
cat("Output: moduleLabels (numeric labels)\n")
cat("Output: moduleColors (color labels)\n")
cat("Output: expr_cd (final expression matrix)\n")
cat("Output: meta_cd (final metadata)\n")

