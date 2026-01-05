########################################
# WGCNA Analysis for Crohn's Disease
# Configuration File
########################################

# Set working directory (adjust as needed)
# The scripts assume they are run from the project root directory
# setwd("/path/to/WGCNA_CD_Analysis")

# Create results directories
results_dir <- "results"
fig_dir <- file.path(results_dir, "figures")
go_dir <- file.path(results_dir, "GO")
hub_dir <- file.path(results_dir, "hubs")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(hub_dir, showWarnings = FALSE, recursive = TRUE)

# Input file paths (relative to project root)
count_matrix_path <- "input/count_matrix.txt"
metadata_path <- "input/combined.csv"

# Analysis parameters
top_var_genes <- 2500          # Number of most variable genes to retain
soft_power <- 10               # Soft thresholding power (adjust based on scale-free topology)
min_module_size <- 30          # Minimum module size
outlier_zcut <- -2             # Z-score cutoff for sample outlier removal

# Module focus for detailed analysis
# Selected based on: (1) strongest associations, (2) biological interpretability,
# (3) combinatorial predictive power
key_modules <- c("Module_3", "Module_6", "Module_13", "Module_14", "Module_8")

# Clinical traits for association analysis
clinical_traits <- list(
  binary = c(
    "Primary_failure_to_therapy_before_baseline",
    "Secondary_failure_to_therapy_before_baseline",
    "anti_tnf_immunogenicity_pre",
    "montreal_perianal_disease"
  ),
  categorical = c(
    "montreal_age",
    "montreal_behavior",
    "montreal_location"
  ),
  continuous = c(
    "Months_to_first_surgery",
    "Surgeries_pre_baseline"
  )
)

cat("=== Configuration loaded ===\n")
cat("Results directory:", results_dir, "\n")
cat("Key modules:", paste(key_modules, collapse = ", "), "\n")

