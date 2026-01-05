########################################
# WGCNA Analysis for Crohn's Disease
# 03: Module Analysis (Eigengenes, Associations, GO, kME)
########################################

source("00_config.R")
source("01_data_preprocessing.R")
source("02_wgcna_network.R")

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(pheatmap)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)

cat("\n=== Step 1: Module eigengenes ===\n")

# Create gene-to-module mapping
gene_ids <- colnames(expr_cd)
gene_module_df <- data.frame(
  ensembl_raw    = gene_ids,
  module_numeric = moduleLabels,
  module_color   = moduleColors,
  stringsAsFactors = FALSE
)

# Remove grey module (module 0)
gene_module_df <- gene_module_df[gene_module_df$module_numeric != 0, ]

# Map numeric labels to module names
module_numeric_unique <- sort(unique(gene_module_df$module_numeric))
module_name_map <- setNames(
  paste0("Module_", seq_along(module_numeric_unique)),
  module_numeric_unique
)
gene_module_df$Module <- module_name_map[as.character(gene_module_df$module_numeric)]

# Create color-to-module mapping
color_to_module <- gene_module_df %>%
  dplyr::distinct(module_color, Module)

# Calculate module eigengenes
MEs0 <- moduleEigengenes(expr_cd, colors = moduleColors)$eigengenes
MEs0 <- orderMEs(MEs0)  # Standard WGCNA ordering

# Remove grey module
me_colors <- gsub("^ME", "", colnames(MEs0))
keep_me <- me_colors != "grey"
MEs <- MEs0[, keep_me, drop = FALSE]
me_colors <- me_colors[keep_me]

# Map colors to module names
module_names <- color_to_module$Module[match(me_colors, color_to_module$module_color)]
if (any(is.na(module_names))) {
  stop("Some eigengene colors could not be mapped to module names")
}
colnames(MEs) <- module_names

# Re-order modules numerically
module_levels <- sort(unique(module_names))
MEs <- MEs[, module_levels, drop = FALSE]

# Save eigengenes
eigengene_out <- file.path(results_dir, "module_eigengenes.csv")
write.csv(MEs, eigengene_out, row.names = TRUE)
cat("Module eigengenes saved to:", eigengene_out, "\n")

# Module eigengene correlation
me_cor <- cor(MEs, use = "pairwise.complete.obs")
pdf(file.path(results_dir, "module_eigengene_correlation.pdf"),
    width = 7, height = 12)
pheatmap(
  me_cor,
  main = "Module eigengene correlations",
  clustering_method = "average",
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.2f",
  fontsize_number = 8
)
dev.off()

cat("\n=== Step 2: GO enrichment analysis ===\n")

gene_module_df$ensembl <- sub("\\..*$", "", gene_module_df$ensembl_raw)
modules_to_test <- sort(unique(gene_module_df$Module))

all_go_table <- list()

for (m in modules_to_test) {
  cat("Processing", m, "...\n")
  
  genes_m <- gene_module_df$ensembl[gene_module_df$Module == m]
  genes_m <- unique(genes_m[!is.na(genes_m) & genes_m != ""])
  
  if (length(genes_m) < 10) {
    cat("  Skipping", m, "(<10 genes)\n")
    next
  }
  
  # GO enrichment
  ego <- tryCatch(
    enrichGO(
      gene          = genes_m,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENSEMBL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    ),
    error = function(e) {
      cat("  enrichGO failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  No significant GO terms\n")
    next
  }
  
  # Simplify redundant terms
  ego_simpl <- tryCatch(
    clusterProfiler::simplify(ego),
    error = function(e) {
      cat("  simplify() failed, using original\n")
      ego
    }
  )
  
  ego_df <- as.data.frame(ego_simpl)
  if (nrow(ego_df) == 0) {
    cat("  No terms after simplification\n")
    next
  }
  
  ego_df$Module <- m
  all_go_table[[m]] <- ego_df
  
  # Create individual dotplots
  n_show <- min(10, nrow(ego_df))
  n_show <- max(n_show, 2)
  
  pdf(file.path(go_dir, paste0(m, "_GO_enrichment.pdf")),
      width = 7, height = 10)
  print(
    dotplot(ego_simpl, showCategory = n_show,
            title = paste0(m, " GO Biological Process"))
  )
  dev.off()
}

# Save combined GO results
if (length(all_go_table) > 0) {
  go_combined <- dplyr::bind_rows(all_go_table) %>%
    dplyr::select(Module, ID, Description, GeneRatio, BgRatio,
                  pvalue, p.adjust, qvalue, Count, everything())
  
  write.csv(go_combined,
            file.path(go_dir, "GO_enrichment_all_modules.csv"),
            row.names = FALSE)
  cat("GO enrichment saved to:", file.path(go_dir, "GO_enrichment_all_modules.csv"), "\n")
}

cat("\n=== Step 3: kME analysis (hub genes) ===\n")

kME_table <- signedKME(expr_cd, MEs)
colnames(kME_table) <- colnames(MEs)

kME_df <- as.data.frame(kME_table) %>%
  tibble::rownames_to_column("ensembl_raw") %>%
  dplyr::left_join(
    gene_module_df %>% dplyr::select(ensembl_raw, ensembl, Module),
    by = "ensembl_raw"
  )

# Calculate kME_self (correlation with own module)
kME_df$kME_self <- NA_real_
for (m in module_levels) {
  idx <- which(kME_df$Module == m)
  if (length(idx) == 0) next
  kME_df$kME_self[idx] <- kME_df[idx, m]
}

# Add gene symbols
sym_map_raw <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(kME_df$ensembl),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "list"
)

sym_map <- vapply(
  sym_map_raw,
  function(x) if (length(x) == 0) NA_character_ else as.character(x[1]),
  FUN.VALUE = character(1)
)

kME_df$symbol <- unname(sym_map[kME_df$ensembl])

kME_df <- kME_df %>%
  dplyr::arrange(Module, dplyr::desc(kME_self))

# Save full kME table
write.csv(kME_df, 
          file.path(results_dir, "kME_all_genes.csv"),
          row.names = FALSE)

# Hub summary per module
hub_summary <- kME_df %>%
  dplyr::filter(!is.na(Module)) %>%
  dplyr::group_by(Module) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    mean_kME = mean(kME_self, na.rm = TRUE),
    median_kME = median(kME_self, na.rm = TRUE),
    top5_symbols = paste(
      symbol[order(-kME_self)][!is.na(symbol)][1:min(5, sum(!is.na(symbol)))],
      collapse = "; "
    ),
    top5_kME = paste(
      round(kME_self[order(-kME_self)][!is.na(kME_self)][1:min(5, sum(!is.na(kME_self)))], 3),
      collapse = "; "
    ),
    .groups = "drop"
  )

write.csv(hub_summary,
          file.path(results_dir, "kME_hub_summary.csv"),
          row.names = FALSE)

cat("\n=== Step 4: Module-trait associations ===\n")

# Prepare metadata
if ("montreal_perianal_disease" %in% names(meta_cd)) {
  meta_cd$montreal_perianal_disease <- factor(
    ifelse(meta_cd$montreal_perianal_disease == "p", "Yes", "No")
  )
}

# Convert variables to appropriate types
binary_vars <- intersect(clinical_traits$binary, names(meta_cd))
meta_cd[binary_vars] <- lapply(meta_cd[binary_vars], factor)

categorical_vars <- intersect(clinical_traits$categorical, names(meta_cd))
meta_cd[categorical_vars] <- lapply(meta_cd[categorical_vars], factor)

numeric_vars <- intersect(clinical_traits$continuous, names(meta_cd))
meta_cd[numeric_vars] <- lapply(
  meta_cd[numeric_vars],
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

# Create trait type mapping
trait_type <- c(
  setNames(rep("binary", length(clinical_traits$binary)), clinical_traits$binary),
  setNames(rep("categorical", length(clinical_traits$categorical)), clinical_traits$categorical),
  setNames(rep("continuous", length(clinical_traits$continuous)), clinical_traits$continuous)
)

clinical_vars <- intersect(names(trait_type), names(meta_cd))

# Association test function
assoc_test <- function(x, y, type) {
  ok <- is.finite(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(y) < 5L) {
    return(data.frame(p_raw = NA_real_, effect = NA_real_, n = length(y)))
  }
  
  if (type == "binary") {
    y <- factor(y)
    if (nlevels(y) != 2L || any(table(y) == 0L)) {
      return(data.frame(p_raw = NA_real_, effect = NA_real_, n = length(y)))
    }
    df <- data.frame(x = x, y = y)
    p_val <- suppressWarnings(wilcox.test(x ~ y, exact = FALSE)$p.value)
    eff <- tryCatch(
      rstatix::wilcox_effsize(df, x ~ y)$effsize[1],
      error = function(e) NA_real_
    )
    return(data.frame(p_raw = p_val, effect = eff, n = length(y)))
  }
  
  if (type == "categorical") {
    y <- factor(y)
    if (nlevels(y) < 2L || any(table(y) == 0L)) {
      return(data.frame(p_raw = NA_real_, effect = NA_real_, n = length(y)))
    }
    p_val <- suppressWarnings(kruskal.test(x ~ y)$p.value)
    eff <- tryCatch(
      rstatix::kruskal_effsize(data.frame(x = x, y = y), x ~ y)$effsize,
      error = function(e) NA_real_
    )
    return(data.frame(p_raw = p_val, effect = eff, n = length(y)))
  }
  
  if (type == "continuous") {
    y <- suppressWarnings(as.numeric(as.character(y)))
    ok2 <- is.finite(y)
    x <- x[ok2]
    y <- y[ok2]
    if (length(y) < 5L || length(unique(y)) < 2L) {
      return(data.frame(p_raw = NA_real_, effect = NA_real_, n = length(y)))
    }
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    return(data.frame(
      p_raw = unname(ct$p.value),
      effect = unname(ct$estimate),
      n = length(y)
    ))
  }
  
  data.frame(p_raw = NA_real_, effect = NA_real_, n = length(y))
}

# Run module-trait associations
res_list <- list()

for (var in clinical_vars) {
  ttype <- trait_type[[var]]
  y <- meta_cd[[var]]
  
  for (m in colnames(MEs)) {
    x <- MEs[, m]
    at <- assoc_test(x, y, ttype)
    at$clinical_var <- var
    at$module <- m
    at$test_type <- ttype
    res_list[[length(res_list) + 1L]] <- at
  }
}

res_tbl <- dplyr::bind_rows(res_list)

# FDR correction within each clinical variable
res_tbl <- res_tbl %>%
  dplyr::group_by(clinical_var) %>%
  dplyr::mutate(p_fdr = p.adjust(p_raw, method = "BH")) %>%
  dplyr::ungroup()

# Create FDR matrix for heatmap
assoc_fdr_mat <- res_tbl %>%
  dplyr::select(clinical_var, module, p_fdr) %>%
  tidyr::pivot_wider(
    names_from = module,
    values_from = p_fdr
  ) %>%
  tibble::column_to_rownames("clinical_var") %>%
  as.matrix()

# Prepare association results table
res_assoc_df <- res_tbl %>%
  dplyr::mutate(
    module_clean = gsub("Module_", "Module ", module)
  ) %>%
  dplyr::relocate(
    clinical_var, test_type, module, module_clean,
    p_raw, p_fdr, effect, n
  )

# Save associations
write.csv(res_assoc_df,
          file.path(results_dir, "module_clinical_associations.csv"),
          row.names = FALSE)

cat("Module-trait associations saved\n")

cat("\n=== Module analysis complete ===\n")

