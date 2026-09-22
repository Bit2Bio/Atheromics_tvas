#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# LIMMA differential analysis on log-transformed metabolomics data.
#
# - Categorical variable: design ~ 0 + variable [+ covariates]
#                         all pairwise contrasts saved as separate CSV
# - Continuous variable:  design ~ variable [+ covariates]
#                         coefficient of variable saved as CSV
#
# Covariates age/gender are auto-excluded if they are the variable of interest.
# Samples with NA in variable or covariates are dropped (pairwise complete).
#
# Usage:
#   Rscript metabolomics_limma.R \
#     --variable    metabolic_status_3g \
#     --type        categorical \
#     --covariates  age,gender \
#     --input       data/multiomics/metabolomics_matched.csv \
#     --metadata    data/multiomics/sample_metadata.csv \
#     --out_dir     results/multiomics/metabolomics/metabolic_status_3g_adj/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
  library(optparse)
  library(readr)
  library(tibble)
})

option_list <- list(
  make_option("--variable",   type = "character", help = "Clinical variable name"),
  make_option("--type",       type = "character", default = "categorical",
              help = "Variable type: categorical or continuous"),
  make_option("--covariates", type = "character", default = "age,gender",
              help = "Comma-separated covariates, or 'none'"),
  make_option("--input",      type = "character", help = "Metabolomics CSV (samples x features)"),
  make_option("--metadata",   type = "character", help = "Sample metadata CSV"),
  make_option("--out_dir",    type = "character", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------
met  <- read_csv(opt$input,    show_col_types = FALSE)
meta <- read_csv(opt$metadata, show_col_types = FALSE)

# metabolomics: samples x features → transpose to features x samples for limma
sample_col <- colnames(met)[1]
met <- as.data.frame(met)
rownames(met) <- met[[sample_col]]
met <- met[, sapply(met, is.numeric), drop = FALSE]   # keep only numeric columns
mat <- t(as.matrix(met))   # features x samples
storage.mode(mat) <- "double"

# metadata: align to matrix columns
meta <- as.data.frame(meta)
rownames(meta) <- meta[[colnames(meta)[1]]]
meta <- meta[colnames(mat), , drop = FALSE]

cat(sprintf("Matrix: %d features x %d samples\n", nrow(mat), ncol(mat)))

# --------------------------------------------------------------------------
# Parse covariates — exclude if same as variable of interest
# --------------------------------------------------------------------------
if (tolower(opt$covariates) == "none") {
  covs <- character(0)
} else {
  covs <- trimws(strsplit(opt$covariates, ",")[[1]])
  covs <- covs[covs != opt$variable]
}
cat(sprintf("Variable: %s  |  Type: %s  |  Covariates: %s\n",
            opt$variable, opt$type,
            if (length(covs) == 0) "none" else paste(covs, collapse = ", ")))

# --------------------------------------------------------------------------
# Filter samples: drop NA in variable + covariates
# --------------------------------------------------------------------------
keep_vars <- c(opt$variable, covs)
complete   <- complete.cases(meta[, keep_vars, drop = FALSE])
meta_f     <- meta[complete, , drop = FALSE]
mat_f      <- mat[, complete, drop = FALSE]
cat(sprintf("Samples after NA removal: %d\n", ncol(mat_f)))

# --------------------------------------------------------------------------
# Output directory
# --------------------------------------------------------------------------
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Build design and run LIMMA
# --------------------------------------------------------------------------
if (opt$type == "categorical") {

  grp <- factor(meta_f[[opt$variable]])
  levels(grp) <- paste0("G", levels(grp))   # ensure valid R names
  cat(sprintf("Groups: %s\n", paste(levels(grp), collapse = ", ")))

  # Design matrix: one column per group
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)

  # Add covariates
  for (cov in covs) {
    cov_vals <- meta_f[[cov]]
    if (is.character(cov_vals) || is.factor(cov_vals)) {
      dummy <- model.matrix(~ factor(cov_vals))[, -1, drop = FALSE]
      colnames(dummy) <- gsub("factor\\(cov_vals\\)", cov, colnames(dummy))
      design <- cbind(design, dummy)
    } else {
      design <- cbind(design, cov_vals)
      colnames(design)[ncol(design)] <- cov
    }
  }

  fit <- lmFit(mat_f, design)

  # All pairwise contrasts
  lvls     <- levels(grp)
  pairs    <- combn(lvls, 2, simplify = FALSE)
  contrast_strings <- sapply(pairs, function(p) paste0(p[2], " - ", p[1]))
  contrast_mat <- makeContrasts(contrasts = contrast_strings, levels = design)

  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)

  for (i in seq_along(contrast_strings)) {
    cs   <- contrast_strings[i]
    pair <- pairs[[i]]
    res  <- topTable(fit2, coef = i, number = Inf, sort.by = "P")
    res  <- rownames_to_column(res, var = "feature")
    fname <- sprintf("contrast_%s_vs_%s.csv", pair[2], pair[1])
    write_csv(res, file.path(opt$out_dir, fname))
    cat(sprintf("  %s: %d FDR<0.05\n", cs, sum(res$adj.P.Val < 0.05, na.rm = TRUE)))
  }

} else if (opt$type == "continuous") {

  cont_var <- meta_f[[opt$variable]]

  # Design matrix: intercept + variable + covariates
  design <- model.matrix(~ cont_var)
  colnames(design) <- c("intercept", opt$variable)

  for (cov in covs) {
    cov_vals <- meta_f[[cov]]
    if (is.character(cov_vals) || is.factor(cov_vals)) {
      dummy <- model.matrix(~ factor(cov_vals))[, -1, drop = FALSE]
      colnames(dummy) <- gsub("factor\\(cov_vals\\)", cov, colnames(dummy))
      design <- cbind(design, dummy)
    } else {
      design <- cbind(design, cov_vals)
      colnames(design)[ncol(design)] <- cov
    }
  }

  fit  <- lmFit(mat_f, design)
  fit  <- eBayes(fit)

  res  <- topTable(fit, coef = opt$variable, number = Inf, sort.by = "P")
  res  <- rownames_to_column(res, var = "feature")
  write_csv(res, file.path(opt$out_dir, paste0(opt$variable, ".csv")))
  cat(sprintf("  %s: %d FDR<0.05\n", opt$variable, sum(res$adj.P.Val < 0.05, na.rm = TRUE)))

} else {
  stop("--type must be 'categorical' or 'continuous'")
}

cat(sprintf("\nDone. Results saved to: %s\n", opt$out_dir))
