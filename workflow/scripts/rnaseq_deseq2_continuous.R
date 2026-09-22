#!/usr/bin/env Rscript
# DESeq2 differential expression for a single continuous clinical variable.
# Wald test with optional covariates. Output: one CSV (all genes, sorted by PValue).
#
# Usage:
#   Rscript rnaseq_deseq2_continuous.R \
#     --variable   age \
#     --covariates gender \
#     --matrix     data/processed/transcriptomics_tvas.csv \
#     --metadata   data/processed/clinical_tvas.csv \
#     --annotation data/raw/transcriptomics/gene_annotation.csv \
#     --out_dir    results/transcriptomics/deseq2/age/

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(readr)
  library(dplyr)
})

option_list <- list(
  make_option("--variable",   type = "character", help = "Continuous clinical variable name"),
  make_option("--covariates", type = "character", default = "age,gender",
              help = "Comma-separated covariates, or 'none' (default: age,gender)"),
  make_option("--matrix",     type = "character", help = "Raw count matrix CSV (genes x samples)"),
  make_option("--metadata",   type = "character", help = "Sample metadata CSV"),
  make_option("--out_dir",    type = "character", help = "Output directory"),
  make_option("--annotation", type = "character", default = "data/raw/transcriptomics/gene_annotation.csv",
              help = "Gene annotation CSV with gene_id and gene_name columns")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Load data ─────────────────────────────────────────────────────────────────
counts <- read.table(opt$matrix, sep = ",", header = TRUE,
                     row.names = 1, check.names = FALSE)
meta   <- read_csv(opt$metadata, show_col_types = FALSE) %>%
          tibble::column_to_rownames(colnames(.)[1])
annot  <- read_csv(opt$annotation, show_col_types = FALSE) %>%
          select(gene_id, gene_name, gene_type) %>%
          distinct(gene_id, .keep_all = TRUE)

# ── Align samples ─────────────────────────────────────────────────────────────
common <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common, drop = FALSE]
meta   <- meta[common, , drop = FALSE]
cat("Matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ── Parse covariates ──────────────────────────────────────────────────────────
if (tolower(opt$covariates) == "none") {
  covs <- character(0)
} else {
  covs <- trimws(strsplit(opt$covariates, ",")[[1]])
  covs <- covs[covs != opt$variable]
}
cat("Variable:", opt$variable, " | Covariates:", paste(covs, collapse = ", "), "\n")

# ── Filter samples: drop NA in variable + covariates ─────────────────────────
keep_vars <- c(opt$variable, covs)
keep_vars <- keep_vars[keep_vars %in% colnames(meta)]
complete  <- complete.cases(meta[, keep_vars, drop = FALSE])
meta   <- meta[complete, , drop = FALSE]
counts <- counts[, rownames(meta), drop = FALSE]
counts <- round(counts)   # DESeq2 requires integers
cat("Samples after NA removal:", ncol(counts), "\n")

# ── Continuous variable ───────────────────────────────────────────────────────
meta[[opt$variable]] <- as.numeric(meta[[opt$variable]])
cat("Variable range:", range(meta[[opt$variable]], na.rm = TRUE), "\n")

# ── Build DESeq2 design ───────────────────────────────────────────────────────
if (length(covs) == 0) {
  design_formula <- as.formula(paste("~", opt$variable))
} else {
  for (cov in covs) meta[[cov]] <- scale(as.numeric(meta[[cov]]))[, 1]
  design_formula <- as.formula(paste("~", opt$variable, "+", paste(covs, collapse = " + ")))
}

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = design_formula
)

dds <- DESeq(dds, quiet = TRUE)
cat("Genes:", nrow(dds), "\n")

# ── Results ───────────────────────────────────────────────────────────────────
res_raw <- results(dds, name = opt$variable, independentFiltering = TRUE)
res <- as.data.frame(res_raw) %>%
  tibble::rownames_to_column("gene_id") %>%
  rename(logFC = log2FoldChange, PValue = pvalue, FDR = padj) %>%
  left_join(annot, by = "gene_id") %>%
  arrange(PValue)

n_sig <- sum(!is.na(res$FDR) & res$FDR < 0.05)
cat("FDR<0.05:", n_sig, "genes\n")

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(res, file.path(opt$out_dir, paste0(opt$variable, "_results.csv")))

cat("\nDone. Results saved to:", opt$out_dir, "\n")
