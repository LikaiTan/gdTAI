#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_boundary_pseudobulk_deg.R <input_dir> <output_dir>")
}
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- t(readMM(gzfile(file.path(input_dir, "pseudobulk_counts.mtx.gz"))))
genes <- readLines(file.path(input_dir, "genes.txt"), warn = FALSE)
meta <- read.csv(file.path(input_dir, "pseudobulk_metadata.csv"), check.names = FALSE)
if (nrow(counts) != length(genes) || ncol(counts) != nrow(meta)) {
  stop("Pseudobulk matrix dimensions do not match genes/metadata")
}
rownames(counts) <- genes
colnames(counts) <- meta$pseudobulk_id
if (anyDuplicated(meta$pseudobulk_id) || anyDuplicated(genes)) {
  stop("Pseudobulk IDs and genes must be unique")
}
if (!all(c("enriched", "depleted") %in% unique(meta$cluster_group))) {
  stop("Both enriched and depleted pseudobulks are required")
}

dge_all <- DGEList(counts = counts)
raw_cpm <- cpm(dge_all)
minimum_profiles <- max(10L, ceiling(0.10 * ncol(counts)))
keep <- rowSums(raw_cpm >= 1) >= minimum_profiles
if (sum(keep) < 100) {
  stop("Expression filter retained fewer than 100 genes")
}
dge <- dge_all[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)

pair_order <- unique(meta$pair_key)
enriched_col <- match(paste0(pair_order, "::enriched"), meta$pseudobulk_id)
depleted_col <- match(paste0(pair_order, "::depleted"), meta$pseudobulk_id)
if (anyNA(enriched_col) || anyNA(depleted_col)) {
  stop("Every pair must contain enriched and depleted pseudobulks")
}
delta <- log_cpm[, enriched_col, drop = FALSE] - log_cpm[, depleted_col, drop = FALSE]
colnames(delta) <- pair_order

fit <- eBayes(lmFit(delta, matrix(1, nrow = ncol(delta), ncol = 1)), trend = TRUE, robust = TRUE)
primary <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
primary$gene <- rownames(primary)
primary$median_log2FC <- apply(delta, 1, median)
primary$fraction_pairs_positive <- rowMeans(delta > 0)
primary$n_pairs <- ncol(delta)
primary <- primary[, c("gene", "logFC", "median_log2FC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "fraction_pairs_positive", "n_pairs")]
names(primary)[names(primary) == "logFC"] <- "mean_log2FC_enriched_vs_depleted"
names(primary)[names(primary) == "P.Value"] <- "paired_p_value"
names(primary)[names(primary) == "adj.P.Val"] <- "paired_fdr"

pair_source <- meta$dataset_name[enriched_col]
source_counts <- table(pair_source)
macro_sources <- names(source_counts[source_counts >= 5])
source_delta <- sapply(macro_sources, function(source) {
  rowMeans(delta[, pair_source == source, drop = FALSE])
})
if (is.null(dim(source_delta))) {
  source_delta <- matrix(source_delta, ncol = 1, dimnames = list(rownames(delta), macro_sources))
}
macro_fit <- eBayes(lmFit(source_delta, matrix(1, nrow = ncol(source_delta), ncol = 1)), trend = TRUE, robust = TRUE)
macro <- topTable(macro_fit, coef = 1, number = Inf, sort.by = "none")
macro$gene <- rownames(macro)
macro$source_median_log2FC <- apply(source_delta, 1, median)
macro$source_sign_consistency <- pmax(rowMeans(source_delta > 0), rowMeans(source_delta < 0))
macro$n_macro_datasets <- ncol(source_delta)
macro <- macro[, c("gene", "logFC", "source_median_log2FC", "P.Value", "adj.P.Val", "source_sign_consistency", "n_macro_datasets")]
names(macro)[names(macro) == "logFC"] <- "dataset_macro_mean_log2FC"
names(macro)[names(macro) == "P.Value"] <- "dataset_macro_p_value"
names(macro)[names(macro) == "adj.P.Val"] <- "dataset_macro_fdr"

result <- merge(primary, macro, by = "gene", sort = FALSE)
result$paired_significant <- result$paired_fdr < 0.05 & abs(result$mean_log2FC_enriched_vs_depleted) >= 0.25
result$robust_cross_dataset <- result$paired_significant & result$dataset_macro_fdr < 0.10 & result$source_sign_consistency >= 0.70
result <- result[order(result$paired_fdr, -abs(result$mean_log2FC_enriched_vs_depleted)), ]
write.csv(result, gzfile(file.path(output_dir, "deg_results.csv.gz")), row.names = FALSE)

source_effect <- data.frame(gene = rownames(source_delta), source_delta, check.names = FALSE)
write.csv(source_effect, gzfile(file.path(output_dir, "source_mean_log2fc.csv.gz")), row.names = FALSE)
write.csv(
  data.frame(dataset_name = names(source_counts), n_paired_samples = as.integer(source_counts)),
  file.path(output_dir, "paired_samples_by_dataset.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    metric = c("n_pseudobulks", "n_pairs", "n_genes_input", "n_genes_tested", "minimum_profiles_cpm_ge_1", "n_macro_datasets", "n_paired_significant", "n_robust_cross_dataset"),
    value = c(ncol(counts), ncol(delta), nrow(counts), nrow(result), minimum_profiles, ncol(source_delta), sum(result$paired_significant), sum(result$robust_cross_dataset))
  ),
  file.path(output_dir, "deg_run_summary.csv"), row.names = FALSE
)

message("PASS boundary pseudobulk DEG: ", ncol(delta), " pairs; ", nrow(result), " genes tested")
