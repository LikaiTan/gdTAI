#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})


parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop(sprintf("Invalid argument sequence near: %s", key))
    }
    out[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}


required_arg <- function(args, key) {
  value <- args[[key]]
  if (is.null(value) || identical(value, "")) {
    stop(sprintf("Missing required argument --%s", key))
  }
  value
}


clean_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  for (name in colnames(metadata)) {
    value <- metadata[[name]]
    if (is.factor(value)) {
      metadata[[name]] <- as.character(value)
    } else if (inherits(value, c("Date", "POSIXct", "POSIXlt"))) {
      metadata[[name]] <- as.character(value)
    } else if (is.list(value)) {
      metadata[[name]] <- vapply(
        value,
        function(item) paste(as.character(item), collapse = ";"),
        character(1)
      )
    }
  }
  metadata
}


write_csv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(data, connection, row.names = FALSE, quote = TRUE, na = "")
}


main <- function() {
  args <- parse_args()
  input_rds <- required_arg(args, "input-rds")
  output_prefix <- required_arg(args, "output-prefix")
  compartment <- required_arg(args, "compartment")
  lineage <- toupper(required_arg(args, "lineage"))
  feature_manifest <- required_arg(args, "feature-manifest")
  overwrite <- identical(tolower(args[["overwrite"]]), "true")

  if (!compartment %in% c("lung", "lymph_node")) {
    stop("--compartment must be lung or lymph_node")
  }
  if (!lineage %in% c("CD4", "CD8")) {
    stop("--lineage must be CD4 or CD8")
  }

  outputs <- c(
    paste0(output_prefix, "_model_gene_counts.mtx.gz"),
    paste0(output_prefix, "_model_genes.tsv"),
    paste0(output_prefix, "_cell_metadata.csv.gz"),
    paste0(output_prefix, "_export_manifest.csv")
  )
  existing <- outputs[file.exists(outputs)]
  if (length(existing) > 0 && !overwrite) {
    stop(paste0(
      "Output exists; pass --overwrite true to replace: ",
      paste(existing, collapse = ", ")
    ))
  }
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  object <- readRDS(input_rds)
  if (!inherits(object, "Seurat")) {
    stop("Expected a Seurat object")
  }
  if (!"RNA" %in% Assays(object)) {
    stop("The Seurat object does not contain an RNA assay")
  }

  feature_table <- read.csv(
    feature_manifest,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_feature_columns <- c("feature_type", "gene")
  missing_feature_columns <- setdiff(required_feature_columns, colnames(feature_table))
  if (length(missing_feature_columns) > 0) {
    stop(paste0(
      "Feature manifest is missing columns: ",
      paste(missing_feature_columns, collapse = ", ")
    ))
  }
  model_genes <- feature_table$gene[
    feature_table$feature_type == "gene_log1p_cp10k" &
      !is.na(feature_table$gene) & feature_table$gene != ""
  ]
  if (length(model_genes) == 0 || anyDuplicated(model_genes)) {
    stop("Model feature manifest has no genes or contains duplicated genes")
  }

  counts <- GetAssayData(object, assay = "RNA", layer = "counts")
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  if (!identical(colnames(counts), rownames(object@meta.data))) {
    stop("RNA count columns and metadata rows are not identically ordered")
  }
  missing_genes <- setdiff(model_genes, rownames(counts))
  if (length(missing_genes) > 0) {
    stop(paste0(
      "RNA counts are missing ", length(missing_genes),
      " model genes: ", paste(head(missing_genes, 30), collapse = ", ")
    ))
  }

  selected_counts <- as(counts[model_genes, , drop = FALSE], "dgCMatrix")
  total_counts <- as.numeric(Matrix::colSums(counts))
  detected_genes <- as.integer(diff(counts@p))
  if (any(!is.finite(total_counts)) || any(total_counts <= 0)) {
    stop("Every exported cell must have a positive finite RNA library size")
  }

  source_cells <- colnames(counts)
  metadata <- clean_metadata(object@meta.data)
  metadata$cell_id_original <- source_cells
  metadata$cell_id <- paste(
    "GSE305372", compartment, lineage, source_cells, sep = ":"
  )
  metadata$gse_id <- "GSE305372"
  metadata$source_compartment <- compartment
  metadata$author_lineage <- lineage
  metadata$author_object <- paste(compartment, lineage, sep = "_")
  metadata$row_sum_counts_layer <- total_counts
  metadata$n_detected_genes_counts_layer <- detected_genes

  if ("umap" %in% Reductions(object)) {
    umap <- Embeddings(object, reduction = "umap")[source_cells, , drop = FALSE]
    metadata$UMAP_1 <- as.numeric(umap[, 1])
    metadata$UMAP_2 <- as.numeric(umap[, 2])
  }

  ncount_max_abs_diff <- NA_real_
  if ("nCount_RNA" %in% colnames(metadata)) {
    ncount_max_abs_diff <- max(
      abs(as.numeric(metadata$nCount_RNA) - total_counts),
      na.rm = TRUE
    )
  }

  temporary_mtx <- paste0(output_prefix, "_model_gene_counts.mtx")
  Matrix::writeMM(selected_counts, temporary_mtx)
  R.utils::gzip(
    temporary_mtx,
    destname = paste0(temporary_mtx, ".gz"),
    overwrite = TRUE,
    remove = TRUE
  )
  writeLines(
    model_genes,
    paste0(output_prefix, "_model_genes.tsv"),
    useBytes = TRUE
  )
  write_csv_gz(metadata, paste0(output_prefix, "_cell_metadata.csv.gz"))

  manifest <- data.frame(
    gse_id = "GSE305372",
    compartment = compartment,
    lineage = lineage,
    input_rds = normalizePath(input_rds),
    input_size_bytes = file.info(input_rds)$size,
    seurat_object_version = as.character(object@version),
    assay = "RNA",
    source_cells = ncol(object),
    exported_cells = ncol(object),
    genes_in_rna_assay = nrow(counts),
    model_gene_count = length(model_genes),
    inclusion_rule = "all_cells_in_author_finalized_T_cell_object",
    ncount_max_abs_diff = ncount_max_abs_diff,
    stringsAsFactors = FALSE
  )
  write.csv(
    manifest,
    paste0(output_prefix, "_export_manifest.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  cat(sprintf(
    "Exported %s %s: %s / %s cells and %s model genes\n",
    compartment,
    lineage,
    format(ncol(object), big.mark = ","),
    format(ncol(object), big.mark = ","),
    format(length(model_genes), big.mark = ",")
  ))
}


main()
