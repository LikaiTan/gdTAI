#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(R.utils)
})


parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for %s", key))
    }
    out[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}


require_arg <- function(args, key) {
  value <- args[[key]]
  if (is.null(value) || identical(value, "")) {
    stop(sprintf("Missing required argument --%s", key))
  }
  value
}


maybe_gunzip_rds <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(path)
  }
  ext <- sub("^.*(\\.[Rr][Dd][Ss])\\.gz$", "\\1", basename(path))
  if (identical(ext, basename(path))) {
    ext <- ".rds"
  }
  tmp <- tempfile(fileext = ext)
  R.utils::gunzip(path, destname = tmp, overwrite = TRUE, remove = FALSE)
  tmp
}


write_sparse_payload <- function(counts, out_prefix, metadata = NULL) {
  if (!inherits(counts, "sparseMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  matrix_path <- paste0(out_prefix, "_counts.mtx")
  genes_path <- paste0(out_prefix, "_genes.tsv.gz")
  barcodes_path <- paste0(out_prefix, "_barcodes.tsv.gz")
  metadata_path <- paste0(out_prefix, "_metadata.csv.gz")

  Matrix::writeMM(counts, matrix_path)

  genes <- data.frame(
    gene_id = rownames(counts),
    gene_name = rownames(counts),
    stringsAsFactors = FALSE
  )
  gz_genes <- gzfile(genes_path, open = "wt")
  write.table(genes, gz_genes, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(gz_genes)

  barcodes <- data.frame(cell_id = colnames(counts), stringsAsFactors = FALSE)
  gz_barcodes <- gzfile(barcodes_path, open = "wt")
  write.table(barcodes, gz_barcodes, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(gz_barcodes)

  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    metadata$cell_id <- rownames(metadata)
    metadata <- metadata[, c("cell_id", setdiff(colnames(metadata), "cell_id"))]
    gz_meta <- gzfile(metadata_path, open = "wt")
    write.csv(metadata, gz_meta, row.names = FALSE, quote = TRUE)
    close(gz_meta)
  }

  cat(sprintf("Wrote payload prefix: %s\n", out_prefix))
}


export_dgcmatrix <- function(input_path, out_prefix) {
  temp_rds <- maybe_gunzip_rds(input_path)
  on.exit({
    if (!identical(temp_rds, input_path) && file.exists(temp_rds)) {
      unlink(temp_rds)
    }
  }, add = TRUE)

  counts <- readRDS(temp_rds)
  if (!inherits(counts, "Matrix")) {
    stop("Expected an R Matrix object for mode=dgcmatrix")
  }
  write_sparse_payload(counts, out_prefix, metadata = NULL)
}


export_seurat <- function(input_path, out_prefix, assay_name = NULL) {
  temp_rds <- maybe_gunzip_rds(input_path)
  on.exit({
    if (!identical(temp_rds, input_path) && file.exists(temp_rds)) {
      unlink(temp_rds)
    }
  }, add = TRUE)

  obj <- readRDS(temp_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Expected a Seurat object for mode=seurat")
  }

  assay_to_use <- assay_name
  if (is.null(assay_to_use) || identical(assay_to_use, "")) {
    assay_to_use <- DefaultAssay(obj)
  }

  counts <- GetAssayData(obj, assay = assay_to_use, layer = "counts")
  metadata <- obj@meta.data
  write_sparse_payload(counts, out_prefix, metadata = metadata)
}


main <- function() {
  args <- parse_args()
  mode <- require_arg(args, "mode")
  input_path <- require_arg(args, "input")
  out_prefix <- require_arg(args, "out-prefix")
  assay_name <- args[["assay"]]

  if (identical(mode, "dgcmatrix")) {
    export_dgcmatrix(input_path, out_prefix)
  } else if (identical(mode, "seurat")) {
    export_seurat(input_path, out_prefix, assay_name = assay_name)
  } else {
    stop(sprintf("Unsupported mode: %s", mode))
  }
}


if (sys.nframe() == 0) {
  main()
}
