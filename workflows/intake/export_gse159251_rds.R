#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

parse_args <- function(args) {
  result <- list(input = NULL, output_dir = NULL, tcr_column = "")
  index <- 1L
  while (index <= length(args)) {
    flag <- args[[index]]
    if (!flag %in% c("--input", "--output-dir", "--tcr-column") || index == length(args)) {
      stop("Usage: export_gse159251_rds.R --input FILE --output-dir DIR [--tcr-column NAME]")
    }
    value <- args[[index + 1L]]
    if (flag == "--input") result$input <- value
    if (flag == "--output-dir") result$output_dir <- value
    if (flag == "--tcr-column") result$tcr_column <- value
    index <- index + 2L
  }
  if (is.null(result$input) || is.null(result$output_dir)) {
    stop("--input and --output-dir are required")
  }
  result
}

is_gzip <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  magic <- readBin(connection, what = "raw", n = 2L)
  length(magic) == 2L && identical(as.integer(magic), c(31L, 139L))
}

copy_connection <- function(input_connection, output_connection) {
  repeat {
    chunk <- readBin(input_connection, what = "raw", n = 8L * 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, output_connection)
  }
}

peel_gzip_layers <- function(path) {
  current <- normalizePath(path, mustWork = TRUE)
  temporary_files <- character()
  layers <- 0L
  while (is_gzip(current)) {
    layers <- layers + 1L
    if (layers > 4L) stop("Refusing input with more than four nested gzip layers")
    destination <- tempfile(pattern = "gse159251_uncompressed_", fileext = ".rds")
    input_connection <- gzfile(current, open = "rb")
    output_connection <- file(destination, open = "wb")
    tryCatch(
      copy_connection(input_connection, output_connection),
      finally = {
        close(input_connection)
        close(output_connection)
      }
    )
    temporary_files <- c(temporary_files, destination)
    current <- destination
  }
  if (layers < 2L) {
    stop(sprintf("Expected a double-gzipped GSE159251 RDS, detected %d gzip layer(s)", layers))
  }
  list(path = current, layers = layers, temporary_files = temporary_files)
}

extract_seurat_payload <- function(object) {
  if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object")
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject is required to extract the RNA counts layer")
  }
  assay_name <- if ("RNA" %in% names(object@assays)) "RNA" else SeuratObject::DefaultAssay(object)
  counts <- tryCatch(
    SeuratObject::LayerData(object[[assay_name]], layer = "counts"),
    error = function(error) NULL
  )
  if (is.null(counts)) {
    counts <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay_name, slot = "counts"),
      error = function(error) NULL
    )
  }
  if (is.null(counts)) stop("Unable to extract the Seurat RNA counts layer/slot")
  metadata <- object[[]]
  list(counts = counts, metadata = metadata, assay = assay_name)
}

normalize_tcr_sentinel <- function(values) {
  values <- as.character(values)
  values[is.na(values)] <- ""
  sentinel <- tolower(trimws(values)) %in% c("notcr", "no_tcr", "no tcr", "unassigned")
  values[sentinel] <- ""
  values
}

find_tcr_column <- function(metadata, requested) {
  if (nzchar(requested)) {
    if (!requested %in% colnames(metadata)) stop(sprintf("Requested TCR column not found: %s", requested))
    return(requested)
  }
  candidates <- character()
  for (column in colnames(metadata)) {
    values <- normalize_tcr_sentinel(metadata[[column]])
    nonblank <- values[nzchar(trimws(values))]
    if (length(nonblank) && all(lengths(regmatches(nonblank, gregexpr("\\|", nonblank))) == 1L)) {
      candidates <- c(candidates, column)
    }
  }
  preferred <- candidates[grepl("tcr|cdr3|tra|trb", candidates, ignore.case = TRUE)]
  if (length(preferred) == 1L) return(preferred)
  if (length(candidates) == 1L) return(candidates)
  stop(sprintf("Could not identify one unambiguous TRA_CDR3|TRB_CDR3 column; candidates=%s", paste(candidates, collapse = ",")))
}

write_gzip_text <- function(path, writer) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection))
  writer(connection)
}

gzip_file <- function(source, destination) {
  input_connection <- file(source, open = "rb")
  output_connection <- gzfile(destination, open = "wb")
  tryCatch(
    copy_connection(input_connection, output_connection),
    finally = {
      close(input_connection)
      close(output_connection)
    }
  )
}

json_escape <- function(value) {
  value <- gsub("\\\\", "\\\\\\\\", value)
  value <- gsub('"', '\\"', value, fixed = TRUE)
  gsub("\n", "\\n", value, fixed = TRUE)
}

main <- function() {
args <- parse_args(commandArgs(trailingOnly = TRUE))
input_path <- normalizePath(args$input, mustWork = TRUE)
output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop(sprintf("Output directory must be empty: %s", output_dir))
}

peeled <- peel_gzip_layers(input_path)
on.exit(unlink(peeled$temporary_files), add = TRUE)
object <- readRDS(peeled$path)
payload <- extract_seurat_payload(object)
counts <- payload$counts
if (!inherits(counts, "sparseMatrix")) stop("Seurat RNA counts are dense; refusing conversion")
counts <- as(counts, "dgCMatrix")
if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != floor(counts@x))) {
  stop("Seurat RNA counts are not finite nonnegative integers")
}
if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) stop("Counts have blank gene names")
if (is.null(colnames(counts)) || any(!nzchar(colnames(counts))) || anyDuplicated(colnames(counts))) {
  stop("Counts have blank or duplicate cell barcodes")
}

metadata <- payload$metadata
if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
  stop("Seurat metadata has blank or duplicate row names")
}
if (!setequal(colnames(counts), rownames(metadata))) {
  stop("Seurat counts and metadata barcodes do not match exactly")
}
metadata <- metadata[colnames(counts), , drop = FALSE]
tcr_column <- find_tcr_column(metadata, args$tcr_column)
tcr_values <- normalize_tcr_sentinel(metadata[[tcr_column]])
nonblank <- tcr_values[nzchar(trimws(tcr_values))]
if (length(nonblank) && any(lengths(regmatches(nonblank, gregexpr("\\|", nonblank))) != 1L)) {
  stop("Embedded TCR values are not consistently TRA_CDR3|TRB_CDR3")
}

matrix_plain <- tempfile(pattern = "gse159251_matrix_", fileext = ".mtx")
on.exit(unlink(matrix_plain), add = TRUE)
Matrix::writeMM(counts, matrix_plain)
gzip_file(matrix_plain, file.path(output_dir, "matrix.mtx.gz"))

feature_ids <- rownames(counts)
write_gzip_text(file.path(output_dir, "features.tsv.gz"), function(connection) {
  write.table(
    data.frame(feature_id = feature_ids, gene_symbol = feature_ids, stringsAsFactors = FALSE),
    connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
})
write_gzip_text(file.path(output_dir, "barcodes.tsv.gz"), function(connection) {
  writeLines(colnames(counts), connection)
})

metadata_export <- metadata
metadata_export$barcode <- rownames(metadata_export)
metadata_export$embedded_tcr <- tcr_values
metadata_export$source_tcr_column <- tcr_column
metadata_export <- metadata_export[, c("barcode", "embedded_tcr", "source_tcr_column", setdiff(colnames(metadata_export), c("barcode", "embedded_tcr", "source_tcr_column"))), drop = FALSE]
write_gzip_text(file.path(output_dir, "metadata.csv.gz"), function(connection) {
  write.csv(metadata_export, connection, row.names = FALSE, na = "")
})

provenance <- c(
  "{",
  sprintf('  "source_path": "%s",', json_escape(input_path)),
  sprintf('  "gzip_layers": %d,', peeled$layers),
  sprintf('  "object_class": "%s",', json_escape(paste(class(object), collapse = ";"))),
  sprintf('  "assay": "%s",', json_escape(payload$assay)),
  sprintf('  "tcr_column": "%s",', json_escape(tcr_column)),
  sprintf('  "genes": %d,', nrow(counts)),
  sprintf('  "cells": %d,', ncol(counts)),
  sprintf('  "nnz": %d,', length(counts@x)),
  '  "count_matrix_state": "raw_sparse_integer_counts",',
  '  "tcr_schema": "partial_embedded_tra_trb_cdr3"',
  "}"
)
writeLines(provenance, file.path(output_dir, "provenance.json"), useBytes = TRUE)
}

main()
