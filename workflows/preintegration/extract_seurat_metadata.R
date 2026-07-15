suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript extract_seurat_metadata.R <input_rds_or_rds.gz> <output_csv_gz>')
infile <- args[[1]]
outfile <- args[[2]]

read_obj <- function(path) {
  # GEO sometimes stores nested gzip streams (*.rds.gz that unzip to *.rds.gz).
  # Try plain, one gzip layer, then two gzip layers.
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(obj)) return(obj)
  if (grepl('\\.gz$', path, ignore.case = TRUE)) {
    con1 <- gzfile(path, open = 'rb')
    on.exit(close(con1), add = TRUE)
    obj1 <- tryCatch(readRDS(con1), error = function(e) NULL)
    if (!is.null(obj1)) return(obj1)
    con2 <- gzcon(gzfile(path, open = 'rb'))
    on.exit(close(con2), add = TRUE)
    obj2 <- tryCatch(readRDS(con2), error = function(e) NULL)
    if (!is.null(obj2)) return(obj2)
  }
  stop(sprintf('Failed to read RDS: %s', path))
}

pick_seurat <- function(x) {
  if (inherits(x, 'Seurat')) return(x)
  if (is.list(x)) {
    for (nm in names(x)) {
      if (inherits(x[[nm]], 'Seurat')) return(x[[nm]])
    }
  }
  return(NULL)
}

obj <- read_obj(infile)
seu <- pick_seurat(obj)
if (is.null(seu)) stop(sprintf('No Seurat object found in %s (class: %s)', infile, paste(class(obj), collapse = ',')))

md <- seu@meta.data
md$cell_barcode <- rownames(md)
md <- md[, c('cell_barcode', setdiff(colnames(md), 'cell_barcode')), drop = FALSE]

outdir <- dirname(outfile)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
con_out <- gzfile(outfile, "wt")
on.exit(close(con_out), add = TRUE)
write.csv(md, con_out, row.names = FALSE)
close(con_out)

cat(sprintf('ok\tinput=%s\trows=%d\tcols=%d\tout=%s\n', infile, nrow(md), ncol(md), outfile))
