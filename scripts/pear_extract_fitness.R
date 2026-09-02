# Extract the PEAR (Zhang 2022, blaCTX-M-14 DMS) per-variant fitness values out of the authors'
# .RData workspaces into plain TSV, so the G6 assay-degeneracy gate can be screened on real numbers.
#
# WHY R AT ALL. The repo ships no plain-text data -- only two .RData workspaces. Python's `rdata` and
# `pyreadr` both fail, and NOT merely at conversion: `rdata.parser.parse_file` itself raises on a
# WEAKREF nested inside serialized bytecode (a closure carried by the ggplot objects). Parsing is
# sequential, so the failure cannot be stepped over without risking silent byte-desync -- which would
# yield numbers that LOOK like data. R reads its own format exactly; that is the whole reason to use it.
#
# Outputs go to D: (this host's project-data drive), never into the repo.

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "D:/dna_decode_cache/pear/extracted"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- c("D:/dna_decode_cache/pear/CTXM-14/Data_for_Figure2.RData",
           "D:/dna_decode_cache/pear/CTXM-14/Data_for_Figure3.RData")

for (f in files) {
  cat("==== ", f, "\n", sep = "")
  e <- new.env()
  nm <- load(f, envir = e)
  for (n in nm) {
    o <- get(n, envir = e)
    cls <- paste(class(o), collapse = "/")
    if (is.data.frame(o)) {
      cat(sprintf("  DF   %-34s %6d x %2d  cols: %s\n", n, nrow(o), ncol(o),
                  paste(names(o), collapse = ",")))
      write.table(o, file.path(outdir, paste0(n, ".tsv")), sep = "\t",
                  row.names = FALSE, quote = FALSE)
    } else if (inherits(o, "ggplot") || (is.list(o) && !is.null(o$data))) {
      # A ggplot object carries the data frame it was built from in $data. These plots are
      # genome-wide tile maps, so that slot is the FULL scan, not a figure-cropped subset.
      d <- o$data
      if (is.data.frame(d)) {
        cat(sprintf("  gg$  %-34s %6d x %2d  cols: %s\n", n, nrow(d), ncol(d),
                    paste(names(d), collapse = ",")))
        write.table(d, file.path(outdir, paste0(n, "__data.tsv")), sep = "\t",
                    row.names = FALSE, quote = FALSE)
      } else {
        cat(sprintf("  obj  %-34s %s ($data is %s)\n", n, cls, paste(class(d), collapse = "/")))
      }
    } else {
      cat(sprintf("  obj  %-34s %s\n", n, cls))
    }
  }
}
cat("\nwrote TSVs to ", outdir, "\n", sep = "")
