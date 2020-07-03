# Make data frame of metadata for hg38-aligned samples from
# https://doi.org/10.1101/120386 aka "BrainEpigenome samples".
# Peter Hickey
# 2018-06-25

# Setup ------------------------------------------------------------------------

library(S4Vectors)
library(SummarizedExperiment)

extdir <- "../extdata"

# Read colData from saved objects ----------------------------------------------

# NOTE: There was no formal workflow for creating the metadata and colData for
#       the BrainEpigenome project. So we just load what was used in practice.

CD_sorted <- colData(
  readRDS(
    file.path(
      extdir,
      "flow-sorted-brain-wgbs",
      "objects",
      "BS.fit.small.sorted.somatic.all",
      "se.rds")
  ))
CD_bulk <- colData(
  readRDS(
    file.path(
      extdir,
      "flow-sorted-brain-wgbs",
      "objects",
      "BS.unsorted.fit.small.somatic.all",
      "se.rds")
  ))

# Normalise columns ------------------------------------------------------------

# NOTE: Some columns added for consistency with sorted data.
CD_bulk$NeuN <- "bulk"
CD_bulk$NeuN_color <- "black"
# NOTE: Drop lty column
CD_sorted$lty <- NULL

# Reorder columns --------------------------------------------------------------

col_order <- sort(colnames(CD_bulk))
CD_bulk <- CD_bulk[, col_order]
CD_sorted <- CD_sorted[, col_order]

# Sanity check -----------------------------------------------------------------

stopifnot(identical(colnames(CD_bulk), colnames(CD_sorted)))

# Combine into single DataFrame ------------------------------------------------

CD <- rbind(CD_bulk, CD_sorted)

# Construct unique sample ID ---------------------------------------------------

CD$Sample_ID <- paste0(CD$Individual, "_", CD$Tissue, "_", CD$NeuN)
rownames(CD) <- CD$Sample_ID

# Save colData -----------------------------------------------------------------

saveRDS(CD, "../objects/colData.BrainEpigenome.rds")
