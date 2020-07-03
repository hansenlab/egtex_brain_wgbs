# Prepare hg38-aligned samples from https://doi.org/10.1101/120386 aka
# "BrainEpigenome samples" for caling CG-DMRs.
# Peter Hickey
# 2019-02-12

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)

# Load BSseq object and perform basic filtering --------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "BrainEpigenome_hg38",
    "BSseq",
    "BrainEpigenome_hg38.mCG.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "BrainEpigenome_hg38",
    "BSseq",
    "BrainEpigenome_hg38.mCG.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

# NOTE: Only using Neun+ samples
bsseq_filtered <- bsseq_filtered[, bsseq_filtered$NeuN == "pos"]

# General CG-DMRs --------------------------------------------------------------

# Our general DMRs consider five groups:
# 1) BA9
# 2) BA24
# 3) HC
# 4) NAcc

bsseq <- bsseq_filtered
assays(bsseq) <- endoapply(assays(bsseq, withDimnames = FALSE), as.matrix)

bsseq$group <- factor(bsseq$Tissue)

zero_coverage <- rowAnys(x = getCoverage(bsseq), value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 88% of loci
#       (remove 3,414,326 loci).
#       n_low_coverage <- colTabulates(x = getCoverage(bsseq), values = c(0L, 1L, 2L))
#       shows that no samples have ~2-fold more zero coverage loci than the
#       median sample.
#       There are 0 samples that are boxplot outliers based on
#       boxplot.stats(n_low_coverage).
#       nn_low_coverage <- rowTabulates(x = getCoverage(bsseq), values = 0L)
#       prop.table(table(nn_low_coverage[nn_low_coverage > 0])) shows that 39%
#       of zero coverage loci are due to 1 sample, 13% due to 2 sample, 7% due
#       to 3 samples, <5% for all other values (except for 7% of loci with zero
#       coverage in all samples).
#       We would re-gain 1,334,940 loci if we relaxed our filter to allow up to
#       1 sample to have zero coverage for a locus.

bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "BrainEpigenome_hg38.mCG.small_smooth.general_DMRs"),
  # NOTE: Row-wise chunking such that a chunk of doubles in-memory will be
  #       100 Mb.
  chunkdim = makeCappedVolumeBox(
    maxvol = 1e8 / 8,
    maxdim = dim(bsseq),
    shape = "last-dim-grows-first"),
  verbose = TRUE)

file.copy(
  from = file.path(
    tempdir(),
    "BrainEpigenome_hg38.mCG.small_smooth.general_DMRs"),
  to = here(
    "DMRs",
    "extdata",
    "BrainEpigenome_hg38",
    "BSseq",
    "BrainEpigenome_hg38.mCG.small_smooth"),
  recursive = TRUE)
