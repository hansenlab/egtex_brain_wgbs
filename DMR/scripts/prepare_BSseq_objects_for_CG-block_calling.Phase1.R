# Prepare hg38-aligned Phase1 samples aka "eGTEx Phase1 samples" BSseq objects
# for calling CG-blocks.
# Peter Hickey
# 2019-02-06

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)

# Load BSseq objects and perform basic filtering -------------------------------

# Brain samples
bsseq_brain_tempdir <- tempfile("BSseq_brain")
dir.create(bsseq_brain_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Bulk_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase1_brain_samples.mCG.large_smooth",
    "se.rds"),
  to = bsseq_brain_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Bulk_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase1_brain_samples.mCG.large_smooth",
    "assays.h5"),
  to = bsseq_brain_tempdir)
bsseq_brain_unfiltered <- loadHDF5SummarizedExperiment(bsseq_brain_tempdir)
bsseq_brain_filtered <- keepSeqlevels(
  x = bsseq_brain_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

# Non-brain samples
bsseq_nonbrain_tempdir <- tempfile("BSseq_nonbrain")
dir.create(bsseq_nonbrain_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "GTEx_nonbrain_hg38",
    "BSseq",
    "eGTEx.Phase1_non-brain_samples.mCG.large_smooth",
    "se.rds"),
  to = bsseq_nonbrain_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "GTEx_nonbrain_hg38",
    "BSseq",
    "eGTEx.Phase1_non-brain_samples.mCG.large_smooth",
    "assays.h5"),
  to = bsseq_nonbrain_tempdir)
bsseq_nonbrain_unfiltered <- loadHDF5SummarizedExperiment(
  bsseq_nonbrain_tempdir)

# Combine
# NOTE: Temporarily disable validity checks to do a fast cbind(); only safe
#       because I know the objects have the identical rowRanges and compatible
#       colData.
S4Vectors:::disableValidity(TRUE)
bsseq_unfiltered <- cbind(bsseq_brain_unfiltered, bsseq_nonbrain_unfiltered)
S4Vectors:::disableValidity(FALSE)

bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")
load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

# General CG-blocks --------------------------------------------------------------

# Our general blocks consider eight groups:
# 1) Brain
# 2) Breast - Mammary Tissue
# 3) Esophagus - Muscularis
# 4) Heart - Left Ventricle
# 5) Lung
# 6) Muscle - Skeletal
# 7) Skin - Sun Exposed (Lower leg)
# 8) Thyroid

bsseq <- bsseq_filtered
assays(bsseq) <- endoapply(assays(bsseq, withDimnames = FALSE), as.matrix)

bsseq$group <- factor(
  x = ifelse(
    grepl("Brain", bsseq_filtered$Tissue), "Brain",
  bsseq_filtered$Tissue))

zero_coverage <- rowAnys(x = getCoverage(bsseq), value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 57% of loci
#       (remove 11,833,083 loci).
#       n_low_coverage <- colTabulates(x = getCoverage(bsseq), values = c(0L, 1L, 2L))
#       shows that GTEX-XUJ4-0926-SM-9VXYD has ~2-fold more zero coverage loci
#       than the median sample.
#       There are 4 samples that are boxplot outliers based on
#       boxplot.stats(n_low_coverage):
#       GTEX-XUJ4-0926-SM-9VXYD, GTEX-13OW8-3026-SM-9VXYK,
#       GTEX-XV7Q-1026-SM-9VXYC, GTEX-XV7Q-2926-SM-9VXYB.
#       nn_low_coverage <- rowTabulates(x = getCoverage(bsseq), values = 0L)
#       prop.table(table(nn_low_coverage[nn_low_coverage > 0])) shows that 29%
#       of zero coverage loci are due to 1 sample, 14% due to 2 sample, 9% due
#       to 3 samples, <6% for all other values.
#       We would re-gain 3,432,566 loci if we relaxed our filter to allow up to
#       1 sample to have zero coverage for a locus.
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase1.mCG.large_smooth.general_blocks"),
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
    "eGTEx.Phase1.mCG.large_smooth.general_blocks"),
  # NOTE: Saving to 'non-brain' directory although this object contains brain
  #       samples as well.
  to = here(
    "DMRs",
    "extdata",
    "GTEx_nonbrain_hg38",
    "BSseq",
    "eGTEx.Phase1.mCG.large_smooth"),
  recursive = TRUE)
