# Prepare hg38-aligned Phase2 brain samples aka "eGTEx Phase2 brain samples"
# BSseq objects for calling CG-blocks.
# Peter Hickey
# 2018-09-30

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)

# Load BSseq object and perform basic filtering --------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "CG-blocks",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "CG-blocks",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

# General CG-blocks ------------------------------------------------------------

# Our general blocks consider five groups:
# 1) Cortical: frontal cortex (BA9) and anterior cingulate cortex (ACC)
# 2) HC
# 3) Amygdala
# 4) Hypothalamus
# 5) Basal ganglia: NAcc, Putamen, and Caudate*

bsseq <- bsseq_filtered
assays(bsseq) <- endoapply(assays(bsseq, withDimnames = FALSE), as.matrix)

bsseq$group <- factor(
  dplyr::case_when(
    bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~ "Cortical",
    bsseq_filtered$Tissue ==
      "Brain - Anterior cingulate cortex (BA24)" ~ "Cortical",
    bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "Hippocampus",
    bsseq_filtered$Tissue == "Brain - Amygdala" ~ "Amygdala",
    bsseq_filtered$Tissue == "Brain - Hypothalamus" ~ "Hypothalamus",
    bsseq_filtered$Tissue ==
      "Brain - Nucleus accumbens (basal ganglia)" ~ "Basal ganglia",
    bsseq_filtered$Tissue ==
      "Brain - Putamen (basal ganglia)" ~ "Basal ganglia",
    bsseq_filtered$Tissue ==
      "Brain - Caudate (basal ganglia)" ~ "Basal ganglia"),
  levels = c(
    "Cortical", "Hippocampus", "Amygdala", "Hypothalamus", "Basal ganglia")
)

zero_coverage <- rowAnys(x = getCoverage(bsseq), value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 77% of loci
#       (remove 6,405,781 loci).
#       n_low_coverage <- colTabulates(x = getCoverage(bsseq), values = c(0L, 1L, 2L))
#       shows that GTEX-RU72-0011-R10a-SM-AFUO8 and
#       GTEX-ZVT3-0011-R10b-SM-AFUP6 have ~2-fold more zero coverage loci than
#       the median sample.
#       There are 7 samples that are boxplot outliers based on
#       boxplot.stats(n_low_coverage):
#       GTEX-RU72-0011-R10a-SM-AFUO8, GTEX-ZVT3-0011-R10b-SM-AFUP6,
#       GTEX-13NZA-0011-R8b-SM-AFUKO, GTEX-T2IS-0011-R5a-SM-AFUOT,
#       GTEX-WZTO-0011-R1b-SM-AFUOR, GTEX-RU72-0011-R6a-SM-AFUOA,
#       GTEX-WZTO-0011-R6b-SM-AFUOQ.
#       nn_low_coverage <- rowTabulates(x = getCoverage(bsseq), values = 0L)
#       prop.table(table(nn_low_coverage[nn_low_coverage > 0])) shows that 35%
#       of zero coverage loci are due to 1 sample, 14% due to 2 sample, 8% due
#       to 3 samples, <5% for all other values.
#       We would re-gain 2,257,972 loci if we relaxed our filter to allow up to
#       1 sample to have zero coverage for a locus.
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.general_blocks"),
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
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.general_blocks"),
  to = here(
    "CG-blocks",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth"),
  recursive = TRUE)

# Hippocampus CG-blocks --------------------------------------------------------

# Two groups of hippocampus samples:
# 1) HC1: SM-AFUMH, SM-AFUKV, SM-AFUNS, SM-AFUMM, SM-AFUK5, SM-AFUKJ, SM-AFUJL
# 2) HC2: The rest of the HC samples

bsseq <- bsseq_filtered[, bsseq_filtered$Tissue == "Brain - Hippocampus"]
assays(bsseq) <- endoapply(assays(bsseq, withDimnames = FALSE), as.matrix)

bsseq$group <- factor(
  dplyr::case_when(
    bsseq$Sample_ID %in% c("SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM",
                           "SM-AFUK5", "SM-AFUKJ", "SM-AFUJL") ~ "HC1",
    TRUE ~ "HC2"),
  levels = c("HC1", "HC2")
)

zero_coverage <- rowAnys(x = getCoverage(bsseq), value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 89% of loci
#       (remove 3,072,602 loci).
#       n_low_coverage <- colTabulates(x = getCoverage(bsseq), values = c(0L, 1L, 2L))
#       shows that no samples has > 1.3-fold more zero coverage loci than the
#       median sample.
#       There are no samples that are boxplot outliers based on
#       boxplot.stats(n_low_coverage).
#       nn_low_coverage <- rowTabulates(x = getCoverage(bsseq), values = 0L)
#       prop.table(table(nn_low_coverage[nn_low_coverage > 0])) shows that 45%
#       of zero coverage loci are due to 1 sample, 14% due to 2 sample, 6% due
#       to 3 samples, <5% for all other values.
#       We would re-gain 1,378,210 loci if we relaxed our filter to allow up to
#       1 sample to have zero coverage for a locus.
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.hippocampus_blocks"),
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
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.hippocampus_blocks"),
  to = here(
    "CG-blocks",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth"),
  recursive = TRUE)

# Basal ganglia CG-blocks ------------------------------------------------------

# Three groups of samples
# 1) NAcc
# 2) Putamen
# 3) Caudate

bsseq <- bsseq_filtered[, bsseq_filtered$Tissue %in% c(
  "Brain - Nucleus accumbens (basal ganglia)",
  "Brain - Putamen (basal ganglia)",
  "Brain - Caudate (basal ganglia)")]
assays(bsseq) <- endoapply(assays(bsseq, withDimnames = FALSE), as.matrix)

bsseq$group <- factor(
  dplyr::case_when(
    bsseq$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~ "NAcc",
    bsseq$Tissue == "Brain - Putamen (basal ganglia)" ~ "Putamen",
    bsseq$Tissue == "Brain - Caudate (basal ganglia)" ~ "Caudate"),
  levels = c("NAcc", "Putamen", "Caudate")
)

zero_coverage <- rowAnys(x = getCoverage(bsseq), value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 84% of loci
#       (remove 4,430,237 loci).
#       n_low_coverage <- colTabulates(x = getCoverage(bsseq), values = c(0L, 1L, 2L))
#       shows that no samples has > 1.5-fold more zero coverage loci than the
#       median sample. There are
#       4 samples that are boxplot (upper) outliers based on n_low_coverage
#       (boxplot.stats(n_low_coverage)):
#       nn_low_coverage <- rowTabulates(x = getCoverage(bsseq), values = 0L)
#       prop.table(table(nn_low_coverage[nn_low_coverage > 0])) shows that 40%
#       of zero coverage loci are due to 1 sample, 14% due to 2 sample, 7% due
#       to 3 samples, <5% for all other values.
#       We would re-gain 1,750,312 loci if we relaxed our filter to allow up to
#       1 sample to have zero coverage for a locus.
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.basal_ganglia_blocks"),
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
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.basal_ganglia_blocks"),
  to = here(
    "CG-blocks",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth"),
  recursive = TRUE)
