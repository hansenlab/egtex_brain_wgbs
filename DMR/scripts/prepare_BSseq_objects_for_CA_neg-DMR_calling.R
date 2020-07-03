# Prepare hg38-aligned Phase2 brain samples aka "eGTEx Phase2 brain samples"
# BSseq objects for calling CA_neg-DMRs.
# Peter Hickey
# 2018-11-06

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
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth",
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
# NOTE: Removing samples based on mCG even though this is mCA (-).
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

# General CA_neg-DMRs ----------------------------------------------------------

# Our general DMRs consider five groups:
# 1) Cortical: frontal cortex (BA9) and anterior cingulate cortex (ACA)
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

zero_coverage <- rowAnys(
  x = assay(bsseq, "Cov", withDimnames = FALSE),
  value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 41% of loci
#       (remove 117,024,866 loci).
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth.general_DMRs"),
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
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth.general_DMRs"),
  to = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth"),
  recursive = TRUE)

# Hippocampus CG-DMRs ----------------------------------------------------------

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

zero_coverage <- rowAnys(
  x = assay(bsseq, "Cov", withDimnames = FALSE),
  value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 66% of loci
#       (remove 68,336,229 loci).
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth.hippocampus_DMRs"),
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
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth.hippocampus_DMRs"),
  to = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth"),
  recursive = TRUE)

# Basal ganglia CG-DMRs --------------------------------------------------------

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

zero_coverage <- rowAnys(
  x = assay(bsseq, "Cov", withDimnames = FALSE),
  value = 0L)
# NOTE: I did some quick checks of how harsh this filter is.
#       prop.table(table(zero_coverage)) shows we retain 54% of loci
#       (remove 92,235,145 loci).
bsseq <- bsseq[!zero_coverage, ]

saveHDF5SummarizedExperiment(
  x = bsseq,
  dir = file.path(
    tempdir(),
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth.basal_ganglia_DMRs"),
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
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth.basal_ganglia_DMRs"),
  to = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_neg.small_smooth"),
  recursive = TRUE)
