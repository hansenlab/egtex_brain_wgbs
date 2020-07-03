# Call CT_neg-DMRs for hg38-aligned Phase2 brain samples aka "eGTEx Phase2
# brain samples".
# Peter Hickey
# 2019-02-25

# NOTE: Try to run with as many cores as possible, e.g.,
#       qrsh -l mem_free=6G,h_vmem=6.1G,cegs2 -pe local 64 bash

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)

source(here("DMRs", "scripts", "DMR_functions.R"))

# General CT_neg-DMRs ----------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth.general_DMRs",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth.general_DMRs",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)

design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20181126)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  # NOTE: Using the same cutoff we used in BrainEpigenome paper for CH-DMRs.
  cutoff = 4 ^ 2,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 300,
  type = "dmrs",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  object = fstat_pipeline,
  file = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CT_neg-DMRs",
    "eGTEx.Phase2_brain_samples.general_CT_neg-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "general_CT_neg-DMRs.rds"))

# Overall CT_neg-DMRs --------------------------------------------------------------

# NOTE: These are like the general DMRs but treating each tissue as its own
#       group.

bsseq$group <- bsseq$Tissue
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20190225)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  # NOTE: Using the same cutoff we used in BrainEpigenome paper for CH-DMRs.
  cutoff = 4 ^ 2,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 300,
  type = "dmrs",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  object = fstat_pipeline,
  file = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CT_neg-DMRs",
    "eGTEx.Phase2_brain_samples.overall_CT_neg-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "overall_CT_neg-DMRs.rds"))

# Hippocampus CT_neg-DMRs ------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth.hippocampus_DMRs",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth.hippocampus_DMRs",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20181130)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses t-stat, not F-stat, so don't square the cutoff.
  cutoff = 4,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 300,
  type = "dmrs",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  object = fstat_pipeline,
  file = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CT_neg-DMRs",
    "eGTEx.Phase2_brain_samples.hippocampus_CT_neg-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "hippocampus_CT_neg-DMRs.rds"))

# Basal ganglia CG-DMRs --------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth.basal_ganglia_DMRs",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth.basal_ganglia_DMRs",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20180730)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 4 ^ 2,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 300,
  type = "dmrs",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  object = fstat_pipeline,
  file = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CT_neg-DMRs",
    "eGTEx.Phase2_brain_samples.basal_ganglia_CT_neg-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "basal_ganglia_CT_neg-DMRs.rds"))
