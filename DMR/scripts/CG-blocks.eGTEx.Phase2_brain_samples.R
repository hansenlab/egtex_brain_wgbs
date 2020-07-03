# Call CG-blocks for hg38-aligned Phase2 brain samples aka "eGTEx Phase2 brain
# samples".
# Peter Hickey
# 2019-02-19

# NOTE: Try to run with as many cores as possible, e.g.,
#       qrsh -l mem_free=6G,h_vmem=6.1G -pe local 64 bash

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)

source(here("DMRs", "scripts", "DMR_functions.R"))

# General CG-blocks ------------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.general_blocks",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.general_blocks",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(666)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 2 ^ 2,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 1000,
  type = "blocks",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  fstat_pipeline,
  here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CG-blocks",
    "eGTEx.Phase2_brain_samples.general_CG-blocks.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50 &
                                 fstat_pipeline$dmrs$maxDiff > 0.1, ],
  file = here("DMRs", "objects", "general_CG-blocks.rds"))

# Overall CG-blocks ------------------------------------------------------------

# NOTE: These are like the general DMRs but treating each tissue as its own
#       group.

bsseq$group <- bsseq$Tissue
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20190220)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 2 ^ 2,
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
  fstat_pipeline,
  here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CG-blocks",
    "eGTEx.Phase2_brain_samples.overall_CG-blocks.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50 &
                                 fstat_pipeline$dmrs$maxDiff > 0.1, ],
  file = here("DMRs", "objects", "overall_CG-blocks.rds"))

# Hippocampus CG-blocks --------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.hippocampus_blocks",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.hippocampus_blocks",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(888)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses t-stat, not F-stat, so don't square the cutoff.
  cutoff = 2,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 1000,
  type = "blocks",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  fstat_pipeline,
  here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CG-blocks",
    "eGTEx.Phase2_brain_samples.hippocampus_CG-blocks.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50 &
                                 fstat_pipeline$dmrs$maxDiff > 0.1, ],
  file = here("DMRs", "objects", "hippocampus_CG-blocks.rds"))

# Basal ganglia CG-blocks ------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.basal_ganglia_blocks",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth.basal_ganglia_blocks",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20010707)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 2 ^ 2,
  fac = bsseq$group,
  nperm = 1000,
  coef = NULL,
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 1000,
  type = "blocks",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  chunksize = chunkdim(assay(bsseq, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

saveRDS(
  fstat_pipeline,
  here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "CG-blocks",
    "eGTEx.Phase2_brain_samples.basal_ganglia_CG-blocks.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50 &
                                 fstat_pipeline$dmrs$maxDiff > 0.1, ],
  file = here("DMRs", "objects", "basal_ganglia_CG-blocks.rds"))
