# Call CG-DMRs for hg38-aligned Phase2 brain samples aka "eGTEx Phase2 brain
# samples".
# Peter Hickey
# 2018-09-30

# NOTE: Try to run with as many cores as possible, e.g.,
#       qrsh -l mem_free=6G,h_vmem=6.1G -pe local 64 bash

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)

source(here("DMRs", "scripts", "DMR_functions.R"))

# General CG-DMRs --------------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth.general_DMRs",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth.general_DMRs",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20180722)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 4.6 ^ 2,
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
    "CG-DMRs",
    "eGTEx.Phase2_brain_samples.general_CG-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "general_CG-DMRs.rds"))

# Overall CG-DMRs --------------------------------------------------------------

# NOTE: These are like the general DMRs but treating each tissue as its own
#       group.

bsseq$group <- bsseq$Tissue
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20190221)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 4.6 ^ 2,
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
    "CG-DMRs",
    "eGTEx.Phase2_brain_samples.overall_CG-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "overall_CG-DMRs.rds"))

# Hippocampus CG-DMRs ----------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth.hippocampus_DMRs",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth.hippocampus_DMRs",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20180930)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses t-stat, not F-stat, so don't square the cutoff.
  cutoff = 4.6,
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
    "CG-DMRs",
    "eGTEx.Phase2_brain_samples.hippocampus_CG-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "hippocampus_CG-DMRs.rds"))

# TODO: Need to check if these plots were made (if not, decide if they should
#       be cause they're gonna be huge ... perhaps just plot top-X and bottom-X)
# pdf(
#   here(
#     "CG-DMRs",
#     "figures",
#     "eGTEx.Phase2_brain_samples.hippocampus_CG-DMRs.pdf"))
# plot.new()
# legend("center", legend = c("HC1", "HC2"), col = bsseq$Brain_color, lty = c(1, 2))
# plotManyRegions(
#   BSseq = bsseq,
#   regions = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
#   extend = 20000,
#   addRegions = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
#   col = bsseq$Brain_color,
#   lty = ifelse(bsseq$group == "HC1", 1, 2))
# dev.off()

# Basal ganglia CG-DMRs --------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth.basal_ganglia_DMRs",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth.basal_ganglia_DMRs",
    "assays.h5"),
  to = bsseq_tempdir)

bsseq <- loadHDF5SummarizedExperiment(bsseq_tempdir)
design <- model.matrix(~ group, colData(bsseq))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
rownames(contrasts) <- colnames(design)

set.seed(20180707)
fstat_pipeline <- fstat.pipeline(
  BSseq = bsseq,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 4.6 ^ 2,
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
    "CG-DMRs",
    "eGTEx.Phase2_brain_samples.basal_ganglia_CG-DMRs.fstat_pipeline.rds"))
saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "basal_ganglia_CG-DMRs.rds"))

# pdf(
#   here(
#     "CG-DMRs",
#     "figures",
#     "eGTEx.Phase2_brain_samples.basal_ganglia_CG-DMRs.pdf"))
# plot.new()
# legend(
#   "center",
#   legend = unique(colData(bsseq)[, c("Brain_color", "group")])[["Brain_color"]],
#   col = unique(colData(bsseq)[, c("Brain_color", "group")])[["group"]],
#   lty = 1)
# plotManyRegions(
#   BSseq = bsseq,
#   regions = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
#   extend = 20000,
#   addRegions = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
#   col = bsseq$Brain_color,
#   lty = 1)
# dev.off()
