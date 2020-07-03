# Prepare eGTEx WGBS data for use in a genome browser. Only autosomal data is
# currently supported.
# Peter Hickey
# 2019-11-07

# Setup ------------------------------------------------------------------------

library(bsseq)
library(HDF5Array)
library(here)
library(DelayedMatrixStats)
library(rtracklayer)
library(dplyr)
library(janitor)
library(BSgenome.Hsapiens.UCSC.hg38)

# NOTE: Don't increase mc.cores (export.bw() craps out)
options("mc.cores" = 8)

DelayedArray::setAutoBlockSize(.Machine$integer.max)
DelayedArray:::set_verbose_block_processing(TRUE)
DelayedArray:::set_verbose_block_processing(TRUE)

# NOTE: chrom.sizes file needed for constructing a bigBED file from a BED file
chrom.sizes_file <- here("genome_browser", "data", "hg38.chrom.sizes")
si <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
write.table(
  x = as.data.frame(si)[, 1, drop = FALSE],
  file = chrom.sizes_file,
  col.names = FALSE,
  row.names = TRUE,
  quote = FALSE,
  sep = "\t")
# NOTE: BED to bigBED requies the chromosomes are sorted lexicographically
si_bed <- si[sort(seqnames(si))]

# Functions --------------------------------------------------------------------

makeBigWig <- function(bsseq, groups, name, mc.cores) {
  mclapply(levels(groups), function(group) {
    message("Converting ", group)
    j <- which(groups == group)
    gr <- rowRanges(bsseq)
    gr$score <- rowMeans2(x = getMeth(bsseq), cols = j, na.rm = TRUE)
    # NOTE: Replace NA/NaN with 0 (not sure how they end up there ...)
    gr$score[is.na(gr$score)] <- 0
    export(
      object = gr,
      con = here(
        "genome_browser",
        "extdata",
        "Sorted_GTEx_Brain_hg38",
        "bigWig",
        paste0(group, ".", name, ".bw")))
  }, mc.cores = mc.cores)
}

# mCG small smooth -------------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCG.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCG.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCG.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCG.small_smooth",
  mc.cores = getOption("mc.cores"))

# mCG large smooth -------------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCG.large_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCG.large_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCG.large_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCG.large_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCG.large_smooth",
  mc.cores = getOption("mc.cores"))

# mCA (+) small smooth ---------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_pos.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCA_pos.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCA_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCA_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCA_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCA_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

# mCA (-) small smooth ---------------------------------------------------------

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

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCA_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCA_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCA_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCA_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

# mCT (+) small smooth ---------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_pos.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_pos.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCT_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCT_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCT_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCT_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

# mCT (-) small smooth ---------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCT_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCT_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCT_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCT_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

# mCC (+) small smooth ---------------------------------------------------------

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCC_pos.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCC_pos.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCC_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCC_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCC_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCC_pos.small_smooth",
  mc.cores = getOption("mc.cores"))

# mCC (-) small smooth ---------------------------------------------------------

# TODO: UP TO HERE for sample-level bigWigs

bsseq_tempdir <- tempfile("BSseq")
dir.create(bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "se.rds"),
  to = bsseq_tempdir)
file.copy(
  from = here(
    "DMRs",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BSseq",
    "eGTEx.Phase2_brain_samples.mCT_neg.small_smooth",
    "assays.h5"),
  to = bsseq_tempdir)
bsseq_unfiltered <- loadHDF5SummarizedExperiment(bsseq_tempdir)

# NOTE: Only support chr1-chr22.
bsseq_filtered <- keepSeqlevels(
  x = bsseq_unfiltered,
  value = paste0("chr", 1:22),
  pruning.mode = "coarse")

load(here("miscellaneous_objects", "QC_failures.rda"))
bsseq_filtered <- bsseq_filtered[, !colnames(bsseq_filtered) %in%
                                   genotype_failures]
bsseq_filtered <- bsseq_filtered[, !bsseq_filtered$Sample_ID %in%
                                   mCG_PCA_outliers]

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Brain_Frontal_Cortex_BA9",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Brain_Anterior_cingulate_cortex_BA24",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~
        "Brain_Hippocampus",
      bsseq_filtered$Tissue == "Brain - Amygdala" ~
        "Brain_Amygdala",
      bsseq_filtered$Tissue == "Brain - Hypothalamus" ~
        "Brain_Hypothalamus",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Brain_Nucleus_accumbens_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Brain_Putamen_basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Brain_Caudate_basal_ganglia")),
  name = "mCC_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Tissue == "Brain - Frontal Cortex (BA9)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Anterior cingulate cortex (BA24)" ~
        "Cortical",
      bsseq_filtered$Tissue == "Brain - Nucleus accumbens (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Putamen (basal ganglia)" ~
        "Basal_ganglia",
      bsseq_filtered$Tissue == "Brain - Caudate (basal ganglia)" ~
        "Basal_ganglia")),
  name = "mCC_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(
    case_when(
      bsseq_filtered$Sample_ID %in% c(
        "SM-AFUMH", "SM-AFUKV", "SM-AFUNS", "SM-AFUMM", "SM-AFUK5", "SM-AFUKJ",
        "SM-AFUJL") ~ "HC1",
      bsseq_filtered$Tissue == "Brain - Hippocampus" ~ "HC2")),
  name = "mCC_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

makeBigWig(
  bsseq = bsseq_filtered,
  groups = factor(colnames(bsseq_filtered)),
  name = "mCC_neg.small_smooth",
  mc.cores = getOption("mc.cores"))

# DMRs ----------------------------------------------------------------------

# TODO: These DMRs are named according to (arbitrary) default row names on
#       the data.frame; should these be removed? (Check BrainEpigenome).

files <- list.files(
  path = here("DMRs", "objects"),
  pattern = glob2rx("*.rds"),
  full.names = TRUE)
names(files) <- gsub(".rds", "", basename(files))
list_of_dmrs <- lapply(files, readRDS)
list_of_dmrs <- lapply(list_of_dmrs, as, "GRanges")
lapply(names(list_of_dmrs), function(n) {
  message(n)
  dmrs <- list_of_dmrs[[n]]
  seqlevels(dmrs) <- seqlevels(intersect(si_bed, seqinfo(dmrs)))
  seqinfo(dmrs) <- merge(seqinfo(dmrs), si_bed)
  bed_file <- here(
    "genome_browser",
    "extdata",
    "Sorted_GTEx_Brain_hg38",
    "BED",
    paste0(n, ".bed"))
  export(
    object = sort(dmrs),
    con = bed_file)
  big_bed_file <- sub("\\.bed", "\\.bb", sub("BED", "bigBED", bed_file))
  system(
    paste0("bedToBigBed ", bed_file, " ", chrom.sizes_file, " ", big_bed_file))
})
