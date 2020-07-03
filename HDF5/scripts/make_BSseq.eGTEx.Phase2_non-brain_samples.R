# Make HDF5-backed BSseq objects for hg38-aligned Phase2 non-brain samples aka
# "eGTEx Phase2 non-brain samples".
# Peter Hickey
# 2018-08-15

# TODO: Currently have to do `export TMPDIR="/sge_scratch/personal/phickey"`
#       at the command line before starting R

# Setup ------------------------------------------------------------------------

# NOTE: Using 'refactor' branch of bsseq
# TODO: Remove comment when merged into master and be sure to update code as
#       needed.
library(bsseq)
library(HDF5Array)
library(BiocParallel)
library(rtracklayer)
library(here)

CD <- readRDS(here("metadata", "objects", "colData.eGTEx.Phase2.rds"))
# NOTE: Only care about non-brain samples.
CD <- CD[!grepl("Brain", CD[["Tissue"]]), ]

# Import FASTA file of reference genome, which will be used to construct 'loci'
# argument of read.bismark().
# NOTE: Cannot do `seqinfo(fasta) <- seqinfo` because there is no dedicated
#       seqinfo<-,DNAStringSet-method.
fasta <- import(here("HDF5", "extdata", "hg38", "hg38_noALT_with_lambda.fasta"))
seqinfo <- seqinfo(fasta)
seqlevels(seqinfo) <- sapply(strsplit(names(fasta), "\\s+"), "[[", 1)
genome(seqinfo) <- ifelse(
  seqlevels(seqinfo) == "chrEBV",
  "Human_herpesvirus_4",
  ifelse(
    seqlevels(seqinfo) == "lambda",
    "ENA|J02459|J02459.1",
    "hg38"))
isCircular(seqinfo) <- ifelse(
  seqlevels(seqinfo) %in% c("chrM", "chrEBV"),
  TRUE,
  FALSE)
names(fasta) <- seqlevels(seqinfo)
# NOTE: Order 'seqinfo' and 'fasta' using sortSeqlevels()
seqinfo <- sortSeqlevels(seqinfo)
fasta <- fasta[seqlevels(seqinfo)]

# mCG --------------------------------------------------------------------------

out_name <- "eGTEx.Phase2_non-brain_samples.mCG"
scratch_dir <- file.path("/sge_scratch", "personal", "phickey", out_name)

cg <- findCytosines(fasta, "CG", seqlevels(fasta))
# NOTE: Need to fix up seqinfo
seqinfo(cg) <- seqinfo

files <- here(
  "HDF5",
  "extdata",
  "Sorted_GTEx_Brain_hg38",
  CD[["HiSeq_Run"]],
  CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  paste0(
    CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
    ".merged.nsorted.deduplicated.CpG_report.txt.CpG_report.txt.gz"))
bsseq <- read.bismark(
  files = files,
  loci = cg,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Collapsing loci by strand.
  strandCollapse = TRUE,
  BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

# Collapse samples sequenced across multiple runs.
# NOTE: There aren't any samples sequenced over multiple runs
#       (i.e. length(unique(CD$Sample_ID)) == nrow(CD)), so this just returns
#       the input object.
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM =  MulticoreParam(workers = 30, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

# mCG small smooth -------------------------------------------------------------

# NOTE: Copy unsmoothed HDF5-backed BSseq object to a new directory in
#       preparation for smoothing.
dir.create(paste0(scratch_dir, ".small_smooth"))
file.copy(
  from = list.files(scratch_dir, full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth"))
bsseq <- loadHDF5SummarizedExperiment(paste0(scratch_dir, ".small_smooth"))
# NOTE: Used `ns = 20` in Nature Neuroscience paper, but decided to revert to
#       `ns = 70` (the default) for this work
#       (https://jhu-genomics.slack.com/archives/GB1UH5DB2/p1528735628000421).
bsseq <- BSmooth(
  BSseq = bsseq,
  ns = 70,
  h = 1000,
  maxGap = 10^8,
  keep.se = FALSE,
  BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Sorted_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCG large smooth -------------------------------------------------------------

dir.create(paste0(scratch_dir, ".large_smooth"))
file.copy(
  from = list.files(scratch_dir, full.names = TRUE),
  to = paste0(scratch_dir, ".large_smooth"))
bsseq <- loadHDF5SummarizedExperiment(paste0(scratch_dir, ".large_smooth"))
# NOTE: Using same smoothing parameters as in Nature Neuroscience paper.
bsseq <- BSmooth(
  BSseq = bsseq,
  ns = 500,
  h = 20000,
  maxGap = 10^8,
  keep.se = FALSE,
  BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".large_smooth"),
  to = here("HDF5", "extdata", "Sorted_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCA (+) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase2_non-brain_samples.mCA_pos"
scratch_dir <- file.path("/sge_scratch", "personal", "phickey", out_name)
# NOTE: dir must not yet exist.
stopifnot(!dir.exists(scratch_dir))

ca_pos <- findCytosines(fasta, "CA", seqlevels(fasta), "+")
# NOTE: Need to fix up seqinfo
seqinfo(ca_pos) <- seqinfo

files <- vapply(
  X = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  FUN = grep,
  FUN.VALUE = character(1L),
  x = list.files(
    path = here("HDF5", "extdata", "Sorted_GTEx_Brain_hg38", "nonCG"),
    pattern = glob2rx("*.CX_report.txt.gz"),
    full.names = TRUE),
  value = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = ca_pos,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 1, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

# Collapse samples sequenced across multiple runs.
# NOTE: There aren't any samples sequenced over multiple runs
#       (i.e. length(unique(CD$Sample_ID)) == nrow(CD)), so this just returns
#       the input object.
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM =  MulticoreParam(workers = 30, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

# mCA (+) small-ish smooth -----------------------------------------------------

# NOTE: Copy unsmoothed HDF5-backed BSseq object to a new directory in
#       preparation for smoothing.
dir.create(paste0(scratch_dir, ".smallish_smooth"))
file.copy(
  from = list.files(dirname(path(assay(bsseq))), full.names = TRUE),
  to = paste0(scratch_dir, ".smallish_smooth"))
bsseq <- loadHDF5SummarizedExperiment(paste0(scratch_dir, ".smallish_smooth"))
# NOTE: Using same smoothing parameters as in Nature Neuroscience paper.
bsseq <- BSmooth(
  BSseq = bsseq,
  ns = 200,
  h = 3000,
  maxGap = 10^8,
  keep.se = FALSE,
  BPPARAM = MulticoreParam(workers = 10, progressbar = TRUE),
  verbose = TRUE)
# TODO: The above smooth yielded an error.
# "Error in serialize(data, node$con, xdr = FALSE) : ignoring SIGPIPE signal
#  Error: failed to stop ‘SOCKcluster’ cluster: error writing to connection"
