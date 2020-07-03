# Make HDF5-backed BSseq objects for hg38-aligned Phase1 brain samples aka
# "eGTEx Phase1 brain samples".
# Peter Hickey
# 2018-12-03

# NOTE: Remember to run qrsh with `bash` and `-l h_fsize=500G`

# Setup ------------------------------------------------------------------------

library(bsseq)
library(BiocParallel)
library(rtracklayer)
library(here)
library(HDF5Array)

DelayedArray::setAutoBlockSize(.Machine$integer.max)
DelayedArray:::set_verbose_block_processing(TRUE)
DelayedArray:::set_verbose_block_processing(TRUE)
Sys.setenv("HDF5_USE_FILE_LOCKING" = FALSE)

# Load metadata and reference genome -------------------------------------------

CD <- readRDS(here("metadata", "objects", "colData.eGTEx.Phase1.rds"))
# NOTE: Only care about brain samples.
CD <- CD[grepl("Brain", CD[["Tissue"]]), ]

# Import FASTA file of reference genome, which will be used to construct 'loci'
# argument of read.bismark().
# NOTE: Cannot do `seqinfo(fasta) <- seqinfo` because there is no dedicated
#       seqinfo<-,DNAStringSet-method (see
#       https://github.com/Bioconductor/Biostrings/issues/18).
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
# NOTE: Ensure `seqinfo` and `fasta` are in same order.
seqinfo <- sortSeqlevels(seqinfo)
fasta <- fasta[seqlevels(seqinfo)]


# Setup ------------------------------------------------------------------------

# NOTE: Using 'refactor' branch of bsseq
# TODO: Remove comment when merged into master and be sure to update code as
#       needed.
library(bsseq)
library(HDF5Array)
library(BiocParallel)
library(rtracklayer)
library(here)

CD <- readRDS(here("metadata", "objects", "colData.eGTEx.Phase1.rds"))
# NOTE: Only care about brain samples.
CD <- CD[grepl("Brain", CD[["Tissue"]]), ]

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

out_name <- "eGTEx.Phase1_brain_samples.mCG"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

cg <- findLoci(
  pattern = "CG",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "*")
seqinfo(cg) <- seqinfo

files <- here(
  "HDF5",
  "extdata",
  "Bulk_GTEx_Brain_hg38",
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
  BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")
file.copy(
  from = scratch_dir,
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCG small smooth -------------------------------------------------------------

# NOTE: Copy unsmoothed HDF5-backed BSseq object to a new directory in
#       preparation for smoothing
dir.create(paste0(scratch_dir, ".small_smooth"))
file.copy(
  from = list.files(dirname(path(assay(bsseq))), full.names = TRUE),
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
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCG large smooth -------------------------------------------------------------

dir.create(paste0(scratch_dir, ".large_smooth"))
file.copy(
  from = list.files(dirname(path(assay(bsseq))), full.names = TRUE),
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
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCA (+) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase1_brain_samples.mCA_pos"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

ca_pos <- findLoci(
  pattern = "CA",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "+")
seqinfo(ca_pos) <- seqinfo

files <- mapply(
  function(sample_id, hiseq_run) {
    list.files(
      path = here(
        "HDF5",
        "extdata",
        "Bulk_GTEx_Brain_hg38",
        hiseq_run,
        "nonCG"),
      pattern = glob2rx(paste0(sample_id, "*.CX_report.txt.gz")),
      full.names = TRUE)
  },
  sample_id = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  hiseq_run = CD[["HiSeq_Run"]],
  SIMPLIFY = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = ca_pos,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

dir.create(file.path(scratch_dir, "collapseBSseq"))
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

file.copy(
  from = file.path(scratch_dir, "collapseBSseq"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCA (+) small smooth ---------------------------------------------------------

# NOTE: For unknown reasons it was crashing when trying to smooth the entire
#       object (even though only a subset of data is loaded into memory). So I
#       had to smooth on a per-chromosome and then rbind() the outputs.
# NOTE: Run with 8 cores, 50 Gb each core.

# TODO: Running on compute-065

dir.create(paste0(scratch_dir, ".small_smooth.per_chromosome"))
file.copy(
  from = list.files(
    path = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq", out_name),
    full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth.per_chromosome"))
bsseq <- loadHDF5SummarizedExperiment(
  dir = paste0(scratch_dir, ".small_smooth.per_chromosome"))

list_of_bsseq <- lapply(seqlevels(bsseq), function(seqlevel) {
  message(seqlevel,
          paste0(
            " (",
            match(seqlevel, seqlevels(bsseq)),
            "/",
            length(seqlevels(bsseq)),
            ")"))
  bsseq <- BSmooth(
    BSseq = keepSeqlevels(x = bsseq, value = seqlevel, pruning.mode = "coarse"),
    ns = 200,
    h = 3000,
    maxGap = 10^8,
    keep.se = FALSE,
    BPPARAM = MulticoreParam(
      workers = as.integer(Sys.getenv("NSLOTS")),
      progressbar = TRUE),
    verbose = TRUE)
  saveHDF5SummarizedExperiment(
    x = bsseq,
    dir = file.path(tempdir(), seqlevel),
    verbose = TRUE)
})

saveHDF5SummarizedExperiment(
  x = do.call(rbind, list_of_bsseq),
  dir = paste0(scratch_dir, ".small_smooth"),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCA (-) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase1_brain_samples.mCA_neg"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

ca_neg <- findLoci(
  pattern = "CA",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "-")
seqinfo(ca_neg) <- seqinfo

files <- mapply(
  function(sample_id, hiseq_run) {
    list.files(
      path = here(
        "HDF5",
        "extdata",
        "Bulk_GTEx_Brain_hg38",
        hiseq_run,
        "nonCG"),
      pattern = glob2rx(paste0(sample_id, "*.CX_report.txt.gz")),
      full.names = TRUE)
  },
  sample_id = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  hiseq_run = CD[["HiSeq_Run"]],
  SIMPLIFY = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = ca_neg,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

dir.create(file.path(scratch_dir, "collapseBSseq"))
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

file.copy(
  from = file.path(scratch_dir, "collapseBSseq"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCA (-) small smooth ---------------------------------------------------------

# NOTE: For unknown reasons it was crashing when trying to smooth the entire
#       object (even though only a subset of data is loaded into memory). So I
#       had to smooth on a per-chromosome and then rbind() the outputs.
# NOTE: Run with 8 cores, 50 Gb each core.

# TODO: Running on compute-072

dir.create(paste0(scratch_dir, ".small_smooth.per_chromosome"))
file.copy(
  from = list.files(
    path = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq", out_name),
    full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth.per_chromosome"))
bsseq <- loadHDF5SummarizedExperiment(
  dir = paste0(scratch_dir, ".small_smooth.per_chromosome"))

list_of_bsseq <- lapply(seqlevels(bsseq), function(seqlevel) {
  message(seqlevel,
          paste0(
            " (",
            match(seqlevel, seqlevels(bsseq)),
            "/",
            length(seqlevels(bsseq)),
            ")"))
  bsseq <- BSmooth(
    BSseq = keepSeqlevels(x = bsseq, value = seqlevel, pruning.mode = "coarse"),
    ns = 200,
    h = 3000,
    maxGap = 10^8,
    keep.se = FALSE,
    BPPARAM = MulticoreParam(
      workers = as.integer(Sys.getenv("NSLOTS")),
      progressbar = TRUE),
    verbose = TRUE)
  saveHDF5SummarizedExperiment(
    x = bsseq,
    dir = file.path(tempdir(), seqlevel),
    verbose = TRUE)
})

saveHDF5SummarizedExperiment(
  x = do.call(rbind, list_of_bsseq),
  dir = paste0(scratch_dir, ".small_smooth"),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCT (+) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase1_brain_samples.mCT_pos"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

ct_pos <- findLoci(
  pattern = "CT",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "+")
seqinfo(ct_pos) <- seqinfo

files <- mapply(
  function(sample_id, hiseq_run) {
    list.files(
      path = here(
        "HDF5",
        "extdata",
        "Bulk_GTEx_Brain_hg38",
        hiseq_run,
        "nonCG"),
      pattern = glob2rx(paste0(sample_id, "*.CX_report.txt.gz")),
      full.names = TRUE)
  },
  sample_id = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  hiseq_run = CD[["HiSeq_Run"]],
  SIMPLIFY = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = ct_pos,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

dir.create(file.path(scratch_dir, "collapseBSseq"))
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

file.copy(
  from = file.path(scratch_dir, "collapseBSseq"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCT (+) small smooth ---------------------------------------------------------

# NOTE: For unknown reasons it was crashing when trying to smooth the entire
#       object (even though only a subset of data is loaded into memory). So I
#       had to smooth on a per-chromosome and then rbind() the outputs.
# NOTE: Run with 8 cores, 50 Gb each core.

# TODO: Running on compute-051

dir.create(paste0(scratch_dir, ".small_smooth.per_chromosome"))
file.copy(
  from = list.files(
    path = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq", out_name),
    full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth.per_chromosome"))
bsseq <- loadHDF5SummarizedExperiment(
  dir = paste0(scratch_dir, ".small_smooth.per_chromosome"))

list_of_bsseq <- lapply(seqlevels(bsseq), function(seqlevel) {
  message(seqlevel,
          paste0(
            " (",
            match(seqlevel, seqlevels(bsseq)),
            "/",
            length(seqlevels(bsseq)),
            ")"))
  bsseq <- BSmooth(
    BSseq = keepSeqlevels(x = bsseq, value = seqlevel, pruning.mode = "coarse"),
    ns = 200,
    h = 3000,
    maxGap = 10^8,
    keep.se = FALSE,
    BPPARAM = MulticoreParam(
      workers = as.integer(Sys.getenv("NSLOTS")),
      progressbar = TRUE),
    verbose = TRUE)
  saveHDF5SummarizedExperiment(
    x = bsseq,
    dir = file.path(tempdir(), seqlevel),
    verbose = TRUE)
})

saveHDF5SummarizedExperiment(
  x = do.call(rbind, list_of_bsseq),
  dir = paste0(scratch_dir, ".small_smooth"),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCT (-) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase1_brain_samples.mCT_neg"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

ct_neg <- findLoci(
  pattern = "CT",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "-")
seqinfo(ct_neg) <- seqinfo

files <- mapply(
  function(sample_id, hiseq_run) {
    list.files(
      path = here(
        "HDF5",
        "extdata",
        "Bulk_GTEx_Brain_hg38",
        hiseq_run,
        "nonCG"),
      pattern = glob2rx(paste0(sample_id, "*.CX_report.txt.gz")),
      full.names = TRUE)
  },
  sample_id = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  hiseq_run = CD[["HiSeq_Run"]],
  SIMPLIFY = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = ct_neg,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

dir.create(file.path(scratch_dir, "collapseBSseq"))
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

file.copy(
  from = file.path(scratch_dir, "collapseBSseq"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCT (-) small smooth ---------------------------------------------------------

# NOTE: For unknown reasons it was crashing when trying to smooth the entire
#       object (even though only a subset of data is loaded into memory). So I
#       had to smooth on a per-chromosome and then rbind() the outputs.
# NOTE: Run with 8 cores, 50 Gb each core.

# TODO: Not yet running.

dir.create(paste0(scratch_dir, ".small_smooth.per_chromosome"))
file.copy(
  from = list.files(
    path = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq", out_name),
    full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth.per_chromosome"))
bsseq <- loadHDF5SummarizedExperiment(
  dir = paste0(scratch_dir, ".small_smooth.per_chromosome"))

list_of_bsseq <- lapply(seqlevels(bsseq), function(seqlevel) {
  message(seqlevel,
          paste0(
            " (",
            match(seqlevel, seqlevels(bsseq)),
            "/",
            length(seqlevels(bsseq)),
            ")"))
  bsseq <- BSmooth(
    BSseq = keepSeqlevels(x = bsseq, value = seqlevel, pruning.mode = "coarse"),
    ns = 200,
    h = 3000,
    maxGap = 10^8,
    keep.se = FALSE,
    BPPARAM = MulticoreParam(
      workers = as.integer(Sys.getenv("NSLOTS")),
      progressbar = TRUE),
    verbose = TRUE)
  saveHDF5SummarizedExperiment(
    x = bsseq,
    dir = file.path(tempdir(), seqlevel),
    verbose = TRUE)
})

saveHDF5SummarizedExperiment(
  x = do.call(rbind, list_of_bsseq),
  dir = paste0(scratch_dir, ".small_smooth"),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCC (+) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase1_brain_samples.mCC_pos"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

cc_pos <- findLoci(
  pattern = "CC",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "+")
seqinfo(cc_pos) <- seqinfo

files <- mapply(
  function(sample_id, hiseq_run) {
    list.files(
      path = here(
        "HDF5",
        "extdata",
        "Bulk_GTEx_Brain_hg38",
        hiseq_run,
        "nonCG"),
      pattern = glob2rx(paste0(sample_id, "*.CX_report.txt.gz")),
      full.names = TRUE)
  },
  sample_id = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  hiseq_run = CD[["HiSeq_Run"]],
  SIMPLIFY = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = cc_pos,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

dir.create(file.path(scratch_dir, "collapseBSseq"))
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

file.copy(
  from = file.path(scratch_dir, "collapseBSseq"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCC (+) small smooth ---------------------------------------------------------

# NOTE: For unknown reasons it was crashing when trying to smooth the entire
#       object (even though only a subset of data is loaded into memory). So I
#       had to smooth on a per-chromosome and then rbind() the outputs.
# NOTE: Run with 8 cores, 50 Gb each core.

# TODO: Not yet running.

dir.create(paste0(scratch_dir, ".small_smooth.per_chromosome"))
file.copy(
  from = list.files(
    path = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq", out_name),
    full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth.per_chromosome"))
bsseq <- loadHDF5SummarizedExperiment(
  dir = paste0(scratch_dir, ".small_smooth.per_chromosome"))

list_of_bsseq <- lapply(seqlevels(bsseq), function(seqlevel) {
  message(seqlevel,
          paste0(
            " (",
            match(seqlevel, seqlevels(bsseq)),
            "/",
            length(seqlevels(bsseq)),
            ")"))
  bsseq <- BSmooth(
    BSseq = keepSeqlevels(x = bsseq, value = seqlevel, pruning.mode = "coarse"),
    ns = 200,
    h = 3000,
    maxGap = 10^8,
    keep.se = FALSE,
    BPPARAM = MulticoreParam(
      workers = as.integer(Sys.getenv("NSLOTS")),
      progressbar = TRUE),
    verbose = TRUE)
  saveHDF5SummarizedExperiment(
    x = bsseq,
    dir = file.path(tempdir(), seqlevel),
    verbose = TRUE)
})

saveHDF5SummarizedExperiment(
  x = do.call(rbind, list_of_bsseq),
  dir = paste0(scratch_dir, ".small_smooth"),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCC (-) ----------------------------------------------------------------------

out_name <- "eGTEx.Phase1_brain_samples.mCC_neg"
scratch_dir <- file.path(tempdir(), out_name)
stopifnot(!dir.exists(scratch_dir))

cc_neg <- findLoci(
  pattern = "CC",
  subject = fasta,
  include = seqlevels(fasta),
  strand = "-")
seqinfo(cc_neg) <- seqinfo

files <- mapply(
  function(sample_id, hiseq_run) {
    list.files(
      path = here(
        "HDF5",
        "extdata",
        "Bulk_GTEx_Brain_hg38",
        hiseq_run,
        "nonCG"),
      pattern = glob2rx(paste0(sample_id, "*.CX_report.txt.gz")),
      full.names = TRUE)
  },
  sample_id = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  hiseq_run = CD[["HiSeq_Run"]],
  SIMPLIFY = TRUE)

bsseq <- read.bismark(
  files = files,
  loci = cc_neg,
  colData = CD,
  rmZeroCov = FALSE,
  # NOTE: Not collapsing loci by strand.
  strandCollapse = FALSE,
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  BACKEND = "HDF5Array",
  dir = scratch_dir,
  verbose = TRUE)

dir.create(file.path(scratch_dir, "collapseBSseq"))
bsseq <- collapseBSseq(
  BSseq = bsseq,
  group = CD[["Sample_ID_for_Data_Sharing_and_Public_Release"]],
  BPPARAM = MulticoreParam(workers = 3, progressbar = TRUE),
  dir = file.path(scratch_dir, "collapseBSseq"),
  chunkdim = NULL,
  level = NULL,
  type = "integer")

file.copy(
  from = file.path(scratch_dir, "collapseBSseq"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)

# mCC (-) small smooth ---------------------------------------------------------

# NOTE: For unknown reasons it was crashing when trying to smooth the entire
#       object (even though only a subset of data is loaded into memory). So I
#       had to smooth on a per-chromosome and then rbind() the outputs.
# NOTE: Run with 8 cores, 50 Gb each core.

# TODO: Not yet running.

dir.create(paste0(scratch_dir, ".small_smooth.per_chromosome"))
file.copy(
  from = list.files(
    path = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq", out_name),
    full.names = TRUE),
  to = paste0(scratch_dir, ".small_smooth.per_chromosome"))
bsseq <- loadHDF5SummarizedExperiment(
  dir = paste0(scratch_dir, ".small_smooth.per_chromosome"))

list_of_bsseq <- lapply(seqlevels(bsseq), function(seqlevel) {
  message(seqlevel,
          paste0(
            " (",
            match(seqlevel, seqlevels(bsseq)),
            "/",
            length(seqlevels(bsseq)),
            ")"))
  bsseq <- BSmooth(
    BSseq = keepSeqlevels(x = bsseq, value = seqlevel, pruning.mode = "coarse"),
    ns = 200,
    h = 3000,
    maxGap = 10^8,
    keep.se = FALSE,
    BPPARAM = MulticoreParam(
      workers = as.integer(Sys.getenv("NSLOTS")),
      progressbar = TRUE),
    verbose = TRUE)
  saveHDF5SummarizedExperiment(
    x = bsseq,
    dir = file.path(tempdir(), seqlevel),
    verbose = TRUE)
})

saveHDF5SummarizedExperiment(
  x = do.call(rbind, list_of_bsseq),
  dir = paste0(scratch_dir, ".small_smooth"),
  verbose = TRUE)
file.copy(
  from = paste0(scratch_dir, ".small_smooth"),
  to = here("HDF5", "extdata", "Bulk_GTEx_Brain_hg38", "BSseq"),
  recursive = TRUE)
