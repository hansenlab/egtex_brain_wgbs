# Make .annot file for use with ldsc containing eGTEx-derived genomic features
# with hg38 co-ordinates.
# Peter Hickey
# 2020-04-20

# NOTE: Run with module load conda_R/3.6.x

# Setup ------------------------------------------------------------------------

library(readr)
library(here)
library(GenomicRanges)
library(janitor)

# Construct eGTEx features -----------------------------------------------------

dmr_files <- list.files(here("DMRs", "objects"), full.names = TRUE)
names(dmr_files) <- gsub("\\.rds", "", basename(dmr_files))
# NOTE: Not running LDSC on sex DMRs of 'Phase1' DMRs.
dmr_files <- dmr_files[!grepl("^sex|Phase1", names(dmr_files))]
list_of_dmrs_gr <- lapply(names(dmr_files), function(n) {
  dmr_file <- dmr_files[[n]]
  gr <- as(readRDS(dmr_file), "GRanges")
  if (grepl("pos-", n)) {
    strand(gr) <- "+"
  } else if (grepl("neg-", n)) {
    strand(gr) <- "-"
  } else {
    strand(gr) <- "*"
  }
  genome(gr) <- "hg38"
  gr
})
names(list_of_dmrs_gr) <- names(dmr_files)

# Add unstranded CH-DMRs.
groups <- c("general", "basal_ganglia", "hippocampus", "overall")
names(groups) <- paste0(groups, "_CH-DMRs")
list_of_dmrs_gr <- c(
  list_of_dmrs_gr,
  lapply(groups, function(g) {
    idx <- grepl(g, names(list_of_dmrs_gr)) &
      !grepl("CG", names(list_of_dmrs_gr))
    gr <- list_of_dmrs_gr[idx]
    gr <- reduce(Reduce(c, gr), ignore.strand = TRUE)
    genome(gr) <- "hg38"
    gr
  }))

# 'Old' VMRs (aka VMRs-095)
load(here("VMR_analyses", "objects", "list_of_autosomal_vmrs.nq.rda"))
list_of_autosomal_vmrs.nq <- lapply(list_of_autosomal_vmrs.nq, function(x) {
  x <- as(x, "GRanges")
  genome(x) <- "hg38"
  x
})
# NOTE: Focus on brain VMRs.
list_of_autosomal_vmrs.nq <-
  list_of_autosomal_vmrs.nq[grep("Brain|HC", names(list_of_autosomal_vmrs.nq))]
# Kludge to clean up the names.
names(list_of_autosomal_vmrs.nq) <- paste0(
  gsub(
    "Brain_",
    "",
    colnames(
      clean_names(
        as.data.frame(
          lapply(list_of_autosomal_vmrs.nq, function(x) vector()),
          check.names = FALSE),
        case = "parsed"))),
  "_CG-VMRs")

# 'New' VMRs (aka VMRs-064)
load(here("VMR_analyses", "objects", "NEW_VMR_lists_for_Pete.rda"))
list_of_new_autosomal_vmrs <- lapply(
  list_of_filtered_autosomal_vmrs_wSDS_gr,
  function(x) {
    genome(x) <- "hg38"
    x
  })
# NOTE: Focus on brain VMRs.
list_of_new_autosomal_vmrs <-
  list_of_new_autosomal_vmrs[grep("Brain|HC", names(list_of_new_autosomal_vmrs))]
# Kludge to clean up the names.
names(list_of_new_autosomal_vmrs) <- paste0(
  "New_",
  gsub(
    "Brain_",
    "",
    colnames(
      clean_names(
        as.data.frame(
          lapply(list_of_new_autosomal_vmrs, function(x) vector()),
          check.names = FALSE),
        case = "parsed"))),
  "_CG-VMRs")

# 'Old' VMRs excluding the general CG-DMRs.
list_of_autosomal_vmrs.nq_subsetByOverlaps_dmrs <- lapply(
  list_of_autosomal_vmrs.nq,
  function(vmrs) {
    subsetByOverlaps(vmrs, list_of_dmrs_gr[["general_CG-DMRs"]], invert = TRUE)
  })
names(list_of_autosomal_vmrs.nq_subsetByOverlaps_dmrs) <- paste0(
  names(list_of_autosomal_vmrs.nq_subsetByOverlaps_dmrs),
  "_subsetByOverlaps_DMRs")

# 'New' VMRs excluding the general CG-DMRs.
list_of_new_autosomal_vmrs_subsetByOverlaps_dmrs <- lapply(
  list_of_new_autosomal_vmrs,
  function(vmrs) {
    subsetByOverlaps(vmrs, list_of_dmrs_gr[["general_CG-DMRs"]], invert = TRUE)
  })
names(list_of_new_autosomal_vmrs_subsetByOverlaps_dmrs) <- paste0(
  names(list_of_new_autosomal_vmrs_subsetByOverlaps_dmrs),
  "_subsetByOverlaps_DMRs")

# NOTE: No need to liftOver any of these because they are already on hg38.
categories <- c(
  list_of_dmrs_gr,
  list_of_autosomal_vmrs.nq,
  list_of_new_autosomal_vmrs,
  list_of_autosomal_vmrs.nq_subsetByOverlaps_dmrs,
  list_of_new_autosomal_vmrs_subsetByOverlaps_dmrs)

# NOTE: This reduces the size of the data by only keeping the necessary parts
#       of the objects and by sorting.
categories <- lapply(categories, function(x) sort(granges(x)))

# Make annotations -------------------------------------------------------------

saveRDS(
  categories,
  here("SLDSR", "objects", "eGTEx_features.rds"), compress = "xz")

mclapply(1:22, function(sl) {
  message(sl)
  cds <- read_tsv(
    here(
      "SLDSR",
      "extdata",
      "ldsc_GRCh38",
      "baselineLD_v2.2/",
      paste0("baselineLD.", sl, ".annot.gz")))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  genome(cds_gr) <- "hg38"
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  annot[names(categories)] <- mclapply(names(categories), function(cn) {
    stopifnot(isDisjoint(categories[[cn]]))
    as.integer(overlapsAny(cds_gr, categories[[cn]]))
  }, mc.cores = 4)
  # 'Marginal' annotation file
  mclapply(names(categories), function(cn) {
    fl <- here(
      "SLDSR",
      "output",
      "ldsc",
      paste0(cn, ".", sl, ".annot.gz"))
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", cn)], fl)
  }, mc.cores = 4)
}, mc.cores = 10)
