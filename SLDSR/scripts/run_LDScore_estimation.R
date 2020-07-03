# Run LD score estimation with eGTEx categories (made in make_LDSC_annot.R).
# Peter Hickey
# 2020-04-20

# Setup ------------------------------------------------------------------------

# NOTE: On JHPCE, requires `module load python/2.7.9`

args <- commandArgs(TRUE)
i <- as.integer(args[1])
message("i = ", i)

library(parallel)
library(here)

options("mc.cores" = 10)

# Load data --------------------------------------------------------------------

categories <- readRDS(here("SLDSR", "objects", "eGTEx_features.rds"))

seqlevels <- 1:22

message("category = ", names(categories)[i])

# Run LD Score estimation ------------------------------------------------------

lapply(names(categories)[i], function(cn) {
  message(cn)
  mclapply(seqlevels, function(sl) {
    cmd <- paste0(
      "python /users/phickey/software/ldsc_hg38/ldsc.py ",
      "--l2 ",
      paste(
        "--bfile ",
        here("SLDSR", "extdata", "ldsc_GRCh38", "plink_files", "1000G.EUR.hg38."),
        sl,
        " ",
        sep = ""),
      "--ld-wind-cm 1 ",
      paste(
        "--annot ",
        here("SLDSR", "output", "ldsc", paste0(cn, ".", sl, ".annot.gz "))),
      paste(
        "--out ",
        here("SLDSR", "output", "ldsc", paste0(cn, ".", sl, " "))),
      paste(
        "--print-snps ",
        here("SLDSR", "extdata", "ldsc_GRCh38", "list.txt")))
      print(cmd)
      system(cmd)
  }, mc.cores = 10)
})
