# Run LDSC with eGTEx categories (made in make_LDSC_annot.R)
# Peter Hickey
# 2020-05-21

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
message("category = ", names(categories)[i])

gwasss <- list.files(
  path = here("SLDSR", "extdata", "ldsc_GRCh38", "munge_sumstats", "Phase1"),
  full.names = TRUE,
  pattern = glob2rx("*.sumstats.gz"))

# Run LDSC (adjusting for baseline) --------------------------------------------

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    cmd <- paste0(
      "python /users/phickey/software/ldsc_hg38/ldsc.py ",
      "--h2 ", x, " ",
      paste(
        "--w-ld-chr ",
        here("SLDSR", "extdata", "ldsc_GRCh38", "weights", "weights.hm3_noMHC. ")),
      paste(
        "--ref-ld-chr ",
        here("SLDSR", "output", "ldsc/"), cn, ".,",
        here("SLDSR", "extdata", "ldsc_GRCh38", "baselineLD_v2.2", "baselineLD. "),
        sep = ""),
      "--overlap-annot ",
      paste(
        "--frqfile-chr",
        here("SLDSR", "extdata", "ldsc_GRCh38", "plink_files", "1000G.EUR.hg38. ")),
      paste(
        "--out ",
        here("SLDSR", "output", "ldsc/"), cn, ".", bn, " ",
        sep = ""),
      "--print-coefficients")
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})
