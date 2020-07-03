# Make CpG island (CGI)related features for hg19
# Peter Hickey
# 2016-06-01
# Edited by Lindsay Rizzardi
# 2018-08-30

library(rtracklayer)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download UCSC CGI definitions for hg38
###

session <- browserSession("UCSC")
genome(session) <- "hg38"
cgis <- track(ucscTableQuery(session, track = "cpgIslandExt"))
# NOTE: Only retaining standard seqlevels
cgis <- keepStandardChromosomes(cgis,pruning.mode="coarse")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Construct CGI shores and shelves
###

# CGI shores defined as CGI +/- 2000bp
shores <- sort(c(flank(cgis, width = 2000, start = TRUE),
            flank(cgis, width = 2000, start = FALSE)))
stopifnot(length(shores) == (2 * length(cgis)))
stopifnot(all(width(shores) == 2000))

# CGI shelves are defined as as +/- 2000-4000bp of CGI
shelves <- sort(c(flank(flank(cgis, width = 2000, start = TRUE),
                        width = 2000, start = TRUE),
                  flank(flank(cgis, width = 2000, start = FALSE),
                        width = 2000, start = FALSE)))
stopifnot(length(shelves) == (2 * length(cgis)))
stopifnot(all(width(shelves) == 2000))

##### make open sea

x <- c(cgis, shores, shelves)
y <- gaps(x)
open_sea  <- y[strand(y) == "*"]

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Save objects
###

save(cgis, shores, shelves, open_sea,
     file = "/dcl01/FB2/data/personal/gtex/gtex/eGTEx_analysis/objects/CGI-related-features-hg38.rda")
