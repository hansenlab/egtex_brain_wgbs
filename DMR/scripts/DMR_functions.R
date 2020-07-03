# Functions used in DMR calling.
# Peter Hickey
# 2018-11-12

# NOTE: These are Modified versions of bsseq functions to support chunking.
#       Brief experimentation shows that lmFit() takes 5-10x as much memory as
#       Y occupies (smaller if there are more columns).
#       E.g., 100 Mb chunks of data means 580 Mb (128 columns) or 890 Mb (8
#       columns), both using 5 groups in a one-way layout.
#       The running time of lmFit() scales linearly in the number of rows of Y.
BSmooth.fstat <- function(BSseq, design, contrasts, chunksize = NULL,
                          verbose = TRUE) {
  stopifnot(is(BSseq, "BSseq"))
  stopifnot(hasBeenSmoothed(BSseq))

  ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
  ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")

  if(verbose) cat("[BSmooth.fstat] fitting linear models ... ")
  ptime1 <- proc.time()
  if (is.null(chunksize)) {
    chunksize <- nrow(BSseq)
  }
  allPs <- getMeth(BSseq, type = "smooth", withDimnames = FALSE)
  chunks <- IRanges(start = seq(1, nrow(BSseq), chunksize), width = chunksize)
  end(chunks) <- pmin(end(chunks), nrow(BSseq))
  stats <- lapply(seq_along(chunks), function(i) {
    chunk <- chunks[i]
    y <- as.matrix(allPs[seq(start(chunk), end(chunk)), , drop = FALSE])
    fit <- limma::lmFit(y, design)
    fitC <- limma::contrasts.fit(fit, contrasts)
    ## Need
    ##   fitC$coefficients, fitC$stdev.unscaled, fitC$sigma, fitC$cov.coefficients
    ## actuall just need
    ##   tstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    ##   rawSds <- fitC$sigma
    ##   cor.coefficients <- cov2cor(fitC$cov.coefficients)
    rawSds <- as.matrix(fitC$sigma)
    cor.coefficients <- cov2cor(fitC$cov.coefficients)
    rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
    names(dimnames(rawTstats)) <- NULL

    list(
      rawSds = rawSds,
      cor.coefficients = cor.coefficients,
      rawTstats = rawTstats)
  })
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if(verbose) cat(sprintf("done in %.1f sec\n", stime))

  parameters <- c(
    BSseq@parameters,
    list(design = design, contrasts = contrasts))
  stats <- list(
    rawSds = do.call(rbind, lapply(stats, "[[", "rawSds")),
    # All cor.coefficients are the same.
    cor.coefficients = stats[[1]][["cor.coefficients"]],
    rawTstats = do.call(rbind, lapply(stats, "[[", "rawTstats")))
  out <- BSseqStat(
    gr = granges(BSseq),
    stats = stats,
    parameters = parameters)
  out
}

getNullDistribution_BSmooth.fstat <- function(BSseq, idxMatrix, design,
                                              contrasts, coef = NULL, cutoff,
                                              maxGap.sd, maxGap.dmr,
                                              chunksize = NULL, mc.cores = 1) {
  message(
    sprintf(
      "[getNullDistribution_BSmooth.fstat] performing %d permutations\n",
      nrow(idxMatrix)))
  nullDist <- mclapply(seq_len(nrow(idxMatrix)), function(ii) {
    ptime1 <- proc.time()
    bstat <- BSmooth.fstat(
      BSseq = BSseq,
      design = design[idxMatrix[ii, ], ],
      contrasts = contrasts,
      chunksize = chunksize)
    bstat <- bsseq:::smoothSds(bstat, maxGap = maxGap.sd)
    bstat <- bsseq:::computeStat(bstat, coef = coef)
    dmrs0 <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    message(
      sprintf(
        "[getNullDistribution_BSmooth.fstat] completing permutation %d in %.1f sec\n",
        ii,
        stime))
    dmrs0
  }, mc.cores = min(nrow(idxMatrix), mc.cores), mc.preschedule = TRUE)
  nullDist
}

fstat.pipeline <- function(BSseq, design, contrasts, cutoff, fac, nperm = 1000,
                           coef = NULL, maxGap.sd = 10^8, maxGap.dmr = 300,
                           type = "dmrs", chunksize = NULL, mc.cores = 1) {
  type <- match.arg(type, c("dmrs", "blocks"))
  bstat <- BSmooth.fstat(
    BSseq = BSseq,
    design = design,
    contrasts = contrasts,
    chunksize = chunksize)
  bstat <- bsseq:::smoothSds(bstat)
  bstat <- bsseq:::computeStat(bstat, coef = coef)
  dmrs <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
  if (is.null(dmrs)) {
    stop("No DMRs identified. Consider reducing the 'cutoff' from (",
         paste0(cutoff, collapse = ", "), ")")
  }
  idxMatrix <- bsseq:::permuteAll(nperm, design)
  nullDist <- getNullDistribution_BSmooth.fstat(
    BSseq = BSseq,
    idxMatrix = idxMatrix,
    design = design,
    contrasts = contrasts,
    coef = coef,
    cutoff = cutoff,
    maxGap.sd = maxGap.sd,
    maxGap.dmr = maxGap.dmr,
    chunksize = chunksize,
    mc.cores = mc.cores)
  fwer <- bsseq:::getFWER.fstat(null = c(list(dmrs), nullDist), type = type)
  dmrs$fwer <- fwer
  meth <- getMeth(BSseq, dmrs, what = "perRegion")
  meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
  dmrs <- cbind(dmrs, meth)
  dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)
  list(bstat = bstat, dmrs = dmrs, idxMatrix = idxMatrix, nullDist = nullDist)
}
