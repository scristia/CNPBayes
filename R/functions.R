posteriorPrecisionConjugateNormal <- function(prior.precision, data.precision) {
  prior.precision+data.precision
}

posteriorMeanConjugateNormal <- function(prior.precision, data.precision,
                                         posterior.precision, prior.mean, data.mean){
  (prior.precision/posterior.precision)*prior.mean + (data.precision/posterior.precision) * data.mean
}

stopif <- function(x) stopifnot(!x)

precision <- function(x) 1/var(x, na.rm=TRUE)

gammaShapeRate <- function(mn, sd){
  ##shape/rate   = mn (1)
  ##shape/rate^2 = sd (2)
  ##shape = mn * rate
  ##mn * rate / rate^2 = sd
  ##mn/rate = sd
  ##mn/sd = rate
  rate <- mn/sd
  shape <- mn*rate
  setNames(c(shape, rate), c("shape", "rate"))
}


#' @export
consensusRegion <- function(g){
  ## Defined the consensus start as the minimum basepair that is
  ## spanned by at least half of all identified CNVs
  ##cat(".")
  grl <- split(g, g$id)
  if(length(grl)/length(g) < 0.9){
    ## more than 10% of individuals have more than one cnv in the region
    ##
    ## assume that these individuals should have one cnv
    grl1 <- grl[elementLengths(grl)==1]
    g1 <- unlist(grl1)
    if(class(g1)=="list") browser()
    g1 <- GRanges(as.character(seqnames(g1)), IRanges(start(g1), end(g1)))
    grl2ormore <- grl[elementLengths(grl) >= 2]
    grl3 <- foreach(g=grl2ormore) %do% reduce(g, min.gapwidth=1e6)
    g2ormore <- unlist(GRangesList(grl3))
    g <- c(g1, g2ormore)
  }
  threshold <- floor(length(g)/2)
  dj <- disjoin(g)
  cnt <- countOverlaps(disjoin(g), g)
  if(all(cnt < threshold)){
    threshold <- floor(threshold/2)
  }
  ## the disjoint intervals are sorted, so we only need to take the first interval that passes this threshold
  index.start <- min(which(cnt > threshold))
  index.end <- max(which(cnt > threshold))
  startg <- dj[index.start]
  endg <- dj[index.end]
  if(startg != endg){
    region <- GRanges(as.character(seqnames(startg)), IRanges(start(startg), end(endg)))
  } else region <- startg
  region
}


#' @export
defineCnpRegions <- function(grl, thr=0.02){
  message("unlist GRangesList...")
  g <- GenomicRanges::unlist(grl)
  names(g) <- NULL
  dj <- disjoin(g)
  ##g2 <- GRanges(seqnames(g), IRanges(start(g), end(g)))
  ##dj <- disjoin(g2)
  prop.cnv <- countOverlaps(dj, g)/length(grl)

  is_cnp <- prop.cnv > thr

  ## throw-out large CNVs and CNVs that do not overlap any of the CNPs
  gsmall <- g[width(g) < 200e3]
  gsmall <- subsetByOverlaps(gsmall, dj[is_cnp])
  regions <- reduce(gsmall) ## 267 regions

  ##library(GenomicRanges)
  ##load("~/MyPapers/aricUricAcid/data/chr4locus.rda")
  ## verify chr4 locus is present
  ##stopifnot(identical(length(subjectHits(findOverlaps(regions, chr4locus))),2L))

  ## Split the cnvs by region
  hits <- findOverlaps(regions, g)
  hitlist <- split(subjectHits(hits), queryHits(hits))
  regionlist <- vector("list", length(hitlist))
  message("find consensus regions...")
  ##browser()
  for(j in seq_along(hitlist)){
    cat(".")
    i <- hitlist[[j]]
    regionlist[[j]] <- consensusRegion(g[i])
  }
  ##regionlist <- GRangesList(foreach(i=hitlist) %do% consensusRegion(g[i]))
  regions <- GenomicRanges::GRangesList(regionlist)
  GenomicRanges::unlist(regions)
}

#' @export
consensusCNP <- function(grl, transcripts, min.width=2e3, min.prevalance=0.02){
  g <- as(unlist(grl), "GRanges")
  names(g) <- NULL
  grl <- split(g, g$id)
  grl <- grl[colnames(views)]
  regions <- defineCnpRegions(grl, thr=min.prevalance)
  regions <- regions[ width(regions) > min.width ]
  regions <- reduce(regions)
  regions <- annotateRegions(regions, transcripts)
  regions
}

annotateRegions <- function(regions, transcripts){
  hits <- findOverlaps(regions, transcripts)
  indices <- split(subjectHits(hits), queryHits(hits))
  hgnc.ids <- sapply(indices, function(i, tx){
    paste(unique(tx$hgnc[i]), collapse=",")
  }, tx=transcripts)
  drivers <- sapply(indices, function(i, tx){
    hgnc <- tx$hgnc[i]
    is.cancer.connection <- tx$cancer_connection[i]
    paste(unique(hgnc[is.cancer.connection]), collapse=",")
  }, tx=transcripts)
  regions$driver <- regions$hgnc <- ""
  i <- as.integer(names(indices))
  regions$driver[i] <- drivers
  regions$hgnc[i] <- hgnc.ids
  regions
}


colLRRMediansSE <- function(object){
  R <- lrr(object)
  colMedians(R, na.rm=TRUE)
}

computeLRRMedians <- function(views, regions){
  hits <- findOverlaps(regions, rowRanges(views))
  R <- lrr(views[subjectHits(hits), ])
  viewsCNP <- views[subjectHits(hits), ]
  indices <- split(seq_len(nrow(viewsCNP)), queryHits(hits))
  lrr.meds <- foreach(i = indices, .combine="rbind", .packages="matrixStats") %do% {
    mat <- R[i, , drop=FALSE]
    colMedians(mat, na.rm=TRUE)
  }
  dimnames(lrr.meds) <- list(names(regions), colnames(views))
  iR <- integerMatrix(lrr.meds, 1000)
  se <- SummarizedExperiment(assays=SimpleList(medr=iR),
                             rowData=regions,
                             colData=colData(views))
  se
}
