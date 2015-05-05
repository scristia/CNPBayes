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

#' @export
cnProbability <- function(prob, K){
  pz2 <- prob[, 1]
  for(i in seq_along(2:K)){
    nonzero.prob <- prob[, i + 1] > 0
    pz2[ nonzero.prob ] <- (prob[, i + 1] + i)[ nonzero.prob ]
  }
  pz2
}

#' @export
imputeFromSampledData <-  function(model, data, index){
  if(is.null(names(data))) stop("data must be a named vector")
  pz <- probz(model)
  cn <- map(model)
  pz2 <- cnProbability(pz, k(model))
  ##r2 <- quantile(data, probs=seq(0, 1, by=0.01))
  df <- data.frame(r=data,
                   cn=rep(NA, length(data)),
                   p=rep(NA, length(data)),
                   row.names=names(data))
  df$p[index] <- pz2
  df$cn[index] <- cn
  tmp <- ntile(df$r, 1000)
  df$quantiles <- ntile(df$r, 1000)

  missing_quantiles <- df$quantiles[is.na(df$cn)]
  complete_quantiles <- df$quantiles[!is.na(df$cn)]
  missing_quantiles <- missing_quantiles[!missing_quantiles %in% complete_quantiles]
  names(missing_quantiles) <- missing_quantiles
  missing_quantiles2 <- unique(missing_quantiles)
  nearest_quantiles <- sapply(missing_quantiles2, function(x, complete){
    complete[which.min(abs(complete-x))]
  }, complete=complete_quantiles)
  names(nearest_quantiles) <- missing_quantiles2
  ## same length as original data.frame
  nearest_quantiles2 <- nearest_quantiles[as.character(df$quantiles)]
  nearest_quantiles2[is.na(nearest_quantiles2)] <- df$quantiles[ is.na(nearest_quantiles2) ]
  df$quantiles2 <- nearest_quantiles2
  ##df[df$quantiles != df$quantiles2, ]

  df.complete <- df[!is.na(df$cn), ]
  df.incomplete <- df[is.na(df$cn), ]
  ## now any missing cn/prob can be 'imputed'
  cn <- setNames(df.complete$cn, df.complete$quantiles2)
  df.incomplete$cn <- cn[as.character(df.incomplete$quantiles2)]
  p <- setNames(df.complete$p, df.complete$quantiles2)
  df.incomplete$p <- p[as.character(df.incomplete$quantiles2)]
  df2 <- rbind(df.incomplete, df.complete)
  df2 <- df2[rownames(df), ]
  df2
}
