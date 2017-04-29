consensusRegion <- function(g){
  ## Defined the consensus start as the minimum basepair that is
  ## spanned by at least half of all identified CNVs
  ##cat(".")
  grl <- split(g, g$id)
  if(length(grl)/length(g) < 0.9){
    ## more than 10% of individuals have more than one cnv in the region
    ##
    ## assume that these individuals should have one cnv
    grl1 <- grl[elementNROWS(grl)==1]
    g1 <- unlist(grl1)
    #if(class(g1)=="list") browser()
    g1 <- GRanges(as.character(seqnames(g1)), IRanges(start(g1), end(g1)))
    grl2ormore <- grl[elementLengths(grl) >= 2]
    grl3 <- unname(lapply(grl2ormore, reduce, min.gapwidth=1e6))
    g2ormore <- unlist(GRangesList(grl3))
    g <- c(g1, g2ormore)
  }
  threshold <- floor(length(g)/2)
  dj <- disjoin(g)
  cnt <- countOverlaps(disjoin(g), g)
  while(all(cnt <= threshold)){
    threshold <- floor(threshold/2)
  }
  ## the disjoint intervals are sorted, so we only need to take the
  ## first interval that passes this threshold
  index.start <- tryCatch(min(which(cnt > threshold)), error=function(e) NULL)
  index.end <- max(which(cnt > threshold))
  startg <- dj[index.start]
  endg <- dj[index.end]
  if(startg != endg){
    region <- GRanges(as.character(seqnames(startg)), IRanges(start(startg), end(endg)))
  } else region <- startg
  region
}

defineCnpRegions <- function(grl, thr=0.02){
  message("unlist GRangesList...")
  g <- unlist(grl)
  names(g) <- NULL
  dj <- disjoin(g)
  prop.cnv <- countOverlaps(dj, g)/length(grl)
  is_cnp <- prop.cnv > thr

  ## throw-out large CNVs and CNVs that do not overlap any of the CNPs
  ## gsmall <- g[width(g) < 200e3]
  gsmall <- subsetByOverlaps(g, dj[is_cnp])
  regions <- reduce(gsmall) ## 267 regions

  ## Split the cnvs by region
  hits <- findOverlaps(regions, g)
  hitlist <- split(subjectHits(hits), queryHits(hits))
  regionlist <- vector("list", length(hitlist))
  message("find consensus regions...")
  for(j in seq_along(hitlist)){
    message(".", appendLF=FALSE)
    i <- hitlist[[j]]
    regionlist[[j]] <- consensusRegion(g[i])
  }
  message("")
  ##regionlist <- GRangesList(foreach(i=hitlist) %do% consensusRegion(g[i]))
  regions <- GRangesList(regionlist)
  unlist(regions)
}

#' Identify consensus start and stop coordinates of a copy number
#' polymorphism
#'
#' The collection of copy number variants (CNVs) identified in a study
#' can be encapulated in a GRangesList, where each element is a
#' GRanges of the CNVs identified for an individual.  (For a study
#' with 1000 subjects, the GRangesList object would have length 1000
#' if each individual had 1 or more CNVs.)  For regions in which CNVs
#' occur in more than 2 percent of study participants, the start and
#' end boundaries of the CNVs may differ because of biological
#' differences in the CNV size as well as due to technical noise of
#' the assay and the uncertainty of the breakpoints identified by a
#' segmentation of the genomic data.  Among subjects with a CNV called
#' at a given locus, the \code{consensusCNP} function identifies the
#' largest region that is copy number variant in half of these
#' subjects.
#'
#' @examples
#' library(GenomicRanges)
#' ##
#' ## Simulate 2 loci at which CNVs are common
#' ##
#' set.seed(100)
#' starts <- rpois(1000, 100) + 10e6L
#' ends <- rpois(1000, 100) + 10.1e6L
#' cnv1 <- GRanges("chr1", IRanges(starts, ends))
#' cnv1$id <- paste0("sample", seq_along(cnv1))
#'
#' starts <- rpois(500, 1000) + 101e6L
#' ends <- rpois(500, 1000) + 101.4e6L
#' cnv2 <- GRanges("chr5", IRanges(starts, ends))
#' cnv2$id <- paste0("sample", seq_along(cnv2))
#'
#' ##
#' ## Simulate a few other CNVs that are less common because they are
#' ## very large, or because they occur in regions that in which copy
#' ## number alerations are not common
#' ##
#' cnv3 <- GRanges("chr1", IRanges(9e6L, 15e6L), id="sample1400")
#' starts <- seq(5e6L, 200e6L, 10e6L)
#' ends <- starts + rpois(length(starts), 25e3L)
#' cnv4 <- GRanges("chr1", IRanges(starts, ends),
#'                 id=paste0("sample", sample(1000:1500, length(starts))))
#'
#' all_cnvs <- suppressWarnings(c(cnv1, cnv2, cnv3, cnv4))
#' grl <- split(all_cnvs, all_cnvs$id)
#' cnps <- consensusCNP(grl)
#'
#' ##
#' ## 2nd CNP is filtered because of its size
#' ##
#' truth <- GRanges("chr1", IRanges(10000100L, 10100100L))
#' seqinfo(truth) <- seqinfo(grl)
#' identical(cnps, truth)
#'
#' ##
#' ## Both CNVs identified
#' ##
#' cnps <- consensusCNP(grl, max.width=500e3)
#' truth <- GRanges(c("chr1", "chr5"),
#'                  IRanges(c(10000100L, 101000999L),
#'                          c(10100100L, 101400999L)))
#' seqlevels(truth, pruning.mode="coarse") <- seqlevels(grl)
#' seqinfo(truth) <- seqinfo(grl)
#' identical(cnps, truth)
#'
#' @param grl  A \code{GRangesList} of all CNVs in a study -- each
#' element is the collection of CNVs for one individual.
#' @param transcripts a \code{GRanges} object containing annotation of
#' genes or transcripts (optional)
#' @param min.width length-one integer vector specifying the minimum width of CNVs
#' @param max.width length-one integer vector specifying the maximum
#' width of CNVs
#' @param min.prevalance a length-one numeric vector specifying the
#' minimum prevalance of a copy number polymorphism.  Must be in the
#' interval [0,1].  If less that 0, this function will return all CNV
#' loci regardless of prevalance.  If greater than 1, this function
#' will return a length-zero GRanges object
#' @return a \code{GRanges} object providing the intervals of all
#' identified CNPs above a user-specified prevalance cutoff.
#' @export
consensusCNP <- function(grl, transcripts, min.width=2e3, max.width=200e3, min.prevalance=0.02){
  g <- as(unlist(grl), "GRanges")
  si <- seqinfo(g)
  names(g) <- NULL
  grl <- split(g, g$id)
  ##grl <- grl[colnames(views)]
  regions <- defineCnpRegions(grl, thr=min.prevalance)
  filters <- width(regions) > min.width  & width(regions) <= max.width
  if(any(filters)){
    message("Dropping CNV regions failing min.width and max.width criteria. See ?consensusCNP to relax these settings.")
    regions <- regions[ filters ]
  }
  regions <- reduce(regions)
  if(!missing(transcripts)){
    regions <- annotateRegions(regions, transcripts)
  }
  seqlevels(regions, pruning.mode="coarse") <- seqlevels(si)
  seqinfo(regions) <- si
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


imputeFromSampledData <-  function(model, data, index){
  if(is.null(names(data))) stop("data must be a named vector")
  pz2 <- probz(model)
  cn <- map(model)
  ##pz2 <- mapCnProbability(pz, k(model))
  ##r2 <- quantile(data, probs=seq(0, 1, by=0.01))
  df <- data.frame(r=data,
                   cn=rep(NA, length(data)),
                   ##p=rep(NA, length(data)),
                   row.names=names(data))
  ##df$p[index] <- pz2
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
  ##p <- setNames(df.complete$p, df.complete$quantiles2)
  ##df.incomplete$p <- p[as.character(df.incomplete$quantiles2)]
  df2 <- rbind(df.incomplete, df.complete)
  df2 <- df2[rownames(df), ]
  df2
}

permnK <- function(k, maxperm){
  if(k < 2) return(matrix(1,1,1))
  kperm <- permn(seq_len(k))
  kperm <- do.call("rbind", kperm)
  kperm.identity <- kperm[1, , drop=FALSE]
  kperm <- kperm[-1, , drop=FALSE]
  neq <- apply(kperm, 1, function(x, y) sum(x != y), y=kperm.identity)
  kperm <- kperm[order(neq, decreasing=TRUE), , drop=FALSE]
  N <- min(maxperm-1, nrow(kperm))
  kperm <- rbind(kperm.identity, kperm[seq_len(N), ])
  kperm
}


#' Create tile labels for each observation
#'
#' @param y vector containing data
#' @param nt the number of tiles in a batch
#' @param batch a vector containing the labels from which batch each observation came from.
#' @return Tile labels for each observation
#' @export
downSampleEachBatch <- function(y, nt, batch){
  ## NULL out these two variables to avoid NOTE about
  ## no visible binding for global variable
  x <- obs.index <- NULL
  yb <- split(y, batch)
  indices <- split(seq_along(y), batch)
  S <- vector("list", length(yb))

  for (i in 1:length(yb)) {
    x <- yb[[i]]
    batch.id <- names(yb)[i]
    obs.index <- indices[[i]]

    tiles <- factor(paste(ntile(x, nt), batch.id, sep="_"))
    xx <- sapply(split(x, tiles), mean)
    tiles <- setNames(as.character(tiles), obs.index)
    S[[i]] <- list(y=xx, 
                   labels=tiles, 
                   batch.id=rep(batch.id, length(xx)))
  }

  ys <- unlist(lapply(S, "[[", 1))
  tiles <- unlist(lapply(S, "[[", 2))
  batch.id <- unlist(lapply(S, "[[", 3))
  ##
  ## tile labels must be in the same order as the original y vector
  ##
  tiles <- tiles[order(as.numeric(names(tiles)))]
  list(y=ys, labels=tiles, batch=batch.id)
}

#' Create tile labels for each observation
#'
#' A wrapper for function downSampleEachBatch. Batches are automatically merged as needed.
#' @param batch.file the name of a file contaning RDS data to be read in.
#' @param plate a vector containing the labels  from which batch each observation came from.
#' @param y in memory data
#' @param ntiles number of tiles in a batch
#' @param THR threshold above which to merge batches in Kolmogorov-Smirnov test.
#' @return Tile labels for each observation
#' @export
downsample <- function(batch.file, plate, y, ntiles=250, THR=0.1){
  if(file.exists(batch.file)){
    batches <- readRDS(batch.file)
  } else {
    message("merging plates... ")
    if(missing(plate)) stop("batch.file does not exist.  Chemistry plate must be specified")
    batches <- collapseBatch(y, plate, THR=THR)
    saveRDS(batches, file=batch.file)
  }
  ds <- downSampleEachBatch(y, ntiles, batches)
  if(length(unique(ds$batch)) > 12){
    batches <- collapseBatch(ds$y, ds$batch, THR=1e-3)
    ds$batch <- batches
  }
  if(any(table(ds$batch) <= 25)){
    tab <- table(ds$batch)
    keep <- ds$batch %in% names(tab)[tab > 25]
    ds$y <- ds$y[keep]
    ds$batch <- ds$batch[keep]
  }
  ds
}

# This function is copied from dplyr! Just copied the code over since this
# is the only function used in dplyr.
ntile <- function(x, n) {
    floor((n * (rank(x, ties.method="first", na.last="keep") - 1)
           / length(x)) + 1)
}

#' Calculate posterior proportion of cases by component
#'
#' @examples
#'      # generate random case control status
#'      case_control <- rbinom(length(y(MarginalModelExample)), 1, 0.5)
#'      case_control_posterior <- posterior_cases(MarginalModelExample,
#'                                                case_control)
#' @param model An instance of a \code{MixtureModel}-derived class.
#' @param case_control A vector of 1's and 0's where a 1 indicates a case and a 0 a control
#' @param alpha prior alpha for the beta
#' @param beta prior beta for the beta
#' @return A matrix of dimension S (MCMC iterations) by K (number of components) where each element i,j indicates the posterior proportion of cases at an iteration and component
#' @export
posterior_cases <- function(model, case_control, alpha=1, beta=1) {
    # model and MCMC params
    z.mat <- z(chains(model))
    S <- nrow(z.mat)

    # empty posterior matrix
    posterior <- matrix(0, nrow=S, ncol=k(model))

    # run model of probabilities
    for (i in seq_len(S)) {
        z <- z.mat[i, ]
        cont.tbl <- table(z, case_control)
        
        cases <- cont.tbl[, 2]
        controls <- cont.tbl[, 1]

        for (j in seq_len(length(cases))) {
            cases.row <- cases[j]
            controls.row <- controls[j]
            posterior[i, j] <- rbeta(1, cases.row + 1, controls.row + 1)
        }
    }

    return(posterior)
}
