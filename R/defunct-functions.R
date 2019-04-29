#' @include Defunct-classes.R
NULL

doKmeans <- function(thetas, modes) length(thetas) > length(modes) && length(modes) > 1

findModes <- function(quantiles, x){ ##quantiles, density
  signs <- sign(diff(x))
  signs <- c(signs[1], signs)
  inflect <- c(0, cumsum(diff(signs) < 0))
  indices <- tapply(seq_along(inflect), inflect, max)
  indices <- indices[-length(indices)]
  modes <- quantiles[indices]
  modes
}

#' Create tile labels for each observation
#'
#' For large datasets (several thousand subjects), the computational burden for fitting Bayesian mixture models can be high.  Downsampling can reduce the computational burden with little effect on inference.  The function tileMedians is useful for putting the median log R ratios for each subject in a bucket. The observations in each bucket are averaged.  This is done independently for each batch and the range of median log R ratios within each bucket is guaranteed to be less than 0.05.  Note this function requires specification of a batch variable. If the study was small enough such that all the samples were processed in a single batch, then downsampling would not be needed.  By summarizing the observations in each bucket by batch, the SingleBatchModels (SB or SBP) and MultiBatchModels (MB or MBP) will be fit to the same data and are still comparable by marginal likelihoods or Bayes Factors.
#'
#' @param y vector containing data
#' @param nt the number of observations per batch
#' @param batch a vector containing the labels from which batch each observation came from.
#' @return A tibble with a tile assigned to each log R ratio
#' @seealso \code{\link[dplyr]{ntile}}
#' @export
#' @rdname tile-functions
tileMedians <- function(y, nt, batch){
  .Defunct("See downSample")
  if(any(is.na(y))) stop("NAs found in y-vector")
  yy <- y
  logratio <- x <- obs.index <- NULL
  ##
  ## split observations by batch
  ##
  yb <- split(y, batch)
  indices <- split(seq_along(y), batch)
  S <- vector("list", length(yb))
  MAX_TILE <- 0
  for (i in seq_along(yb)) {
    x <- yb[[i]]
    batch.id <- names(yb)[i]
    obs.index <- indices[[i]]
    ##
    ## observations can be quite different within a tile (e.g., homozygous deletions)
    ##
    tiles <- ntile(x, nt) %>% as_tibble %>%
      mutate(x=x, obs.index) %>%
      set_colnames(c("tile", "logratio", "obs.index"))
    tiles2 <- tiles %>%
      group_by(tile) %>%
      summarize(spread=abs(diff(range(logratio))),
                n=n()) %>%
      arrange(-spread)
    ## Split tiles with large spread into multiple tiles
    tiles.keep <- tiles2 %>%
      filter(spread < 0.01)
    tiles.drop <- tiles2 %>%
      filter(spread >= 0.01)
    newtiles <- tiles %>% filter(tile %in% tiles.drop$tile)
    newtiles$tile <- seq_len(nrow(newtiles)) + max(tiles$tile)
    tiles3 <- filter(tiles, tile %in% tiles.keep$tile) %>%
      bind_rows(newtiles) %>%
      mutate(tile=tile + MAX_TILE)
    tiles3$batch.var <- batch.id
    ##tiles3$tile <- paste(batch.id, tiles3$tile, sep="_")
    S[[i]] <- tiles3
    MAX_TILE <- max(tiles3$tile)
  }
  tiles <- do.call(bind_rows, S)
  tiles$batch <- as.integer(factor(tiles$batch.var))
  tiles
}

