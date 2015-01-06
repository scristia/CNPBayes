#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @import GenomicRanges
#' @import mixtools
#' @importFrom gtools rdirichlet
#' @importFrom msm rtnorm
#' @importFrom truncnorm rtruncnorm
#' @import Hmisc
#' @import Rcpp
#' @import RcppArmadillo
#' @import GenomeInfoDb
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom matrixStats colSds
#' @importFrom HardyWeinberg HWChisq
#' @useDynLib CNPBayes
NULL

## what are we iimporting in gtools?
## HMisc rMultinom
## gtools rdirichlet
