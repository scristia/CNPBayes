#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @import GenomicRanges
#' @import S4Vectors
#' @import mixtools
#' @importFrom gtools rdirichlet
#' @importFrom msm rtnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom Hmisc rMultinom
#' @import Rcpp
#' @import RcppArmadillo
#' @import GenomeInfoDb
#' @import BiocGenerics
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds
#' @importFrom HardyWeinberg HWChisq
#' @importFrom oligoClasses copyNumber batch chromosome integerMatrix
#' @importFrom coda effectiveSize
#' @importFrom combinat permn
#' @useDynLib CNPBayes
NULL

## what are we iimporting in gtools?
## HMisc rMultinom
## gtools rdirichlet

## S4Vectors is needed for the  SimpleList constructor
