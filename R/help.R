#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom msm rtnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom Hmisc rMultinom
#' @import Rcpp
#' @import RcppArmadillo
#' @import GenomeInfoDb
#' @import BiocGenerics
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom HardyWeinberg HWChisq
#' @importFrom oligoClasses copyNumber batch chromosome integerMatrix
#' @importMethodsFrom oligoClasses lrr
#' @importFrom coda effectiveSize
#' @importFrom combinat permn
#' @importFrom dplyr ntile
#' @useDynLib CNPBayes
NULL

## what are we iimporting in gtools?
## HMisc rMultinom
## gtools rdirichlet ddirichlet

## S4Vectors is needed for the  SimpleList constructor
#  mixtools
