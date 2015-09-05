#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom gtools rdirichlet ddirichlet
#' @import Rcpp
#' @import RcppArmadillo
#' @import GenomeInfoDb
#' @import BiocGenerics
#' @import IRanges
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% %dopar% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom oligoClasses copyNumber batch chromosome integerMatrix
#' @importMethodsFrom oligoClasses lrr baf
#' @useDynLib CNPBayes
NULL
