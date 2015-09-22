#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet
#' @import Rcpp
#' @import GenomeInfoDb
#' @import IRanges
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom GenomicRanges rowRanges SummarizedExperiment colData GRanges GRangesList assays
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom oligoClasses copyNumber batch chromosome integerMatrix
#' @importMethodsFrom oligoClasses lrr baf
#' @useDynLib CNPBayes
NULL
