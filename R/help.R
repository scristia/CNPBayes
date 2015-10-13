#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @importFrom Rcpp evalCpp
#' @importFrom S4Vectors SimpleList DataFrame
#' @import GenomicRanges 
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importMethodsFrom oligoClasses lrr baf
#' @useDynLib CNPBayes
NULL
