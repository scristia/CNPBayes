#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods BiocGenerics S4Vectors IRanges GenomicRanges 
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @importFrom Rcpp evalCpp
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importMethodsFrom oligoClasses lrr baf
#' @useDynLib CNPBayes
NULL
