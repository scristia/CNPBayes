#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @importFrom Rcpp evalCpp
#' @importFrom S4Vectors SimpleList DataFrame Rle
#' @import GenomicRanges 
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom foreach %do% foreach
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm kmeans ks.test plot.ts qgamma rbeta rgamma rgeom rnorm runif setNames
#' @importMethodsFrom oligoClasses lrr baf
#' @useDynLib CNPBayes
NULL
