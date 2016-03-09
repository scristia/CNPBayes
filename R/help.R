#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods BiocGenerics S4Vectors IRanges GenomicRanges 
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @importFrom Rcpp evalCpp
#' @importFrom S4Vectors SimpleList DataFrame Rle
#' @import GenomicRanges 
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm kmeans ks.test plot.ts qgamma rbeta rgamma rgeom rnorm runif setNames
#' @importFrom coda effectiveSize
#' @useDynLib CNPBayes
NULL
