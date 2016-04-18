#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @importFrom Rcpp evalCpp
#' @import S4Vectors
#' @import GenomicRanges 
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSds colMedians rowCumsums rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm kmeans ks.test plot.ts qgamma rbeta rgamma rgeom rnorm runif setNames
#' @importFrom coda effectiveSize
#' @importMethodsFrom SummarizedExperiment assays SummarizedExperiment
#' @useDynLib CNPBayes
NULL
