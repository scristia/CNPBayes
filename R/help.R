#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @import GenomicRanges
#' @importFrom Rcpp evalCpp
#' @import S4Vectors
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSds colVars colMedians rowCumsums rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm kmeans ks.test plot.ts qgamma rbeta rgamma rgeom rnorm runif setNames
#' @importFrom coda effectiveSize
#' @importFrom mclust Mclust mclustBIC
#' @importFrom reshape2 melt
#' @importMethodsFrom SummarizedExperiment assays
#' @useDynLib CNPBayes
NULL
