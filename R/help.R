#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @importFrom utils capture.output
#' @import GenomicRanges
#' @importFrom Rcpp evalCpp
#' @import S4Vectors
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSds colVars colMedians rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm qnorm kmeans ks.test plot.ts qgamma rbeta rgamma rgeom rnorm runif setNames rpois
#' @importFrom coda effectiveSize mcmc.list gelman.diag mcmc as.mcmc.list
#' @importFrom mclust Mclust mclustBIC
#' @importFrom reshape2 melt
#' @importMethodsFrom SummarizedExperiment assays SummarizedExperiment
#' @import ggplot2
#' @importFrom magrittr set_colnames "%>%"
#' @importFrom purrr map map2 map_lgl map_dbl map_chr
#' @importFrom tidyr gather
#' @importFrom dplyr mutate bind_rows group_by summarize arrange
#' @importFrom tibble as.tibble
#' @useDynLib CNPBayes
NULL
