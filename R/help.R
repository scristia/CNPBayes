#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @importFrom utils capture.output
#' @importFrom scales rescale
#' @import GenomicRanges
#' @importFrom Rcpp evalCpp
#' @import S4Vectors
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSds colVars colMedians rowProds colMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm qnorm kmeans ks.test plot.ts qgamma rbeta rgamma dbeta
#' @importFrom stats rgeom rnorm runif setNames rpois rchisq dgamma df
#' @importFrom coda effectiveSize mcmc.list gelman.diag mcmc as.mcmc.list
#' @importFrom mclust Mclust mclustBIC
#' @importFrom reshape2 melt
#' @importMethodsFrom SummarizedExperiment assays SummarizedExperiment assays<-
#' @import ggplot2
#' @importFrom magrittr set_colnames "%>%"
#' @importFrom purrr map map2 map_lgl map_dbl map_chr
#' @importFrom tidyr gather spread
#' @importFrom dplyr mutate bind_rows group_by summarize arrange filter n left_join select ungroup
#' @import tibble
#' @useDynLib CNPBayes
NULL
