#' Bayesian mixture models for copy number estimation
#'
#' @docType package
#' @name CNPBayes
#' @import methods
#' @importFrom gtools rdirichlet ddirichlet permutations smartbind
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqlevels seqinfo<- seqnames
#' @import IRanges
#' @importFrom utils capture.output
#' @importFrom scales rescale
#' @import GenomicRanges
#' @importFrom Rcpp evalCpp
#' @import S4Vectors
#' @importFrom combinat permn
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSds colVars colMedians rowProds colMaxs rowMaxs
#' @importFrom BiocGenerics unlist
#' @importFrom graphics lines par
#' @importFrom stats dnorm qnorm kmeans ks.test plot.ts qgamma rbeta rgamma dbeta pnorm prcomp rbinom t.test
#' @importFrom stats rgeom rnorm runif setNames rpois rchisq dgamma df
#' @importFrom coda effectiveSize mcmc.list gelman.diag mcmc as.mcmc.list
#' @import SummarizedExperiment
#' @import ggplot2
#' @importFrom magrittr set_colnames "%>%" "%$%"
#' @importFrom purrr map map2 map_lgl map_dbl map_chr
#' @importFrom tidyr gather spread complete
#' @importFrom stringr str_pad
#' @importFrom dplyr mutate bind_rows group_by summarize arrange filter n pull left_join select ungroup rowwise slice_sample
#' @import grid
#' @import tibble
#' @useDynLib CNPBayes
NULL
