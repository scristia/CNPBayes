% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiBatch.R
\name{revertToMultiBatch}
\alias{revertToMultiBatch}
\title{Convert a single-batch object to a multi-batch model}
\usage{
revertToMultiBatch(sb)
}
\arguments{
\item{sb}{a `MultiBatch` object}
}
\description{
We often evaluate single-batch models even when `kolmogorov_batches` identifies multiple batches.  In particular, if the batch effects are modest and the overlap of mixture components is minimal in the single-batch model, the more parsimonious single-batch model is preferred.  We can convert the single-batch instance back to a multi-batch instance by using the `batch-labels` field. The vignette provides an example.
}
