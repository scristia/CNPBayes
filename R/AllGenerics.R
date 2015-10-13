#' @include help.R
NULL

#' Number of components.
#'
#' This function retrieves the number of a priori components.
#' @param object see \code{showMethods(k)}
#' @return The number of components
#' @export
#' @docType methods
#' @rdname k-method
setGeneric("k", function(object) standardGeneric("k"))

#' Change the number of components.
#'
#' Updates the number of components and erases chains from a previous
#' posteriorSimulation (if one was performed). Draws from prior to guess
#' new starting values.
#' @examples
#' k(MarginalModelExample) <- 2
#' @param value An integer for the new number of components.
#' @export
#' @docType methods
#' @rdname k-method
setGeneric("k<-", function(object, value) standardGeneric("k<-"))

setGeneric("z<-", function(object, value) standardGeneric("z<-"))

setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))
setGeneric("sigma2<-", function(object, value) standardGeneric("sigma2<-"))
setGeneric("p<-", function(object, value) standardGeneric("p<-"))

#' Retrieve overall mean
#'
#' @examples
#'      mu(MarginalModelExample)
#' @param object see \code{showMethods(mu)}
#' @return A vector containing 'mu'
#' @export
#' @docType methods
#' @rdname mu-method
setGeneric("mu", function(object) standardGeneric("mu"))
setGeneric("mu<-", function(object, value) standardGeneric("mu<-"))

#' Accessor for the tau2 parameter in the hierarchical mixture model
#'
#' The interpretation of \code{tau2} depends on whether \code{object}
#' is a \code{MarginalModel} or a \code{BatchModel}.  For
#' \code{BatchModel}, \code{tau2} is a vector with length equal to the
#' number of components.  Each element of the \code{tau2} vector can
#' be interpreted as the within-component variance of the batch means
#' (\code{theta}).  For objects of class \code{MarginalModel}
#' (assumes no batch effect), \code{tau2} is a length-one vector that
#' describes the variance of the component means between batches.  The
#' hyperparameters of \code{tau2} are \code{eta.0} and \code{m2.0}. See the
#' following examples for setting the hyperparameters, accessing the current
#' value of \code{tau2} from a \code{MixtureModel}-derived object, and
#' for plotting the chain of \code{tau2} values.
#'
#' @examples
#' k(BatchModelExample)
#' tau2(BatchModelExample)
#' plot.ts(tau2(chains(BatchModelExample)))
#'
#' @param object see \code{showMethods(tau2)}
#' @return A vector of variances
#' @export
#' @docType methods
#' @rdname tau2-method
#' @seealso \code{Hyperparameters}
setGeneric("tau2", function(object) standardGeneric("tau2"))
setGeneric("tau2<-", function(object, value) standardGeneric("tau2<-"))
setGeneric("nu.0<-", function(object, value) standardGeneric("nu.0<-"))
setGeneric("sigma2.0<-", function(object, value) standardGeneric("sigma2.0<-"))
setGeneric("dataMean<-", function(object, value) standardGeneric("dataMean<-"))
setGeneric("dataPrec<-", function(object, value) standardGeneric("dataPrec<-"))

setGeneric("chains<-", function(object, value) standardGeneric("chains<-"))

#' Retrieve simulated chains from model object.
#'
#' The method \code{chains} applied to a \code{MixtureModel}-derived
#' class will return an object of class \code{McmcChains} that
#' contains the chains for all simulated parameters. Typically, chains
#' is called in conjunction with an accessor for one of these
#' parameters.
#' @examples
#' theta.chain <- theta(chains(MarginalModelExample))
#' dim(theta.chain)
#' plot.ts(theta.chain, plot.type="single",
#'         col=seq_len(k(MarginalModelExample)))
#' @param object \code{showMethods(chains)}
#' @return The simulated chains.
#' @export
#' @docType methods
#' @rdname chains-method
setGeneric("chains", function(object) standardGeneric("chains"))

#' Accessor for Hyperparameters object for a MixtureModel-derived object
#'
#' @examples
#' \dontrun{
#'     hyperParams(MarginalModelExample)
#' }
#' @param object see \code{showMethods(hyperParams)}
#' @return The Hyperparameters of a MixtureModel
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams", function(object) standardGeneric("hyperParams"))

#' Replace the hyperparameters for a \code{MixtureModel}-derived object
#'
#' @examples
#' hypp <- Hyperparameters(type="marginal",
#'                         k=k(MarginalModelExample),
#'                         alpha=c(9, 9, 10))
#' hyperParams(MarginalModelExample) <- hypp
#'
#' @param value an object of class 'Hyperparameters'
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams<-", function(object,value) standardGeneric("hyperParams<-"))

setGeneric("McmcChains", function(object) standardGeneric("McmcChains"))

setGeneric("hist")

#' Plot the densities estimated from a mixture model for a copy number polymorphism
#'
#' Plot estimates of the posterior density for each component and the
#' overall, marginal density.  For batch models, one can additionally
#' plot batch-specific density estimates.
#'
#' @param x a \code{DensityModel}-derived object, or a
#' \code{MixtureModel}-derived object.
#' @param y If \code{x} is a \code{DensityModel}, \code{y} is a
#' numeric vector of the one-dimensional summaries for a given copy
#' number polymorphism. If \code{x} is a \code{MixtureModel}, \code{y}
#' is ignored.
#' @param show.batch a logical. If true, batch specific densities
#' will be plotted.
#' @param ... Additional arguments passed to \code{hist}.
#' @return A plot showing the density estimate
#' @examples
#'   set.seed(100)
#'   truth <- simulateData(N=2500,
#'                         theta=c(-2, -0.4, 0),
#'                         sds=c(0.3, 0.15, 0.15),
#'                         p=c(0.05, 0.1, 0.8))
#'
#'   mcmcp <- McmcParams(iter=500, burnin=500, thin=2)
#'   model <- MarginalModel(y(truth), k=3, mcmc.params=mcmcp)
#'   model <- CNPBayes:::startAtTrueValues(model, truth)
#'   model <- posteriorSimulation(model)
#'   par(mfrow=c(1,2), las=1)
#'   plot(truth)
#'   plot(model)
#' @export
setGeneric("plot")


#' Retrieve batches from object.
#'
#' The batches are represented as a vector of integers.
#' @examples
#'      batch(BatchModelExample)
#' @param object see \code{showMethods(batch)}
#' @return The batch of each data element.
#' @export
#' @docType methods
#' @rdname batch-method
setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("batch<-", function(object,value) standardGeneric("batch<-"))

setGeneric("startingValues", function(object, params, zz) standardGeneric("startingValues"))

setGeneric("computeMeans", function(object) standardGeneric("computeMeans"))
setGeneric("computeVars", function(object) standardGeneric("computeVars"))

setGeneric("computePotential", function(object) standardGeneric("computePotential"))

setGeneric("initializeSigma2", function(object) standardGeneric("initializeSigma2"))

setGeneric("initializeTheta", function(object) standardGeneric("initializeTheta"))

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))

#' Calculate BIC of a model
#'
#' @examples
#'      bic(BatchModelExample)
#' @param object see \code{showMethods(bic)}
#' @return The BIC of the model.
#' @docType methods
#' @rdname bic-method
#' @export
setGeneric("bic", function(object) standardGeneric("bic"))

#' Accessor for the theta parameter in the hierarchical mixture model
#'
#' The interpretation of \code{theta} depends on whether \code{object}
#' is a \code{MarginalModel} or a \code{BatchModel}.  For
#' \code{BatchModel}, \code{theta} is a matrix of size B x K, where B is
#' the number of batches and K is the number of components.
#' Each column of the \code{theta} matrix can be interpreted as the
#' batch means for a particular component. For objects of class
#' \code{MarginalModel} (assumes no batch effect), \code{theta} is a
#' vector of length K. Each element of \code{theta} can be interpreted
#' as the mean for a component. See the following examples for accessing
#' the current value of \code{theta} from a \code{MixtureModel}-derived
#' object, and for plotting the chain of \code{theta} values.
#' @examples
#' ## MarginalModel
#' k(MarginalModelExample)
#' theta(MarginalModelExample)
#' plot.ts(theta(chains(MarginalModelExample)))
#' ## BatchModel
#' k(BatchModelExample)
#' length(unique(batch(BatchModelExample)))
#' theta(BatchModelExample)
#' ## Plot means for batches in one component
#' plot.ts(theta(chains(BatchModelExample))[, 1:3])
#' @param object see \code{showMethods(theta)}
#' @return A vector of length number of components or a matrix of size 
#' number of batches x number of components
#' @export
#' @docType methods
#' @rdname theta-method
setGeneric("theta", function(object) standardGeneric("theta"))

#' Retrieve the variances of each component and batch distribution
#'
#' For a MarginalModel, this function returns a vector of variances. For a BatchModel, returns a matrix of size number of batches by number of components.
#' @examples
#'      sigma2(MarginalModelExample)
#' @param object see \code{showMethods(sigma2)}
#' @return A vector of length number of components or a matrix of size 
#' number of batches x number of components
#' @export
#' @docType methods
#' @rdname sigma2-method
setGeneric("sigma2", function(object) standardGeneric("sigma2"))

#' Retrieve the probability of latent variable membership by observation.
#'
#' @examples
#'      probz(MarginalModelExample)
#' @param object see \code{showMethods(probz)}
#' @return A matrix of size number of observations x number of components
#' @export
#' @docType methods
#' @rdname probz-method
setGeneric("probz", function(object) standardGeneric("probz"))

setGeneric("probz<-", function(object, value) standardGeneric("probz<-"))

#' Retrieve the shape parameter for the sigma.2 distribution.
#'
#' @examples
#'      nu.0(MarginalModelExample)
#' @param object see \code{showMethods(nu.0)}
#' @return An integer
#' @export
#' @docType methods
#' @rdname nu.0-method
setGeneric("nu.0", function(object) standardGeneric("nu.0"))

#' Retrieve the rate parameter for the sigma.2 distribution.
#'
#' @examples
#'      sigma2.0(MarginalModelExample)
#' @param object see \code{showMethods(sigma2.0)}
#' @return A length 1 numeric
#' @export
#' @docType methods
#' @rdname sigma2.0-method
setGeneric("sigma2.0", function(object) standardGeneric("sigma2.0"))

#' Retrieve data.
#'
#' @examples
#'      y(MarginalModelExample)
#' @param object see \code{showMethods(y)}
#' @return A vector containing the data
#' @export
#' @docType methods
#' @rdname y-method
setGeneric("y", function(object) standardGeneric("y"))

setGeneric("y<-", function(object, value) standardGeneric("y<-"))

## TODO: oned is for 1-dimensional summary. Perhaps more informative
## than 'y'.

#' Retrieve data.
#'
#' @param object see \code{showMethods(oned)}
#' @return A vector the length of the data
#' @export
#' @docType methods
#' @rdname oned-method
#' @export
setGeneric("oned", function(object) standardGeneric("oned"))

#' Retrieve latent variable assignments.
#'
#' Retrieves the simulated latent variable assignments of each observation at each MCMC simulation.
#' @examples
#'      z(MarginalModelExample)
#' @param object see \code{showMethods(z)}
#' @return A vector the length of the data
#' @export
#' @docType methods
#' @rdname z-method
setGeneric("z", function(object) standardGeneric("z"))

setGeneric("modalParameters", function(object) standardGeneric("modalParameters"))

setGeneric("computeModes", function(object) standardGeneric("computeModes"))

#' Retrieve the modes from a model.
#'
#' The iteration which maximizes log likelihood and log prior is found. The estimates for each parameter at this iteration are retrieved.
#' @examples
#'      modes(MarginalModelExample)
#' @param object a \code{MixtureModel}-derived class
#' @return A list of the modes of each parameter
#' @export
#' @docType methods
#' @rdname modes-method
setGeneric("modes", function(object) standardGeneric("modes"))


#' Replacement method for modes
#'
#' For a mixture model with K components, there are K! possible modes.
#' One can permute the ordering of the modes and assign the permuted
#' order to a MixtureModel derived class by this method.
#'
#' @param value a \code{list} of the modes.  See \code{mode(object)}
#' to obtain the correct format of the list.
#'
#' @export
#' @docType methods
#' @rdname modes-method
setGeneric("modes<-", function(object,value) standardGeneric("modes<-"))


setGeneric("mu.0", function(object) standardGeneric("mu.0"))
setGeneric("tau2.0", function(object) standardGeneric("tau2.0"))

#' Retrieve the rate parameter for the tau2 distribution.
#'
#' @examples
#'      eta.0(MarginalModelExample)
#' @param object see \code{showMethods(eta.0)}
#' @return eta.0 of a 'MixtureModel'
#' @export
#' @docType methods
#' @rdname eta.0-method
setGeneric("eta.0", function(object) standardGeneric("eta.0"))
setGeneric("eta.0<-", function(object,value) standardGeneric("eta.0<-"))

#' Retrieve the shape parameter for the tau2 distribution.
#'
#' @examples
#'      m2.0(MarginalModelExample)
#' @param object see \code{showMethods(m2.0)}
#' @return m2.0 for a model
#' @export
#' @docType methods
#' @rdname m2.0-method
setGeneric("m2.0", function(object) standardGeneric("m2.0"))

setGeneric("m2.0<-", function(object,value) standardGeneric("m2.0<-"))

setGeneric("showMeans", function(object) standardGeneric("showMeans"))
setGeneric("showSigmas", function(object) standardGeneric("showSigmas"))

## JC this function could use a better name. See examples
##
#' Estimate batch from a collection of chemistry plates or some other
#' variable that captures the time in which the arrays were processed.
#'
#' In high-throughput assays, low-level summaries of copy number at
#' copy number polymorphic loci (e.g., the mean log R ratio for each
#' sample, or a principal-component derived summary) often differ
#' between groups of samples due to technical sources of variation
#' such as reagents, technician, or laboratory.  Technical (as opposed
#' to biological) differences between groups of samples are referred
#' to as batch effects.  A useful surrogate for batch is the chemistry
#' plate on which the samples were hybridized. In large studies, a
#' Bayesian hierarchical mixture model with plate-specific means and
#' variances is computationally prohibitive.  However, chemistry
#' plates processed at similar times may be qualitatively similar in
#' terms of the distribution of the copy number summary statistic.
#' Further, we have observed that some copy number polymorphic loci
#' exhibit very little evidence of a batch effect, while other loci
#' are more prone to technical variation.  We suggest combining plates
#' that are qualitatively similar in terms of the Kolmogorov-Smirnov
#' two-sample test of the distribution and to implement this test
#' independently for each candidate copy number polymophism identified
#' in a study.  The \code{collapseBatch} function is a wrapper to the
#' \code{ks.test} implemented in the \code{stats} package that
#' compares all pairwise combinations of plates.  The \code{ks.test}
#' is performed recursively on the batch variables defined for a given
#' CNP until no batches can be combined.
#' @examples
#' bt <- collapseBatch(y(BatchModelExample), batch(BatchModelExample))
#' newBatchModel <- BatchModel(y(BatchModelExample), k(BatchModelExample),
#'                             bt, hyperParams(BatchModelExample),
#'                             mcmcParams(BatchModelExample))
#' @param object see \code{showMethods(collapseBatch)}
#' @param plate a vector labelling from which batch each observation came from.
#' @param THR threshold below which the null hypothesis should be rejected and batches are collapsed.
#' @return The new batch value.
#' @export
#' @docType methods
#' @rdname collapseBatch-method
setGeneric("collapseBatch", function(object, plate, THR=0.1) standardGeneric("collapseBatch"))

setGeneric("thetac", function(object) standardGeneric("thetac"))

setGeneric("thetaMean", function(object) standardGeneric("thetaMean"))

setGeneric("sigmaMean", function(object) standardGeneric("sigmaMean"))

setGeneric("pMean", function(object) standardGeneric("pMean"))

#' Create a trace plot of a parameter estimated by MCMC.
#'
#' @examples
#' tracePlot(BatchModelExample, "theta")
#' tracePlot(BatchModelExample, "sigma")
#' @param object see \code{showMethods(tracePlot)}
#' @param name the name of the parameter for which to plot values. Can be 'theta', 'sigma', 'p', 'mu', or 'tau'.
#' @param ... Other argument to pass to plot.
#' @return A traceplot of a parameter value
#' @export
#' @docType methods
#' @rdname tracePlot-method
setGeneric("tracePlot", function(object, name, ...) standardGeneric("tracePlot"))

setGeneric("tablez", function(object) standardGeneric("tablez"))

#' Number of MCMC chains.
#'
#' This function retrieves the number of chains used for an MCMC simulation.
#' @examples
#' number_of_chains <- nStarts(MarginalModelExample)
#' @param object see \code{showMethods(nStarts)}
#' @return An integer of the number of different starts.
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts", function(object) standardGeneric("nStarts"))

#' Reset number of MCMC chains in simulation.
#'
#' This function changes the number of chains used for an MCMC simulation.
#' @examples
#' number_of_chains <- 3
#' nStarts(MarginalModelExample) <- number_of_chains
#' @param value new number of chains
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts<-", function(object, value) standardGeneric("nStarts<-"))

setGeneric("alpha", function(object) standardGeneric("alpha"))

#' Retrieve log likelihood.
#'
#' @examples
#' ## retrieve log likelihood at each MCMC iteration
#' log_lik(chains(MarginalModelExample))
#' ## retrieve log likelihood at last MCMC iteration
#' log_lik(MarginalModelExample)
#' @param object see showMethods(log_lik)
#' @return The log likelihood
#' @export
#' @docType methods
#' @rdname log_lik-method
setGeneric("log_lik", function(object) standardGeneric("log_lik"))

setGeneric("log_lik<-", function(object,value) standardGeneric("log_lik<-"))

setGeneric("computeLoglik", function(object) standardGeneric("computeLoglik"))

#' Number of burnin iterations.
#'
#' This function retrieves the number of burnin simulations to be discarded.
#' @examples
#' burnin(MarginalModelExample)
#' mp <- mcmcParams(MarginalModelExample)
#' burnin(mp)
#' @param object see \code{showMethods(burnin)}
#' @return The number of burnin simulations.
#' @export
#' @docType methods
#' @rdname burnin-method
setGeneric("burnin", function(object) standardGeneric("burnin"))

#' Reset number of burnin iterations.
#'
#' This function changes the number of burnin simulations to be discarded.
#' @param value new number of burnin iterations
#' @export
#' @docType methods
#' @rdname burnin-method
setGeneric("burnin<-", function(object, value) standardGeneric("burnin<-"))

#' Reset number of iterations.
#'
#' This function changes the number of simulations.
#' @param force Allow changing of the size of the elements?
#' @param value new number of iterations
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter<-", function(object, force=FALSE, value) standardGeneric("iter<-"))

#' Number of MCMC iterations.
#'
#' This function retrieves the number of iterations of an MCMC simulation.
#' @examples
#'      iter(MarginalModelExample)
#' @param object see \code{showMethods(iter)}
#' @return The number of MCMC iterations
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter", function(object) standardGeneric("iter"))

#' Number of thinning intervals.
#'
#' This function retrieves the number of thinning intervals used for an MCMC simulation.
#' @examples
#'      thin(MarginalModelExample)
#' @param object see showMethods(thin)
#' @return An integer of the number of thinning intervals
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin", function(object) standardGeneric("thin"))


#' Run the MCMC simulation.
#'
#' nStarts chains are run. b burnin iterations are run and then discarded.
#' Next, s iterations are run in each train. The user can also specify
#' an alternative number of components.
#' The mode of the MCMC simulation is also calculated.
#' @param object see showMethods(posteriorSimulation)
#' @param k The number of a priori components. This is optional and if not
#' specified, the stored k model components are used. This parameters is
#' useful for running multiple models of varying components.
#' @return An object of class 'MarginalModel' or 'BatchModel'
#' @export
#' @docType methods
#' @rdname posteriorSimulation-method
setGeneric("posteriorSimulation", function(object, k) standardGeneric("posteriorSimulation"))

setGeneric("isMarginalModel", function(object) standardGeneric("isMarginalModel"))

setGeneric("computePrior", function(object) standardGeneric("computePrior"))

setGeneric("relabel", function(object, zindex) standardGeneric("relabel"))

setGeneric("runBurnin", function(object, mcmcp) standardGeneric("runBurnin"))
setGeneric("runMcmc", function(object, mcmcp) standardGeneric("runMcmc"))

setGeneric("paramUpdates<-", function(x, value) standardGeneric("paramUpdates<-"))

#' Calculates a frequency table of latent variable assigments by observation.
#'
#' @examples
#'      zFreq(MarginalModelExample)
#' @param object see \code{showMethods(zfreq)}
#' @return An integer vector of length the number of components
#' @export
#' @docType methods
#' @rdname zfreq-method
setGeneric("zFreq", function(object) standardGeneric("zFreq"))
setGeneric("zFreq<-", function(object,value) standardGeneric("zFreq<-"))

#' Retrieve MCMC parameters from model.
#'
#' View number of iterations, burnin, etc.
#' @examples
#'      mcmcParams(MarginalModelExample)
#' @param object see \code{showMethods(mcmcParams)}
#' @return An object of class 'McmcParams'
#' @export
#' @docType methods
#' @rdname mcmcParams-method
setGeneric("mcmcParams", function(object) standardGeneric("mcmcParams"))

#' Replace MCMC parameters of model.
#'
#' Replace number of iterations, burnin, etc. Any update of the MCMC parameters will trigger an update of the chains. However, if iter (the number of MCMC iterations) is set to a nonpositive value, the chains will not be updated and kept as is.
#' @param force logical value. If false (default) the update will not proceed.
#' @param value an object of class 'McmcParams' containing the new number of iterations, etc.
#' @export
#' @docType methods
#' @rdname mcmcParams-method
setGeneric("mcmcParams<-", function(object, force=FALSE, value) standardGeneric("mcmcParams<-"))

#' Calculate log likelihood of prior for model
#'
#' @examples
#'      logPrior(MarginalModelExample)
#' @param object see \code{showMethods(logPrior)}
#' @return log likelihood of the prior.
#' @export
#' @docType methods
#' @rdname logPrior-method
setGeneric("logPrior", function(object) standardGeneric("logPrior"))

setGeneric("logPrior<-", function(object,value) standardGeneric("logPrior<-"))

setGeneric("paramUpdates", function(x) standardGeneric("paramUpdates"))

setGeneric("computePrec", function(object) standardGeneric("computePrec"))

setGeneric("marginal", function(object, batch, mcmc.params, K=1:4, maxperm=5, ...)
  standardGeneric("marginal"))

setGeneric("zChain", function(object) standardGeneric("zChain"))

setGeneric("updateMultinomialProb", function(object) standardGeneric("updateMultinomialProb"))

setGeneric("component", function(object) standardGeneric("component"))
setGeneric("overall", function(object) standardGeneric("overall"))

setGeneric("densities", function(object) standardGeneric("densities"))
setGeneric("densitiesCluster", function(object) standardGeneric("densitiesCluster"))

#' Accessor for extracting the kmeans clusters from a DensityModel
#' instance
#'
#' @param object an instance of class 'DensityModel'
#' @seealso \code{\link{DensityModel-class}}
#' @examples
#' truth <- simulateData(N=2500, p=rep(1/3, 3),
#'                       theta=c(-1, 0, 1),
#'                       sds=rep(0.1, 3))
#' dm <- DensityModel(truth)
#' clusters(dm)
#' @export
#' @docType methods
#' @rdname clusters-method
setGeneric("clusters", function(object) standardGeneric("clusters"))
setGeneric("quantiles", function(object) standardGeneric("quantiles"))

#' Compute the marginal likelihood of a converged model.
#' @examples
#'      marginalLikelihood(MarginalModelExample)
#' @param model An object of class \code{MarginalModel}, or a list of
#'        \code{MarginalModel}'s. Can also be an object of \code{BatchModel} or
#'        a list of such models.
#' @param niter The number of iterations for the reduced Gibb's sampler.
#' @return A vector of the marginal likelihood of the model(s)
#' @export
#' @docType methods
#' @rdname marginalLikelihood-method
setGeneric("marginalLikelihood",
           function(model, niter) standardGeneric("marginalLikelihood"))

#' Extract character vector of sequence names
#'
#' Short cut for \code{as.character(seqnames(g))} where g is a
#' \code{GRanges} object.
#' @param object a \code{GRanges} instance
#' @param ... currently ignored
#' @return A character vector
#' @examples
#' \dontrun{
#'    g <- GRanges("chr1", IRanges(10, 15))
#'    chromosome(g)
#' }
#' @export
setGeneric("chromosome", function(object, ...) standardGeneric("chromosome"))
