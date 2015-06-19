#' @include help.R
NULL

#' Number of components.
#'
#' This function retrieves the number of a priori components.
#' @param object see \code{showMethods(k)}
#' @export
#' @docType methods
#' @rdname k-method
setGeneric("k", function(object) standardGeneric("k"))

setGeneric("k<-", function(object, value) standardGeneric("k<-"))

setGeneric("z<-", function(object, value) standardGeneric("z<-"))

setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))
setGeneric("sigma2<-", function(object, value) standardGeneric("sigma2<-"))
setGeneric("p<-", function(object, value) standardGeneric("p<-"))

#' Retrieve overall mean
#'
#' @param object see \code{showMethods(mu)}
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
#' data(BatchModelExample)
#' k(BatchModelExample)
#' tau2(BatchModelExample)
#' plot.ts(tau2(chains(BatchModelExample)))
#'
#' @param object see \code{showMethods(tau2)}
#' @export
#' @docType methods
#' @rdname tau2-method
#' @seealso \code{Hyperparameters}
setGeneric("tau2", function(object) standardGeneric("tau2"))
setGeneric("tau2<-", function(object, value) standardGeneric("tau2<-"))
setGeneric("nu.0<-", function(object, value) standardGeneric("nu.0<-"))
setGeneric("sigma2.0<-", function(object, value) standardGeneric("sigma2.0<-"))
setGeneric("logpotential<-", function(object, value) standardGeneric("logpotential<-"))
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
#' data(MarginalModelExample)
#' theta.chain <- theta(chains(MarginalModelExample))
#' dim(theta.chain)
#' plot.ts(theta.chain, plot.type="single",
#'         col=seq_len(k(MarginalModelExample)))
#' @param object \code{showMethods(chains)}
#' @export
#' @docType methods
#' @rdname chains-method
setGeneric("chains", function(object) standardGeneric("chains"))

#' Accessor for Hyperparameters object for a MixtureModel-derived object
#'
#' @examples
#' \dontrun{
#'     data(MarginalModelExample)
#'     hyperParams(MixtureModelExample)
#' }
#' @param object see \code{showMethods(hyperParams)}
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams", function(object) standardGeneric("hyperParams"))

#' Replace the hyperparameters for a \code{MixtureModel}-derived object
#'
#' @examples
#' data(MarginalModelExample)
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
#' @param ... Additional arguments passed to \code{hist}.
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
#' \dontrun{
#'      data(BatchModelExample)
#'      batch(BatchModelExample)
#' }
#' @param object see \code{showMethods(batch)}
#' @export
#' @docType methods
#' @rdname batch-method
setGeneric("batch", function(object) standardGeneric("batch"))
setGeneric("batch<-", function(object,value) standardGeneric("batch<-"))

setGeneric("startingValues", function(object, params, zz) standardGeneric("startingValues"))

setGeneric("computeMeans", function(object) standardGeneric("computeMeans"))
setGeneric("computeVars", function(object) standardGeneric("computeVars"))

setGeneric("updateZ", function(object) standardGeneric("updateZ"))

setGeneric("computePotential", function(object) standardGeneric("computePotential"))

setGeneric("dat", function(object) standardGeneric("dat"))
setGeneric("dat<-", function(object,value) standardGeneric("dat<-"))

setGeneric("initializeSigma2.0", function(object) standardGeneric("initializeSigma2.0"))

setGeneric("initializeSigma2", function(object) standardGeneric("initializeSigma2"))

setGeneric("initializeMu", function(object) standardGeneric("initializeMu"))

setGeneric("initializeTheta", function(object) standardGeneric("initializeTheta"))

setGeneric("initializeTau2", function(object) standardGeneric("initializeTau2"))

setGeneric("posteriorMultinomial", function(object) standardGeneric("posteriorMultinomial"))

setGeneric("simulateY", function(object, N) standardGeneric("simulateY"))
setGeneric("batchCorrect", function(object) standardGeneric("batchCorrect"))

setGeneric("moveChain", function(object, s) standardGeneric("moveChain"))

setGeneric("updateThetaCpp", function(object, constrain) standardGeneric("updateThetaCpp"))
setGeneric("updateSigma2Cpp", function(object) standardGeneric("updateSigma2Cpp"))
setGeneric("updateTheta", function(object, constrain) standardGeneric("updateTheta"))
setGeneric("updateSigma2", function(object) standardGeneric("updateSigma2"))
setGeneric("updateSigma2.0", function(object) standardGeneric("updateSigma2.0"))

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))

#' Calculate BIC of a model
#'
#' @param object see \code{showMethods(bic)}
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
#' data(MarginalModelExample)
#' k(MarginalModelExample)
#' theta(MarginalModelExample)
#' plot.ts(theta(chains(MarginalModelExample)))
#' ## BatchModel
#' data(BatchModelExample)
#' k(BatchModelExample)
#' length(unique(batch(BatchModelExample)))
#' theta(BatchModelExample)
#' ## Plot means for batches in one component
#' plot.ts(theta(chains(BatchModelExample))[, 1:3])
#' @param object see \code{showMethods(theta)}
#' @export
#' @docType methods
#' @rdname theta-method
setGeneric("theta", function(object) standardGeneric("theta"))

#' Retrieve the variances of each component and batch distribution
#'
#' For a MarginalModel, this function returns a vector of variances. For a BatchModel, returns a matrix of size number of batches by number of components.
#' @param object see \code{showMethods(sigma2)}
#' @export
#' @docType methods
#' @rdname sigma2-method
setGeneric("sigma2", function(object) standardGeneric("sigma2"))

setGeneric("reorderComponents", function(object, new_levels) standardGeneric("reorderComponents"))

setGeneric("hwe", function(object) standardGeneric("hwe"))

#' Retrieve the probability of latent variable membership by observation.
#'
#' @param object see \code{showMethods(probz)}
#' @export
#' @docType methods
#' @rdname probz-method
setGeneric("probz", function(object) standardGeneric("probz"))

setGeneric("probz<-", function(object, value) standardGeneric("probz<-"))

setGeneric("fitMixtureModels", function(object, mcmcp, K=1:5, batch) standardGeneric("fitMixtureModels"))

#' Retrieve the shape parameter for the sigma.2 distribution.
#'
#' @param object see \code{showMethods(nu.0)}
#' @export
#' @docType methods
#' @rdname nu.0-method
setGeneric("nu.0", function(object) standardGeneric("nu.0"))

#' Retrieve the rate parameter for the sigma.2 distribution.
#'
#' @param object see \code{showMethods(sigma2.0)}
#' @export
#' @docType methods
#' @rdname sigma2.0-method
setGeneric("sigma2.0", function(object) standardGeneric("sigma2.0"))

#' Retrieve data.
#'
#' @param object see \code{showMethods(y)}
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
#' @export
#' @docType methods
#' @rdname oned-method
#' @export
setGeneric("oned", function(object) standardGeneric("oned"))

#' Retrieve latent variable assignments.
#'
#' Retrieves the simulated latent variable assignments of each observation at each MCMC simulation.
#' @param object see \code{showMethods(z)}
#' @export
#' @docType methods
#' @rdname z-method
setGeneric("z", function(object) standardGeneric("z"))

setGeneric("modalParameters", function(object) standardGeneric("modalParameters"))

setGeneric("computeModes", function(object) standardGeneric("computeModes"))

setGeneric("switchLabels", function(object) standardGeneric("switchLabels"))

setGeneric("computeDistance", function(object) standardGeneric("computeDistance"))

#' Retrieve the modes from a model.
#'
#' The iteration which maximizes log likelihood and log prior is found. The estimates for each parameter at this iteration are retrieved.
#' @param object a \code{MixtureModel}-derived class
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
#' @param object see \code{showMethods(eta.0)}
#' @export
#' @docType methods
#' @rdname eta.0-method
setGeneric("eta.0", function(object) standardGeneric("eta.0"))
setGeneric("eta.0<-", function(object,value) standardGeneric("eta.0<-"))

#' Retrieve the shape parameter for the tau2 distribution.
#'
#' @param object see \code{showMethods(m2.0)}
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
#' data(BatchModelExample)
#' bt <- collapseBatch(y(BatchModelExample), batch(BatchModelExample))
#' newBatchModel <- BatchModel(y(BatchModelExample), k(BatchModelExample),
#'                             bt, Hyperparameters(BatchModelExample),
#'                             mcmcParams(BatchModelExample))
#' @param object see \code{showMethods(collapseBatch)}
#' @param plate a vector labelling from which batch each observation came from.
#' @param THR threshold below which the null hypothesis should be rejected and batches are collapsed.
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
#' data(BatchModelExample)
#' tracePlot(BatchModelExample, "theta")
#' tracePlot(BatchModelExample, "sigma")
#' @param object see \code{showMethods(tracePlot)}
#' @param name the name of the parameter for which to plot values. Can be 'theta', 'sigma', 'p', 'mu', or 'tau'.
#' @param ... Other argument to pass to plot.
#' @export
#' @docType methods
#' @rdname tracePlot-method
setGeneric("tracePlot", function(object, name, ...) standardGeneric("tracePlot"))

setGeneric("tablez", function(object) standardGeneric("tablez"))

#' Number of MCMC chains.
#'
#' This function retrieves the number of chains used for an MCMC simulation.
#' @examples
#' data(MarginalModelExample)
#' number_of_chains <- nStarts(MarginalModelExample)
#' @param object see \code{showMethods(nStarts)}
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts", function(object) standardGeneric("nStarts"))

#' Reset number of MCMC chains in simulation.
#'
#' This function changes the number of chains used for an MCMC simulation.
#' @examples
#' data(MarginalModelExample)
#' number_of_chains <- 3
#' nStarts(MarginalModelExample) <- number_of_chains
#' @param value new number of chains
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts<-", function(object, value) standardGeneric("nStarts<-"))

setGeneric("alpha", function(object) standardGeneric("alpha"))

setGeneric("orderTheta<-", function(object, value) standardGeneric("orderTheta<-"))
setGeneric("orderTheta", function(object) standardGeneric("orderTheta"))

#' Retrieve log likelihood.
#'
#' @examples
#' data(MarginalModelExample)
#' ## retrieve log likelihood at each MCMC iteration
#' log_lik(chains(MarginalModelExample))
#' ## retrieve log likelihood at last MCMC iteration
#' log_lik(MarginalModelExample)
#' @param object see showMethods(log_lik)
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
#' data(MarginalModelExample)
#' burnin(MarginalModelExample)
#' mp <- mcmcParams(MarginalModelExample)
#' burnin(mp)
#' @param object see \code{showMethods(burnin)}
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
#' @param object see \code{showMethods(iter)}
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter", function(object) standardGeneric("iter"))

#' Number of thinning intervals.
#'
#' This function retrieves the number of thinning intervals used for an MCMC simulation.
#' @param object see showMethods(thin)
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin", function(object) standardGeneric("thin"))


#' Run the MCMC simulation.
#'
#' nStarts chains are run. b burnin iterations are run and then discarded. Next, s iterations are run in each train. The mode of the MCMC simulation is also calculated.
#' @param object see showMethods(posteriorSimulation)
#' @export
#' @docType methods
#' @rdname posteriorSimulation-method
setGeneric("posteriorSimulation", function(object) standardGeneric("posteriorSimulation"))

setGeneric("initializeModel", function(params, hypp) standardGeneric("initializeModel"))

setGeneric("posteriorTheta", function(object, mcmcp) standardGeneric("posteriorTheta"))
setGeneric("posteriorSigma2", function(object, mcmcp) standardGeneric("posteriorSigma2"))

setGeneric("post1", function(object) standardGeneric("post1"))
setGeneric("post2", function(object) standardGeneric("post2"))
setGeneric("post3", function(object) standardGeneric("post3"))
setGeneric("isMarginalModel", function(object) standardGeneric("isMarginalModel"))

setGeneric("computePrior", function(object) standardGeneric("computePrior"))
setGeneric("computeLogLikxPrior", function(object) standardGeneric("computeLogLikxPrior"))

setGeneric("postFiles", function(object) standardGeneric("postFiles"))

setGeneric("relabel", function(object, zindex) standardGeneric("relabel"))

setGeneric("HyperParameterList", function(hypp, K) standardGeneric("HyperParameterList"))

#' Creates a list of ModelParams for differing number of components
#'
#' Simplifies creating multiple objects for model parameters for multiple component sizes. Used when number of components is not known. Can be used with Marginal and Batch models.
#' @param hypp An object of class 'Hyperparameters'
#' @param K An integer vector containing the different component sizes for which you wish to create ModelParams.
#' @param data A vector containing the data.
#' @param mcmcp An object of class 'ModelParams'. Can have varying numbers of iterations, etc for each component size.
#' @param batch A vector indicating the batches from which each observation came from.
#' @export
#' @docType methods
#' @rdname ModelParamList-method
setGeneric("ModelParamList", function(hypp, K, data, mcmcp, batch) standardGeneric("ModelParamList"))

setGeneric("runBurnin", function(object, mcmcp) standardGeneric("runBurnin"))
setGeneric("runMcmc", function(object, mcmcp) standardGeneric("runMcmc"))

setGeneric("paramUpdates<-", function(x, value) standardGeneric("paramUpdates<-"))

#' Calculates a frequency table of latent variable assigments by observation.
#'
#' @param object see \code{showMethods(zfreq)}
#' @export
#' @docType methods
#' @rdname zfreq-method
setGeneric("zFreq", function(object) standardGeneric("zFreq"))
setGeneric("zFreq<-", function(object,value) standardGeneric("zFreq<-"))

#' Retrieve MCMC parameters from model.
#'
#' View number of iterations, burnin, etc.
#' @param object see \code{showMethods(mcmcParams)}
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
#' @param object see \code{showMethods(logPrior)}
#' @export
#' @docType methods
#' @rdname logPrior-method
setGeneric("logPrior", function(object) standardGeneric("logPrior"))

setGeneric("logPrior<-", function(object,value) standardGeneric("logPrior<-"))

setGeneric("fullGibbs", function(object, mcmcp) standardGeneric("fullGibbs"))

setGeneric("paramUpdates", function(x) standardGeneric("paramUpdates"))

setGeneric("computePrec", function(object) standardGeneric("computePrec"))

setGeneric("marginal", function(object, batch, mcmc.params, K=1:4, maxperm=5, ...)
  standardGeneric("marginal"))

setGeneric("modelList", function(object) standardGeneric("modelList"))

setGeneric("modelList<-", function(object, value) standardGeneric("modelList<-"))


setGeneric("computeMarginalEachK", function(object, K=1:4, hypp, mcmcp=McmcParams(), maxperm=5)
  standardGeneric("computeMarginalEachK"))

setGeneric("computeMarginalEachK2",
           function(object, batch, maxperm=3,
                    K=1:4, mcmcp=McmcParams(),
                    hypp) standardGeneric("computeMarginalEachK2"))

setGeneric("best", function(object) standardGeneric("best"))

setGeneric("posteriorP", function(object) standardGeneric("posteriorP"))

setGeneric("zChain", function(object) standardGeneric("zChain"))
setGeneric("zChain<-", function(object,value) standardGeneric("zChain<-"))

setGeneric("permuteModes", function(object, ix) standardGeneric("permuteModes"))

setGeneric("pThetaStar", function(kmod, maxperm=5, T2) standardGeneric("pThetaStar"))

setGeneric("pSigma2", function(object) standardGeneric("pSigma2"))

setGeneric("pMixProb", function(object) standardGeneric("pMixProb"))

setGeneric("reducedGibbsThetaFixed", function(object)
  standardGeneric("reducedGibbsThetaFixed"))

setGeneric("reducedGibbsThetaSigmaFixed", function(object)
  standardGeneric("reducedGibbsThetaSigmaFixed"))

setGeneric("reducedGibbsZThetaFixed", function(object) standardGeneric("reducedGibbsZThetaFixed"))

setGeneric("pTheta_Zfixed", function(object) standardGeneric("pTheta_Zfixed"))

setGeneric("pTheta", function(object) standardGeneric("pTheta"))

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
