#' @include help.R
NULL

#' Number of observations
#'
#' @param model a MixtureModel-derived object
#' @examples
#' numberObs(SingleBatchModelExample)
#' @export
#' @rdname numberObs-method
setGeneric("numberObs", function(model) standardGeneric("numberObs"))

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
#' k(SingleBatchModelExample) <- 2
#' @param value An integer for the new number of components.
#' @export
#' @docType methods
#' @rdname k-method
setGeneric("k<-", function(object, value) standardGeneric("k<-"))

setGeneric("z<-", function(object, value) standardGeneric("z<-"))

setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))

setGeneric("sigma2<-", function(object, value) standardGeneric("sigma2<-"))
setGeneric("p<-", function(object, value) standardGeneric("p<-"))
setGeneric("pp<-", function(object, value) standardGeneric("pp<-"))

#' @export
#' @rdname sigma2-method
setGeneric("sigma<-", function(object, value) standardGeneric("sigma<-"))

#' Retrieve overall mean
#'
#' @examples
#'      mu(SingleBatchModelExample)
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
#' k(MultiBatchModelExample)
#' tau2(MultiBatchModelExample)
#' plot.ts(tau2(chains(MultiBatchModelExample)))
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

setGeneric("dataMean", function(object) standardGeneric("dataMean"))
setGeneric("dataPrec", function(object) standardGeneric("dataPrec"))
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
#'     theta.chain <- theta(chains(MultiBatchModelExample))
#'     dim(theta.chain)
#'     plot.ts(theta.chain, plot.type="single",
#'             col=seq_len(k(SingleBatchModelExample)))
#' @param object \code{showMethods(chains)}
#' @return The simulated chains.
#' @export
#' @docType methods
#' @rdname chains-method
setGeneric("chains", function(object) standardGeneric("chains"))

#' Accessor for Hyperparameters object for a MixtureModel-derived object
#'
#' @examples
#'     hyperParams(SingleBatchModelExample)
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
#'                         k=k(SingleBatchModelExample),
#'                         alpha=c(9, 9, 10))
#' hyperParams(SingleBatchModelExample) <- hypp
#'
#' @param value an object of class 'Hyperparameters'
#' @export
#' @docType methods
#' @rdname hyperParams-method
setGeneric("hyperParams<-", function(object,value) standardGeneric("hyperParams<-"))

setGeneric("McmcChains", function(object) standardGeneric("McmcChains"))

setGeneric("hist")




#' Retrieve batches from object.
#'
#' The batches are represented as a vector of integers.
#' @examples
#'      batch(MultiBatchModelExample)
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
#' mb <- MultiBatchModelExample
#' mcmcParams(mb) <- McmcParams(iter=100, burnin=50)
#' mb <- posteriorSimulation(mb)
#' bic(mb)
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
#' \dontrun{
#' k(SingleBatchModelExample)
#' theta(SingleBatchModelExample)
#' plot.ts(theta(chains(SingleBatchModelExample)))
#' ## BatchModel
#' k(MultiBatchModelExample)
#' length(unique(batch(MultiBatchModelExample)))
#' theta(MultiBatchModelExample)
#' ## Plot means for batches in one component
#' plot.ts(theta(chains(MultiBatchModelExample))[, 1:3])
#' }
#' @param object see \code{showMethods(theta)}
#' @return A vector of length number of components or a matrix of size 
#' number of batches x number of components
#' @export
#' @docType methods
#' @rdname theta-method
setGeneric("theta", function(object) standardGeneric("theta"))

#' @rdname theta-method
#' @export
setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))

#' Retrieve the variances of each component and batch distribution
#'
#' For a MarginalModel, this function returns a vector of variances. For a BatchModel, returns a matrix of size number of batches by number of components.
#' @examples
#'      sigma2(SingleBatchModelExample)
#' @param object see \code{showMethods(sigma2)}
#' @return A vector of length number of components or a matrix of size 
#' number of batches x number of components
#' @export
#' @docType methods
#' @rdname sigma2-method
setGeneric("sigma2", function(object) standardGeneric("sigma2"))
setGeneric("sigma_", function(object) standardGeneric("sigma_"))

#' Retrieve the probability of latent variable membership by observation.
#'
#' @examples
#'      head(probz(SingleBatchModelExample))
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
#'      nu.0(SingleBatchModelExample)
#' @param object see \code{showMethods(nu.0)}
#' @return An integer
#' @export
#' @docType methods
#' @rdname nu.0-method
setGeneric("nu.0", function(object) standardGeneric("nu.0"))

#' Retrieve the rate parameter for the sigma.2 distribution.
#'
#' @examples
#'      sigma2.0(SingleBatchModelExample)
#' @param object see \code{showMethods(sigma2.0)}
#' @return A length 1 numeric
#' @export
#' @docType methods
#' @rdname sigma2.0-method
setGeneric("sigma2.0", function(object) standardGeneric("sigma2.0"))

#' Retrieve data.
#'
#' @examples
#'      head(y(SingleBatchModelExample))
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

setGeneric("oned<-", function(object, value) standardGeneric("oned<-"))

#' Retrieve latent variable assignments.
#'
#' Retrieves the simulated latent variable assignments of each observation at each MCMC simulation.
#' @examples
#'      head(z(SingleBatchModelExample))
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
#'      modes(SingleBatchModelExample)
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
#'      eta.0(SingleBatchModelExample)
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
#'      m2.0(SingleBatchModelExample)
#' @param object see \code{showMethods(m2.0)}
#' @return m2.0 for a model
#' @export
#' @docType methods
#' @rdname m2.0-method
setGeneric("m2.0", function(object) standardGeneric("m2.0"))

setGeneric("m2.0<-", function(object,value) standardGeneric("m2.0<-"))

setGeneric("showMeans", function(object) standardGeneric("showMeans"))
setGeneric("showSigmas", function(object) standardGeneric("showSigmas"))

## Combine batches if distribution of one-dimensional summary is similar by Kolmogorov-Smirnov test statistic
##

#' Estimate batch from any sample-level surrogate variables that capture aspects of sample processing, such as the PCR experiment (e.g., the 96 well chemistry plate), laboratory, DNA source, or DNA extraction method.
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
#' CNP until no batches can be combined. For smaller values of THR, plates are more likely to be judged as similar and combined.
#'
#' @examples
#' mb.ex <- MultiBatchModelExample
#' batches <- batch(mb.ex)
#' bt <- collapseBatch(y(mb.ex), batches)
#' batches <- as.integer(factor(bt))
#' hp <- hpList(k=k(mb.ex))[["MB"]]
#' model <- MB(dat=y(mb.ex),
#'             hp=hp,
#'             batch=batches,
#'             mp=mcmcParams(mb.ex))
#' @param object see \code{showMethods(collapseBatch)}
#' @param provisional_batch a vector labelling from which batch each observation came from.
#' @param THR p-value threshold below which the null hypothesis should be rejected and batches are collapsed
#' @param nchar integer specifying the maximum number of characters in the batch labels
#' @return The new batch value.
#' @export
#' @docType methods
#' @rdname collapseBatch-method
setGeneric("collapseBatch", function(object, provisional_batch, THR=0.1, nchar=8) standardGeneric("collapseBatch"))



## #' Combine chemistry plates into batches
## #'
## #' In high-throughput assays, low-level summaries of copy number at
## #' copy number polymorphic loci (e.g., the mean log R ratio for each
## #' sample, or a principal-component derived summary) often differ
## #' between groups of samples due to technical sources of variation
## #' such as reagents, technician, or laboratory.  Technical (as opposed
## #' to biological) differences between groups of samples are referred
## #' to as batch effects.  A useful surrogate for batch is the chemistry
## #' plate on which the samples were hybridized. In large studies, a
## #' Bayesian hierarchical mixture model with plate-specific means and
## #' variances is computationally prohibitive.  However, chemistry
## #' plates processed at similar times may be qualitatively similar in
## #' terms of the distribution of the copy number summary statistic.
## #' Further, we have observed that some copy number polymorphic loci
## #' exhibit very little evidence of a batch effect, while other loci
## #' are more prone to technical variation.  We suggest combining plates
## #' that are qualitatively similar in terms of the Kolmogorov-Smirnov
## #' two-sample test of the distribution and to implement this test
## #' independently for each candidate copy number polymophism identified
## #' in a study.  The \code{combinePlates} function is a wrapper to the
## #' \code{ks.test} implemented in the \code{stats} package that
## #' compares all pairwise combinations of plates.  The \code{ks.test}
## #' is performed recursively on the batch variables defined for a given
## #' CNP until no batches can be combined.
## #' @examples
## 
## #' @param object see \code{showMethods(combinePlates)}
## #' @param plate a vector labelling from which batch each observation came from.
## #' @param THR threshold below which the null hypothesis should be rejected and batches are collapsed.
## #' @return The new batch value.
## #' @export
## #' @docType methods
## #' @rdname combinePlates-method
## setGeneric("combinePlates", function(object, plate, THR=0.1) standardGeneric("combinePlates"))

setGeneric("thetac", function(object) standardGeneric("thetac"))

setGeneric("thetaMean", function(object) standardGeneric("thetaMean"))

setGeneric("sigmaMean", function(object) standardGeneric("sigmaMean"))

setGeneric("pMean", function(object) standardGeneric("pMean"))

setGeneric("tablez", function(object) standardGeneric("tablez"))

#' Number of MCMC chains.
#'
#' This function retrieves the number of chains used for an MCMC simulation.
#' @examples
#' number_of_chains <- nStarts(SingleBatchModelExample)
#' @param object see \code{showMethods(nStarts)}
#' @return An integer of the number of different starts.
#' @export
#' @docType methods
#' @rdname nStarts-method
setGeneric("nStarts", function(object) standardGeneric("nStarts"))

#' Reset number of starting values 
#'
#' @details Simulating starting values from the priors makes it imperative to
#'   run a large nubmer of simulations for burnin and to carefully evaluate the
#'   chains following burning for convergence. The adequacy of the burnin is
#'   difficult to assess in high-dimensional settings with a large number of
#'   CNPs. To avoid starting in regions of low posterior probabilitiy, we use
#'   existing EM-based methods in the package \code{{mclust}} to select starting
#'   values from \code{N} bootstrap sample of the observed data, where \code{N}
#'   is specificed as in the example below. For each bootstrap sample, starting
#'   values for the model are estimated. For each set of simulated starting
#'   values, the log likelihood of the full data is evaluated. The starting
#'   values with the largest log likelihood are used as initial values for the
#'   MCMC simulations.
#'
#' @examples
#' number_of_chains <- 10
#' nStarts(SingleBatchModelExample) <- number_of_chains
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
#' head(log_lik(chains(SingleBatchModelExample)))
#' ## retrieve log likelihood at last MCMC iteration
#' log_lik(SingleBatchModelExample)
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
#' burnin(SingleBatchModelExample)
#' mp <- mcmcParams(SingleBatchModelExample)
#' burnin(mp)
#' @param object see \code{showMethods(burnin)}
#' @return The number of burnin simulations.
#' @export
#' @docType methods
#' @rdname burnin-method
setGeneric("burnin", function(object) standardGeneric("burnin"))

setGeneric("min_GR", function(object) standardGeneric("min_GR"))
setGeneric("min_effsize", function(object) standardGeneric("min_effsize"))
setGeneric("max_burnin", function(object) standardGeneric("max_burnin"))
setGeneric("min_chains", function(object) standardGeneric("min_chains"))

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
setGeneric("iter<-", function(object, value) standardGeneric("iter<-"))

#' Number of MCMC iterations.
#'
#' This function retrieves the number of iterations of an MCMC simulation.
#' @examples
#'      iter(SingleBatchModelExample)
#' @param object see \code{showMethods(iter)}
#' @return The number of MCMC iterations
#' @export
#' @docType methods
#' @rdname iter-method
setGeneric("iter", function(object) standardGeneric("iter"))

#' Get or set the number of thinning intervals.
#'
#' This function gets or sets the number of thinning intervals used for an MCMC
#' simulation.
#' @examples
#'      thin(SingleBatchModelExample)
#' @param object see showMethods(thin)
#' @return An integer of the number of thinning intervals
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin", function(object) standardGeneric("thin"))

#' @examples
#'      thin(SingleBatchModelExample) <- 10L
#' @export
#' @docType methods
#' @rdname thin-method
setGeneric("thin<-", function(object, value) standardGeneric("thin<-"))

#' Run MCMC simulation.
#'
#' nStarts chains are run. b burnin iterations are run and then discarded.
#' Next, s iterations are run in each train. The user can also specify
#' an alternative number of components.
#' The mode of the MCMC simulation is also calculated.
#' @examples
#' # Fit model with pre-specified number of components (k=3)
#' set.seed(123)
#' ## specify small number of iterations so that the example runs quickly
#' mp <- McmcParams(iter=2, burnin=0, nStarts=3)
#' sb <- SingleBatchModelExample
#' mcmcParams(sb) <- mp
#' posteriorSimulation(sb)
#'
#' # Run additional iterations, but set nStart = 0 so that the last value of the
#' # chain is the first value of the next chain
#' mcmcParams(sb) <- McmcParams(iter=5, nStarts=0, burnin=0)
#' posteriorSimulation(sb)
#'
#' # Fit batch models of different sizes (k=1 and 2)
#' mb <- MultiBatchModelExample
#' mcmcParams(mb) <- mp
#' yy <- sample(y(mb), 300)
#' batches <- rep(1:3, length.out=length(yy))
#' mp <- McmcParams(iter=1000, burnin=500, thin=1, nStarts=4)
#' \dontrun{
#'   mlist <- gibbs(model="MB", k_range=c(1, 2), dat=yy, batches=batches)
#' }
#' @param object see showMethods(posteriorSimulation)
#' @param k The number of a priori components. This is optional and if not
#' specified, the stored k model components are used. This parameters is
#' useful for running multiple models of varying components.
#' @return An object of class 'MarginalModel' or 'BatchModel'
#' @seealso  \code{\link{ggChains}} for diagnosing convergence.  See \code{\link{ggMixture}} for plotting the model-based densities.
#' @export
#' @docType methods
#' @rdname posteriorSimulation-method
setGeneric("posteriorSimulation", function(object, k) standardGeneric("posteriorSimulation"))

## setGeneric("posteriorSimulation2", function(object, params)  {
##   standardGeneric("posteriorSimulation2")
## })

setGeneric("isSB", function(object) standardGeneric("isSB"))

setGeneric("computePrior", function(object) standardGeneric("computePrior"))

setGeneric("relabel", function(object, zindex) standardGeneric("relabel"))

setGeneric("runBurnin", function(object, mcmcp) standardGeneric("runBurnin"))
setGeneric("runMcmc", function(object, mcmcp) standardGeneric("runMcmc"))

setGeneric("paramUpdates<-", function(x, value) standardGeneric("paramUpdates<-"))

#' Calculates a frequency table of latent variable assigments by observation.
#'
#' @examples
#'      zFreq(SingleBatchModelExample)
#' @param object see \code{showMethods(zfreq)}
#' @return An integer vector of length the number of components
#' @export
#' @docType methods
#' @rdname zfreq-method
setGeneric("zFreq", function(object) standardGeneric("zFreq"))
setGeneric("zFreq<-", function(object,value) standardGeneric("zFreq<-"))

#' Calculates a frequency table of latent variable assigments for parents by observation.
#' 
#' @examples
#'      zfreqpar(TrioBatchModelExample)
#' @param object see \code{showMethods(zfreqpar)}
#' @return An integer vector of length the number of components
#' @export
#' @docType methods
#' @rdname zfreqpar-method
setGeneric("zFreqPar", function(object) standardGeneric("zFreqPar"))
setGeneric("zFreqPar<-", function(object,value) standardGeneric("zFreqPar<-"))

#' Retrieves intensity data from trios
#' 
#' @examples
#'      triodata_lrr(TrioBatchModelExample)
#' @param object see \code{showMethods(triodata_lrr)}
#' @return An integer vector of length the number of components
#' @export
#' @docType methods
#' @rdname triodata_lrr-method
setGeneric("triodata_lrr", function(object) standardGeneric("triodata_lrr"))

#' Retrieve MCMC parameters from model.
#'
#' View number of iterations, burnin, etc.
#' @examples
#'      mcmcParams(SingleBatchModelExample)
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
setGeneric("mcmcParams<-", function(object, value) standardGeneric("mcmcParams<-"))

#' Calculate log likelihood of prior for model
#'
#' @examples
#'      logPrior(SingleBatchModelExample)
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

setGeneric("quantiles", function(object) standardGeneric("quantiles"))

#' Compute the marginal likelihood of a converged model.
#'
#' The recommended function for fitting mixture models and evaluating convergence is through the `gibbs` function. This function will return a list of models ordered by the marginal likelihood. The marginal likelihood is computed using the Chib's estimator (JASA, Volume 90 (435), 1995).
#'
#'
#' @examples
#' ## In practice, run a much longer burnin and increase the number of
#' ## iterations to save after burnin
#'    mm <- SingleBatchModelExample
#'    mcmcParams(mm) <- McmcParams(iter=50, burnin=0, nStarts=0)
#'    mm <- posteriorSimulation(mm)
#'    marginalLikelihood(mm)
#' @param model An object of class \code{MarginalModel}, or a list of
#'        \code{MarginalModel}'s. Can also be an object of \code{BatchModel} or
#'        a list of such models.
#' @param params A list containing parameters for marginalLikelihood computation. See \code{mlParams} for details.
#'
#' @seealso See \code{\link{mlParams}} for parameters related to computing the log marginal likelihood via Chib's estimator. See \code{\link{gibbs}} for fitting multiple mixture models and returning a list sorted by the marginal likelihood.  See \code{\link{marginal_lik}} for the accessor.
#'
#' Note: currently thinning of the reduced MCMC chains is not allowed.
#'
#' @return A vector of the marginal likelihood of the model(s)
#' @export
#' @docType methods
#' @rdname marginalLikelihood-method
setGeneric("marginalLikelihood", function(model, params=mlParams()) standardGeneric("marginalLikelihood"))

#' Extract character vector of sequence names
#'
#' Short cut for \code{as.character(seqnames(g))} where g is a
#' \code{GRanges} object.
#' @param object a \code{GRanges} instance
#' @param ... currently ignored
#' @return A character vector
#' @examples
#'    g <- GRanges("chr1", IRanges(10, 15))
#'    chromosome(g)
#' @rdname chromosome
#' @export
setGeneric("chromosome", function(object, ...) standardGeneric("chromosome"))

#' @rdname chromosome
#' @aliases chromosome,GenomicRanges-method
setMethod("chromosome", "GenomicRanges", function(object, ...){
  as.character(seqnames(object))
})



#' Accessor for determing whether label switching occurred during MCMC
#'
#'
#' @param object MixtureModel-derived class
#' @export
#' @examples
#' label_switch(SingleBatchModelExample)
setGeneric("label_switch", function(object) standardGeneric("label_switch"))


setGeneric("label_switch<-", function(object, value) standardGeneric("label_switch<-"))

#' Accessor for the log marginal likelihood of a SB, SBP, MB, or MBP model
#'
#' The marginal likelihood is computed by Chib's estimator (JASA, Volume 90 (435), 1995).
#'
#' @seealso See \code{\link{marginalLikelihood}} for computing the marginal likelihood of a mixture model.
#' @param object a SB, SBP, MB, or MBP model
#' @export
#' @examples
#' sb <- SingleBatchModelExample
#' marginal_lik(sb)
setGeneric("marginal_lik", function(object) standardGeneric("marginal_lik"))

#' @export
#' @rdname marginal_lik
setGeneric("marginal_lik<-", function(object, value) standardGeneric("marginal_lik<-"))

setGeneric("updateZ", function(object) standardGeneric("updateZ"))

setGeneric("gatherChains", function(object) standardGeneric("gatherChains"))

setGeneric("sortComponentLabels", function(model) standardGeneric("sortComponentLabels"))

setGeneric("isOrdered", function(object) standardGeneric("isOrdered"))

#' Trace plots of MCMC chains and mixture model densities
#'
#' The \code{ggChains} method provides a convenient wrapper for plotting the chains of all parameters in the various mixture model implementations.  In addition to the estimated number of independent MCMC draws (effective sample size) and Gelman-Rubin convergence diagnostics implemented in \code{gibbs}, visualization of the chains is helpful for assessing convergence.
#'
#' The \code{ggMixture} method overlays the density of the posterior predictive distribution of the Gaussian mixture on the empirical data. \code{ggMixture} assumes that you have already run the Gibbs sampler either by the \code{gibbs} function or by the \code{posteriorSimulation} function.
#'
#' @seealso \code{\link{gibbs}} 
#'
#' @param model A SB, MB, SBP, or MBP model
#' @rdname ggplot-functions
#' @return A \code{gg} object
#' @export
#' @examples
#'   sb <- SingleBatchModelExample
#'   iter(sb) <- 1000
#'   burnin(sb) <- 100
#'   sb <- posteriorSimulation(sb)
#'   fig.chains <- ggChains(sb)
#'   ## component-specific chains
#'   fig.chains[["comp"]]
#'   ## single-parameter chains and log-likelihood
#'   fig.chains[["single"]]
#'
#'   ## plot the mixture
#'   fig.mix <- ggMixture(sb)
setGeneric("ggChains", function(model) standardGeneric("ggChains"))

#' @rdname ggplot-functions
#' @param bins a length-one numeric vector indicating the number of bins -- passed to \code{geom_hist}
#' @export
#' @rdname ggplot-functions
setGeneric("ggMixture", function(model, bins=100) standardGeneric("ggMixture"))

## #' @export
## #' @rdname ggplot-functions
## setGeneric("ggMultiBatch", function(model, bins=100) standardGeneric("ggMultiBatch"))

## #' @param model a SB, MB, SBP, or MBP model
## #' @examples
## #'
## #' @export
## #' @return a SB, MB, SBP, or MBP model
## #' @rdname tile-functions
## setGeneric("upSample", function(model, tiles) standardGeneric("upSample"))

#' Constructs a CopyNumberModel from SB, SBP, MB, or MBP models
#'
#' The mixture components do not necessarily reflect distinct copy number
#' states, possibly due to skewed (non-Gaussian) log R ratios. While easy to fit skewed data with a finite mixture of Gaussians, additional steps are needed to assess whether the components correspond to distinct copy number states.  This accessor \code{copyNumber} returns the copy number states -- i.e., the result after mapping mixture components to copy number states.
#'
#' @param model a SB, SBP, MB, or MBP model
#' @param params a list of parameters used for mapping mixture components to copy number states.
#' @seealso \code{\link{copyNumber}}
#' @export
#' @examples
#' sb <- SingleBatchModelExample
#' cn.model <- CopyNumberModel(sb, mapParams())
#' @rdname CopyNumber-methods 
setGeneric("CopyNumberModel", function(model, params=mapParams()) standardGeneric("CopyNumberModel"))

#' Map mixture components to copy number states
#'
#'
#' @export
#' @examples
#' cn.model <- CopyNumberModel(SingleBatchModelExample)
#' ## manually remap first two components to the same copy number state
#' \dontrun{
#'  mapping(cn.model) <- c(1, 1, 2)
#'  ggMixture(cn.model)
#' }
#' @param object a SB, SBP, MB, or MBP model
#' @seealso \code{\link{CopyNumberModel}}
#' @rdname mapping
setGeneric("mapping", function(object) standardGeneric("mapping"))

#' @param value a k-length numeric vector with values in {1, 2, ..., k}, where k is the number of mixture components
#' @export
#' @rdname mapping
setGeneric("mapping<-", function(object, value) standardGeneric("mapping<-"))

setGeneric("numberStates", function(model) standardGeneric("numberStates"))

#' Posterior probabilities for copy number states
#'
#' In contrast to posterior probabilities for mixture components, this function
#' returns posterior probabilities for distinct copy number states.
#' a \code{SingleBatchCopyNumber} or \code{MultiBatchCopyNumber} instance
#' @param model a SB, SBP, MB, or MBP model
#' @rdname probCopyNumber
#' @seealso \code{\link{CopyNumberModel}}
#' @export
setGeneric("probCopyNumber", function(model) standardGeneric("probCopyNumber"))

#' Extract copy number estimates from a `CopyNumberModel`
#'
#'
#' @param object a \code{SingleBatchCopyNumber} or \code{MultiBatchCopyNumber} object
#' @examples
#' sb <- SingleBatchModelExample
#' cn.model <- CopyNumberModel(sb)
#' head(copyNumber(cn.model))
#'
#' ## here is an identity mapping
#' \dontrun{
#'     mapping(cn.model) <- 1:3
#'     identical(copyNumber(cn.model), z(cn.model))
#'     table(copyNumber(cn.model))
#'
#'     ## here, we map the first two mixture components to one copy number state
#'     mapping(cn.model) <- c(1, 1, 2)
#'     table(copyNumber(cn.model))
#' }
#' @seealso \code{\link{CopyNumberModel}}
#' @export
#' @rdname copyNumber
setGeneric("copyNumber", function(object) standardGeneric("copyNumber"))

setGeneric("mergeComponents", function(model, j) standardGeneric("mergeComponents"))

### t-distribution stuff

#' Accessor for degrees of freedom
#'
#' @param object a Hyperparameters- or MixtureModel-derived class
#' @export
#' @examples
#'   hp <- Hyperparameters()
#'   dfr(hp)
#'   dfr(hp) <- 10
#'   dfr(hp)
#' @rdname dfr-method
setGeneric("dfr", function(object) standardGeneric("dfr"))

#' @export
#' @rdname dfr-method
setGeneric("dfr<-", function(object, value) standardGeneric("dfr<-"))

setGeneric("u", function(object) standardGeneric("u"))


#' Retrieve the probability of latent variable membership by observation for parents.
#'
#' @examples
#'      probzpar(TrioBatchModelExample)
#' @param object see \code{showMethods(probzpar)}
#' @return A matrix of size number of observations x number of components
#' @export
#' @docType methods
#' @rdname probzpar-method
setGeneric("probzpar", function(object) standardGeneric("probzpar"))

setGeneric("probzpar<-", function(object, value) standardGeneric("probzpar<-"))

setGeneric("p", function(object) standardGeneric("p"))

setGeneric("sigma_<-", function(object, value) standardGeneric("sigma_<-"))

setGeneric("flags", function(object) standardGeneric("flags"))
setGeneric("flags<-", function(object, value) standardGeneric("flags<-"))
setGeneric("current_values<-", function(object, value) standardGeneric("current_values<-"))
setGeneric("current_values", function(object, value) standardGeneric("current_values"))
setGeneric("u<-", function(object, value) standardGeneric("u<-"))
setGeneric("parameters", function(object) standardGeneric("parameters"))
setGeneric("parameters<-", function(object, value) standardGeneric("parameters<-"))
setGeneric("summaries", function(object) standardGeneric("summaries"))
setGeneric("summaries<-", function(object, value) standardGeneric("summaries<-"))
setGeneric("dataSd", function(object) standardGeneric("dataSd"))

#' @export
setGeneric("findSurrogates", function(object, THR=0.1) standardGeneric("findSurrogates"))
setGeneric("isSampled", function(object, THR=0.1) standardGeneric("isSampled"))
setGeneric("down_sample", function(object) standardGeneric("down_sample"))
setGeneric("down_sample<-", function(object, value) standardGeneric("down_sample<-"))
setGeneric("downSampledData", function(object) standardGeneric("downSampledData"))
setGeneric("downSampledData<-", function(x, value) standardGeneric("downSampledData<-"))
setGeneric("probability_z", function(object) standardGeneric("probability_z"))
setGeneric("mcmc2", function(object, guide) standardGeneric("mcmc2"))
setGeneric("setModes", function(object) standardGeneric("setModes"))
setGeneric("upsample_z", function(object) standardGeneric("upsample_z"))
setGeneric("specs", function(object) standardGeneric("specs"))
setGeneric("specs<-", function(object, value) standardGeneric("specs<-"))
setGeneric("downSampleModel", function(object, N, i) standardGeneric("downSampleModel") )
setGeneric("upSampleModel", function(object) standardGeneric("upSampleModel"))
setGeneric("listChains", function(object) standardGeneric("listChains"))
setGeneric("predictive", function(object) standardGeneric("predictive"))
setGeneric("zstar", function(object) standardGeneric("zstar"))
setGeneric("predictive<-", function(object, value) standardGeneric("predictive<-"))
setGeneric("revertBack", function(object, mbm) standardGeneric("revertBack"))
setGeneric("numBatch", function(object) standardGeneric("numBatch"))
setGeneric("numBatch<-", function(object, value) standardGeneric("numBatch<-"))
setGeneric("max_burnin<-", function(object, value) standardGeneric("max_burnin<-"))
setGeneric("useModes", function(object) standardGeneric("useModes"))
setGeneric("compute_marginal_lik", function(object, params) standardGeneric("compute_marginal_lik"))
setGeneric("modelName", function(object) standardGeneric("modelName"))
setGeneric("singleBatchGuided", function(x, guide) standardGeneric("singleBatchGuided"))

setGeneric("convergence", function(object) standardGeneric("convergence"))
