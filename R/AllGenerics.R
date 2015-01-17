#' @include help.R

##setGeneric("rv", function(object) standardGeneric("rv"))

#' @export
setGeneric("k", function(object) standardGeneric("k"))
##setGeneric("params", function(object) standardGeneric("params"))

##setGeneric("posteriorMean", function(object) standardGeneric("posteriorMean"))
##setGeneric("posteriorPrec", function(object) standardGeneric("posteriorPrec"))
##
##setGeneric("posteriorMean<-", function(object, value) standardGeneric("posteriorMean<-"))
##setGeneric("posteriorPrec<-", function(object, value) standardGeneric("posteriorPrec<-"))

setGeneric("z<-", function(object, value) standardGeneric("z<-"))

setGeneric("theta<-", function(object, value) standardGeneric("theta<-"))
setGeneric("sigma2<-", function(object, value) standardGeneric("sigma2<-"))
setGeneric("p<-", function(object, value) standardGeneric("p<-"))
setGeneric("mu", function(object) standardGeneric("mu"))
setGeneric("mu<-", function(object, value) standardGeneric("mu<-"))
setGeneric("tau2", function(object) standardGeneric("tau2"))
setGeneric("tau2<-", function(object, value) standardGeneric("tau2<-"))
setGeneric("nu.0<-", function(object, value) standardGeneric("nu.0<-"))
setGeneric("sigma2.0<-", function(object, value) standardGeneric("sigma2.0<-"))
setGeneric("logpotential<-", function(object, value) standardGeneric("logpotential<-"))
setGeneric("dataMean<-", function(object, value) standardGeneric("dataMean<-"))
setGeneric("dataPrec<-", function(object, value) standardGeneric("dataPrec<-"))

#' @export
setGeneric("mcmcChains<-", function(object, value) standardGeneric("mcmcChains<-"))

#' @export
setGeneric("mcmcChains", function(object) standardGeneric("mcmcChains"))

#' @export
setGeneric("hyperParams", function(object) standardGeneric("hyperParams"))

#' @export
setGeneric("hyperParams<-", function(object,value) standardGeneric("hyperParams<-"))

#' @export
setGeneric("McmcChains", function(object, mcmc.params) standardGeneric("McmcChains"))

setGeneric("hist")
setGeneric("plot")

##setGeneric("BatchModel", function(object, ncomp, batch) standardGeneric("BatchModel"))


##setGeneric("batch", function(object) standardGeneric("batch"))

setGeneric("batch<-", function(object,value) standardGeneric("batch<-"))

setGeneric("startingValues", function(object) standardGeneric("startingValues"))

setGeneric("computeMeans", function(object) standardGeneric("computeMeans"))
setGeneric("computeVars", function(object) standardGeneric("computeVars"))

setGeneric("updateZ", function(object) standardGeneric("updateZ"))

setGeneric("computePotential", function(object) standardGeneric("computePotential"))

setGeneric("dat", function(object) standardGeneric("dat"))
setGeneric("dat<-", function(object,value) standardGeneric("dat<-"))

setGeneric("updateMu", function(object) standardGeneric("updateMu"))

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
setGeneric("updateNu.0", function(object) standardGeneric("updateNu.0"))
setGeneric("updateTau2", function(object) standardGeneric("updateTau2"))

setGeneric("alpha<-", function(object, value) standardGeneric("alpha<-"))

setGeneric("updateWithPosteriorMeans", function(object) standardGeneric("updateWithPosteriorMeans"))

#' @export
setGeneric("bic", function(object, ...) standardGeneric("bic"))

#' @export
setGeneric("theta", function(object) standardGeneric("theta"))

#' @export
setGeneric("sigma2", function(object) standardGeneric("sigma2"))

setGeneric("reorderComponents", function(object) standardGeneric("reorderComponents"))

#' @export
setGeneric("hwe", function(object) standardGeneric("hwe"))

#' @export
setGeneric("probz", function(object) standardGeneric("probz"))

setGeneric("probz<-", function(object, value) standardGeneric("probz<-"))

#' @export
setGeneric("fitMixtureModels", function(object, mcmcp, K=1:5) standardGeneric("fitMixtureModels"))

#' @export
setGeneric("nu.0", function(object) standardGeneric("nu.0"))

#' @export
setGeneric("sigma2.0", function(object) standardGeneric("sigma2.0"))

#' @export
setGeneric("y", function(object) standardGeneric("y"))

#' @export
setGeneric("z", function(object) standardGeneric("z"))

setGeneric("modalParameters", function(object) standardGeneric("modalParameters"))

setGeneric("computeModes", function(object) standardGeneric("computeModes"))

setGeneric("switchLabels", function(object) standardGeneric("switchLabels"))

setGeneric("computeDistance", function(object) standardGeneric("computeDistance"))

setGeneric("modes", function(object) standardGeneric("modes"))
setGeneric("modes<-", function(object,value) standardGeneric("modes<-"))


setGeneric("mu.0", function(object) standardGeneric("mu.0"))
setGeneric("tau2.0", function(object) standardGeneric("tau2.0"))
setGeneric("eta.0", function(object) standardGeneric("eta.0"))
setGeneric("m2.0", function(object) standardGeneric("m2.0"))

setGeneric("showMeans", function(object) standardGeneric("showMeans"))
setGeneric("showSigmas", function(object) standardGeneric("showSigmas"))

#' @export
setGeneric("collapseBatch", function(object) standardGeneric("collapseBatch"))

#' @export
setGeneric("thetac", function(object) standardGeneric("thetac"))

#' @export
setGeneric("thetaMean", function(object) standardGeneric("thetaMean"))

#' @export
setGeneric("sigmaMean", function(object) standardGeneric("sigmaMean"))

#' @export
setGeneric("pMean", function(object) standardGeneric("pMean"))

#' @export
setGeneric("tracePlot", function(object, name, ...) standardGeneric("tracePlot"))

setGeneric("tablez", function(object) standardGeneric("tablez"))

setGeneric("nStarts", function(object) standardGeneric("nStarts"))
setGeneric("checkLabels", function(object) standardGeneric("checkLabels"))
setGeneric("nStartIter", function(object) standardGeneric("nStartIter"))

setGeneric("updateAlpha", function(object) standardGeneric("updateAlpha"))
setGeneric("alpha", function(object) standardGeneric("alpha"))

setGeneric("orderTheta<-", function(object, value) standardGeneric("orderTheta<-"))
setGeneric("orderTheta", function(object) standardGeneric("orderTheta"))
