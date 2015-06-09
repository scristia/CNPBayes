.test_profiling <- function(){
  ##
  ## example (from test_marginal.R)
  ##
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  if(FALSE) plot(truth, use.current=TRUE)
  params <- ModelParams("marginal", y=y(truth), k=3)
  mcmcp <- McmcParams(iter=150, burnin=0)
  ##trace(initializeModel, browser)
  model <- CNPBayes:::initializeModel(params)
  
  Rprof("m2.prof")
  model <- posteriorSimulation(model, mcmcp)
  Rprof(NULL)
  m1 <- summaryRprof("m1.prof")
  m2 <- summaryRprof("m2.prof")
  
  ##
  ## example (from test_batch.R)
  ##
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                    -0.2, 0, 0.2,
                    0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  truth <- simulateBatchData(N=2500,
                             batch=rep(letters[1:3], length.out=2500),
                             theta=means,
                             sds=sds,
                             p=c(1/5, 1/3, 1-1/3-1/5))
  mcmcp <- McmcParams(iter=100, burnin=0)
  ## use the true k
  k <- 3
  ## marginal
  params <- ModelParams("marginal", y=y(truth), k=k,
                        batch=rep(letters[1:3], length.out=2500),
                        mcmc.params=mcmcp)
  model <- CNPBayes:::initializeModel(params)
  #model <- posteriorSimulation(model, mcmcp)
  updateTheta(model, constrain=FALSE)
  updateThetaCpp(model, constrain=TRUE)
  updateSigma2Cpp(model)
  ## Need to profile posteriorSimulation.
  Rprof()
  posteriorSimulation(model, mcmcp)
  
  ##
  ## Batch model
  ##
  outdir <- '~/Software'
  se472 <- readRDS("~/Software/CNPBayes/inst/extdata/se_cnp472_EA.rds")
  B <- getFiles(outdir, rownames(se472), "batch")
  
  batch.files <- paste0(dirname(model(B)), "/", rownames(se472), "_batch.rds")
  object <- BatchModel(copyNumber(se472)[1,], k=3, batch= batch.files)
  
  load_all()
  Rprof()
  batchExperiment(se472, outdir, test=TRUE)
  Rprof(NULL)
  summaryRprof()
  prof <- summaryRprof()
  tot <- prof$by.total
  self <- prof$by.self
  
  load_all()
  batchExperiment(se472, outdir, test=TRUE)
  outdir <- tempdir()
  Rprof("marginal.prof", interval=0.1)
  marginalExperiment(se472, outdir)
  Rprof(NULL)
  mp2 <- summaryRprof("marginal.prof")
  
  load_all()
  marginalExperiment(se472, outdir, test=TRUE)
  
  
  
  #mcmcp <- McmcParams(iter=c(500, 2000, 2000, 2000),
  #                      burnin=c(50, 100, 100, 100),
  #                      thin=c(1, 2, 2, 2),
  #                      nStarts=20,
  #                      nStartIter=200)
  #  hypp <- HyperparametersBatch()
  #  hplist <- HyperParameterList(hypp, K=1:4)
  #  mplist <- ModelParamList(hypp, K=1:4, data=copyNumber(se472)[1, ],
  #                           batch=se472$batch, mcmcp=mcmcp)
  #  bmodlist <- foreach(hypp=hplist, param=mplist) %do% {
  #    initializeBatchModel(params=param, hypp=hypp)
  #  }
  #  bmodels <- foreach(k=1:4, model=bmodlist) %do% posteriorSimulation(model, mcmcp[k])
  
  ##### Rcpp S4 testing stuff
  library(inline)
  src <- '
  S4 foo(x) ; foo.slot(".Data") = "bar" ; return(foo);
  '
  fun <- cxxfunction(signature(x="any"), src,
                     plugin="Rcpp")
  setClass( "S4ex", contains = "character",
  representation( x = "numeric" ) )
  x <- new( "S4ex", "bla", x = 10 )
  fun(x)
  str(fun(x))
}
