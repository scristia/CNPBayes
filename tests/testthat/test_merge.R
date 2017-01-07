context("merging")

test_that("simple merge", {
  set.seed(1)
  ##
  ## IDEA: There is only one copy number state, but data is not quite normal
  ##
  ##  - simulate two components, but with substantial overlap
  ##
  ##  - verify best model should identify two components
  ##
  ##  - check whether best model should be merged
  ##
  ##  - merge the components, keeping all probabilistic estimates intact
  ##
  truth <- simulateData(N=100, p=c(0.9, 0.1),
                        theta=c(0, 0.25), sds=c(0.2, 0.2))
  mp <- McmcParams(iter = 500, burnin = 150, nStarts = 25)
  model <- MarginalModel(data = y(truth), k = 2, mcmc.params = mp)
  model <- posteriorSimulation(model)

  ## why doesn't probz get updated for first iteration
  m <- MarginalModel(data = y(truth), k = 2,
                     mcmc.params = McmcParams(iter=2, burnin=0))
  m <- posteriorSimulation(m)

  ## p is k-length vector of posterior probabilites for one observation (sample)
  expect_identical(.merge_prob(c(0.99, 0.01, 0)), 1:3)
  expect_identical(.merge_prob(c(0.8, 0.2, 0)), 1:3)
  expect_identical(.merge_prob(c(0.8, 0.2, 0)), 1:3)
  expect_identical(.merge_prob(c(0.8, 0.2, 0), threshold=0), rep(1L, 3))
  expect_identical(.merge_prob(c(0.8, 0.2, 0), threshold=0.15),
                   c(1L, 1L, 3L))


  ##
  ## Another way to frame this is that we'd like to define a new mixture model
  ## with fewer components such that the majority of posterior probabilities are
  ## very high
  ##
  trace(CNPBayes:::assessMerge, browser)
  ## merge components 1 and 2
  expect_identical(assessMerge(probz(model), params=mergeParams()), 1:2)

  p <- cbind(rep(1, 10), rep(0, 10), rep(0, 10))
  ## no merging
  expect_true(is.null(assessMerge(p, params=mergeParams())))

  ##trace(CNPBayes:::assessMerge, browser)
  ## merge components 2 and 3
  p <- cbind(rep(0.01, 10), rep(0.45, 10), rep(0.54, 10))
  expect_identical(assessMerge(p, params=mergeParams()), 2:3)

  ## merge all components
  p <- cbind(rep(0.01, 10), rep(0.45, 10), rep(0.54, 10))
  trace(CNPBayes:::assessMerge, browser)
  expect_identical(assessMerge(p, params=mergeParams(threshold=0)), 1:2)
  expect_true(is.null(assessMerge(p, params=mergeParams(threshold=1))))

  params <- mergeParams(threshold=0.1)
  model2 <- .merge_components(model, assessMerge(probz(model), params))
  expect_equal(probz(model2), cbind(rep(0, length(y(model))), 1))

  model2 <- mergeComponents(model, params)
  expect_equal(probz(model2), cbind(rep(0, length(y(model))), 1))

  ## three compomnents -- merge all
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  trace(CNPBayes:::singleBatchDensities, browser)
  CNPBayes:::singleBatchDensities(truth)

  trace(CNPBayes::ggSingleBatch, browser)
  CNPBayes::ggSingleBatch(truth)

  mp <- McmcParams(iter = 500, burnin = 150, nStarts = 25)
  model <- MarginalModel(data = y(truth), k = 2, mcmc.params = mp)
  model <- posteriorSimulation(model)




  expect_identical(p(model2), 1)
  trace(CNPBayes:::singleBatchDensities, browser)
  tmp <- CNPBayes:::singleBatchDensities(model2)

  expect_true(all(z(model2) == 1L))
  expect_true(all(probz(model2) == 1))

    .merge_prob(frac.obs.uncertain, cutoff)
    new.labels <- mergeLabels(p.z, params)
    ##
    ##
    ##
    fraction.relabeled <- 


    posterior.prob <- colMeans(p.z)
    if(all(posterior.prob > threshold)){
      return(model)
    }
    ## combine components with substantial overlap
    zz <- z(chains(model))
    ## suppose 4-comp model
    ## - components 1 and 2 overlap
    ## - components 3 and 4 overlap
    ## posterior.probs:
    ## 0.34, 0.66, 0.55, 0.45
    ## 
    ## Can tell apart by the p.z matrix
    ##  y[1]  0.34, 0.66, 0, 0
    ##  y[2]  0,  0, 0.55, 0.45
    ##
    ## Need function that merges for each z.
    apply(p.z, 1)
    components.to.combine <- which(posterior.prob < threshold)
    ## recursive?
    remove.component <- max(components.to.combine)
    expand.component <- components.to.combine[]


  }

})



##
## This is a tough example. Best approach is unclear.
##
smallPlates <- function(x){
  tab <- table(x)
  names(tab)[tab < 20]
}

readLocalHapmap <- function(){
  ddir <- "~/Dropbox/labs/cnpbayes"
  lrr <- readRDS(file.path(ddir, "data/EA_198_lrr.rds"))
  lrr1 <- lapply(lrr, function(x) x/1000)
  batch.id <- c(rep(0,8), rep(1, 8))
  avg.lrr <- unlist(lapply(lrr1, colMeans, na.rm=TRUE))
  plate <- substr(names(avg.lrr), 1, 5)
  avg.lrr <- avg.lrr[!plate %in% smallPlates(plate)]
  plate <- plate[!plate %in% smallPlates(plate)]
  names(avg.lrr) <- plate
  avg.lrr
}

mclustMeans <- function(y, batch){
  ylist <- split(y, plates2)
  .mclust <- function(y){
    Mclust(y)$parameters$mean
  }
  mns <- lapply(ylist, .mclust)
  L <- sapply(mns, length)
  collections <- split(names(L), L)
}

.test_that <- function(expr, name) NULL

.test_that("hapmap", {
  set.seed(134)
  dat <- readLocalHapmap()
  b <- collapseBatch(dat, names(dat))
  mp <- McmcParams(iter=1000, burnin=500, nStarts=20)
  ml <- BatchModelList(dat, k=2:5, batch=b, mcmc.params=mp)
  ml <- posteriorSimulation(ml)
  ggMultiBatchChains(ml[[4]])[["batch"]]

  sb <- MarginalModelList(dat, k=4:8, mcmc.params=mp)
  sb <- posteriorSimulation(sb)
  tmp <- sample(dat, length(dat), replace=TRUE)
  ggSingleBatchChains(sb[[2]])[["comp"]]
  ggSingleBatch(sb[[3]])
  ggSingleBatch(sb[[4]])
  ggSingleBatch(sb[[5]])

  ggSingleBatch(model)
  ## evaluate merging for k=4
  m4 <- mlist[[3]]
  ggSingleBatch(m4)
  ##
  ## here, component 2 has a large variance
  ##
  ggSingleBatchChains(m4)[["comp"]]





  model <- mlist[[select]]
  d <- densities(model)
  dc <- densitiesCluster(model, merge=TRUE)
  dmlist <- lapply(mlist, DensityModel, merge=TRUE)
  n.comp <- sapply(dmlist, function(x) length(modes(x)))
  ## remove merge models where number components are duplicated
  mlist <- mlist[!duplicated(n.comp)]
  m.y <- marginalLikelihood(mlist)##, params=params)
  argmax <- which.max(m.y)
  expect_true(argmax == 2L)
  if(FALSE){
    plist <- ggSingleBatchChains(mlist[[2]])
    plist[["comp"]]

    plist3 <- ggSingleBatchChains(mlist[[3]])
    plist3[["comp"]]

    ggSingleBatch(mlist[[3]])
    ggSingleBatch(mlist[[2]])

    pstar <- marginal_theta(mlist[[2]])
  }
})
