context("Down sampling")

## For large datasets, downsample observations. This makes it easier to run a
## large number of MCMC simulations for big datasets.

.test_that <- function(name, expr) NULL

.test_that("downsample", {
  x <- .read_hapmap()
  xx <- do.call(cbind, x)/1000
  sds <- matrixStats::rowSds(xx)
  meds <- matrixStats::rowMedians(xx)
  ## bad probes
  index <- which(abs(meds) > 0.01 )
  xx <- xx[-index, ]

  pc <- prcomp(xx)
  pc1 <- pc$rotation[, "PC1"]
  hist(pc1, breaks=100)

  ##meds <- colMedians(xx)
  dat <- readLocalHapmap()

  mp <- McmcParams(iter=1000, burnin=1000,
                   nStarts=4, thin=10)
  batches <- collapseBatch(dat, names(dat), THR=0.05) %>%
    factor %>% as.integer
  hp.list <- hyperparam_list()

  hp <- hp.list[["single_batch"]]
  k(hp) <- 4
  nStarts(mp) <- 1L
  burnin(mp) <- 10000L
  iter(mp) <- 0L
  set.seed(123)
  sb <- replicate(20, SingleBatchModel(hp=hp,
                                       mp=mp,
                                       dat=dat))
  tmp <- map(sb, posteriorSimulation)
  set_mp <- function(model, mp){
    mcmcParams(model) <- mp
    model
  }
  iter(mp) <- 1000L
  burnin(mp) <- 0L
  tmp2 <- map(tmp, set_mp, mp)
  tmp3 <- map(tmp2, posteriorSimulation)

  mp <- McmcParams(iter=1000L,
                   thin=10L,
                   nStarts=4L,
                   burnin=10000L)
  sb <- gibbs(mp=mp, dat=dat, hp=hp)


  tmp <- posteriorSimulation(sb)
  nStarts(mp) <- 4L

  tmp2 <- gibbs(dat=dat, mp=mp, hp=hp)
  models <- gibbs_all(dat=dat,
                      batches=batches,
                      mp=mp,
                      hp.list=hp.list)
  cn.model <- SingleBatchCopyNumber(models[[1]])
  p1 <- ggSingleBatch(models[[1]])
  p2 <- ggSingleBatch(cn.model)
  library(grid.arrange)
  grid.arrange(p1, p2)

  hp <- hp.list[["single_batch"]]
  k(hp) <- 4
  mp2 <- mp
  iter(mp2) <- 0L
  mm.list <- replicate(4, SingleBatchModel2(dat=dat, mp=mp, hp=hp))
  mm.list <- map(mm.list, posteriorSimulation)
  mm <- posteriorSimulation(mm)

  sb.models <- gibbs(dat=dat, mp=mp2, hp=hp)



  sb <- SingleBatchModelList(dat, k=2:6, mcmc.params=mp)
  sb <- posteriorSimulation(sb)
  ml <- marginalLikelihood(sb)
  select <- which.max(ml)
  select <- useModes(sb[[select]])
  ggSingleBatch(select)
  cn.model <- SingleBatchCopyNumber(select)
  ##
  ## Should we have a merging parameter that detects 'outlier' components
  ## (components with large variance) and prevent merging the outlier component
  ## with other (copy number) components? This would prevent assigning the same
  ## copy number to extreme values (negative and positive)
  ##
  map <- mapComponents(cn.model)
  mapping(cn.model) <- map
  if(FALSE)
    ggSingleBatch(cn.model)

  ##trace(CNPBayes:::mapCopyNumber, browser)
  params <- mapParams()
  CNPBayes:::mapCopyNumber(cn.model)
  map <- mapCopyNumber(cn.model)
  mapping(cn.model) <- map
  if(FALSE){
    trace(.gg_singlebatch_copynumber, browser)
    .gg_singlebatch_copynumber(cn.model)
    ggSingleBatch(cn.model)
  }

  ##cn.model <- SingleBatchCopyNumber(sb[[5]])
  ##map <- mapComponents(cn.model)
  ##mapping(cn.model) <- map
  ##ggSingleBatch(cn.model)
  ml <- marginalLikelihood(sb[[3]],
                           mlParams(ignore.effective.size=TRUE,
                                    ignore.small.pstar=TRUE))
  ml <- marginalLikelihood(sb[[4]],
                           mlParams(ignore.effective.size=TRUE,
                                    ignore.small.pstar=TRUE))
  ml <- marginalLikelihood(sb[[5]],
                           mlParams(ignore.effective.size=TRUE,
                                    ignore.small.pstar=TRUE))

  set.seed(134)
  dat <- readLocalHapmap()
  plt <- names(dat)
  dat2 <- data.frame(y=dat, plate=plt)
  ggplot(dat2, aes(y)) + geom_histogram(bins=50) +
    facet_wrap(~plate)

  b <- collapseBatch(dat, names(dat))
  mp <- McmcParams(iter=1000, burnin=1000, nStarts=20, thin=10)
  dat2 <- downSampleEachBatch(dat, 50, b)
  sb <- SingleBatchModelList(dat2$y, k=2:5, mcmc.params=mp)
  sb <- posteriorSimulation(sb)
  ggSingleBatch(sb[[2]])
  ml <- marginalLikelihood(sb)


  sb2 <- sb[[2]]
  sb.cn <- SingleBatchCopyNumber(sb2)
  map <- mapComponents(sb.cn, mapParams(threshold=0.2))
  mapping(sb.cn) <- map
  ##
  ## observations on both tails captured by 3rd component
  ## -- need to merge 3 with 2 and 1, but not merge 1 and 2
  ggSingleBatch(sb.cn)


  ##
  ## pseudocode
  ##
  model.fulldata <- SingleBatchModel(data)
  model.partdata <- downSample(model.fulldata)
  model.partdata <- posteriorSimulation(model.partdata)
  model.fulldata1 <- upSample(model.partdata)

  model.fulldata2 <- posteriorSimulation(model.fulldata)

  ## compare z's
  ## compare thetas
  ## compare likelihood
  identical(model.fulldata1, model.fulldata2)
})
