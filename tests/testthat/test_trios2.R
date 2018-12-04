context("new trios class")

test_that("Michael's example", {
  skip("Michael's code requires triostat package")
  library(devtools)
  library(coda)
  library(tidyverse)
  library(bindrcpp)
  library(triostat)
  p <- c(0.09, 0.42, 0.49)
  theta <- c(-3, 0.15, 1.2)
  sigma2 <- c(0.2, 0.2, 0.2)
  ##sigma2 <- c(0.1, 0.1, 0.1)
  ##sigma2 <- c(0.05, 0.05, 0.05)
  N <- 1000
  params <- data.frame(cbind(p, theta, sigma2))
  seed <- 4990917
  set.seed(seed)
  # note maplabel and hp defined manually here for now
  maplabel <- c(0,1,2)
  k <- length(maplabel)
  hp <- HyperparametersTrios(k = 3)
  nbatch <- 1
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                length.out = 3*N),
                                error=0, mprob,
                                maplabel)
  mp <- McmcParams(iter=300, burnin=200, thin=1)
  base <- TBM(triodata=as.tibble(truth$data),
           hp=hp,
           mp=mp,
           mprob=mprob,
           maplabel=maplabel)
  tdat <- truth$data %>%
    mutate(copy_number=factor(copy_number))
  ggplot(tdat, aes(copy_number, log_ratio)) +
    geom_jitter(data=tdat[is_parent(base), ],
                width=0.2, color="gray") +
    geom_jitter(data=tdat[!is_parent(base), ],
                width=0.2, color="steelblue")

  expect_true(validObject(base))
  z1 <- z(base)
  zz <- update_zchild(base)
  is_offspr <- is_child(base)
  expect_identical(z1[ !is_offspr], zz[!is_offspr])
  expect_true(!identical(z1[ is_offspr], zz[is_offspr]))

  m <- posteriorSimulation(base)
  if(FALSE){
    ggMixture(m)
    prob.mendelian <- probMendelian(m)
    dat <- tibble(prob=prob.mendelian, iter=seq_along(prob.mendelian))
    ggplot(dat, aes(iter, prob)) +
      geom_point() +
      xlab("Index for trio") +
      ylab("Prob transmission is Mendelian") +
      ylim(c(0, 1)) +
      geom_hline(yintercept=0.9)
  }


  model <- base
  ptrio <- update_trioPr2(model)
  expect_identical(dim(ptrio), c(1000L, 3L))
  expect_true(all(rowSums(ptrio) == 1))
  expect_identical(is_father(model), triodata(model)$family_member=="f")
  expect_identical(is_mother(model), triodata(model)$family_member=="m")
  p2 <- update_multinomialPrChild(model)
  expect_equal(rowSums(p2), rep(1, 1000))
  zo <- update_offspring(model)
  expect_true(all(zo > 0))
  zz <- update_zchild(model)
  expect_true(all(zz > 0))
  m <- trios_burnin(base)
  expect_true(validObject(m))
  m <- trios_mcmc(m)
  expect_true(validObject(m))
  m2 <- sortComponentLabels(m)
  expect_true(validObject(m2))
  m2 <- posteriorSimulation(m2)
  expect_true(validObject(m2))
  m <- m2
  ##
  ## These tests fail in current implementation
  ##
  expect_equal(theta(m)[1, ], theta, tolerance=0.05)
  expect_equal(sigma2(m)[1, ], sigma2, tolerance=0.1)
  expect_equal(p(m), p, tolerance=0.05)
  pm <- probMendelian(m)
  ## Since observations were simulated under mendelian transmission, we would
  ## expect this to be high
  expect_true(mean(pm > 0.9) > 0.75)


  ##
  ## force the independence model by not updating is_mendelian
  ## and commenting relevant updates
  m@is_mendelian <- rep(0L, nTrios(m))
  m <- burnin_nomendelian_update(m)
  m <- mcmc_nomendelian_update(m)
  expect_true(all(isMendelian(m) == 0))
  expect_equal(theta(m)[1, ], theta, tolerance=0.05) ## WORKS!
  ## if offspring were generated under mendelian transmission
  ## and parents are updated in the usual way, why do we do worse?
  triodat <- select(truth$data, c(id, copy_number, family_member)) %>%
    spread("family_member", "copy_number") %>%
    ## find the unique combinations
    unite("of", c("o", "f")) %>%
    unite("ofm", c("of", "m")) %>%
    group_by(ofm) %>%
    summarize(n=n())

  ## force mendelian updates
  m3 <- m
  burnin(m3) <- 1000
  m3@is_mendelian <- rep(1L, nTrios(m))
  m3 <- burnin_nomendelian_update(m3)
  ##m3 <- mcmc_nomendelian_update(m3)
  expect_true(all(isMendelian(m3) == 1))
  ## Fails
  expect_equal(theta(m3)[1, ], theta, tolerance=0.05)
  truez <- truth$data$copy_number + 1
  table(z(m3), truez)
  if(FALSE) ggMixture(m3)
  ## posterior predictive distribution looks OK
  ## however, sigma2 values are way off
  model <- m
  model@is_mendelian[] <- 1L
  for(i in 1:10){
    z(model) <- update_z(model)
    ##z(model) <- update_zchild(model)
    model@zfreq_parents <- tableZpar(model)
    sigma2(model) <- update_sigma2(model)
    print(sigma2(model))
    nu.0(model) <- update_nu0(model)
    sigma2.0(model) <- update_sigma20(model)
    ##model@is_mendelian = update_mendelian(model)
    model@zfreq = tableZ(k(model), z(model))
    theta(model) <- update_theta(model)
    mu(model) <- update_mu(model)
    tau2(model) <- update_tau2(model)
    p(model) <- update_p(model)
    model@pi_parents <- update_pp(model)
    model@u <- rchisq(3000, 100)
  }
  table(z(model), truez)
  table(z(model)[is_child(model)], truez[is_child(model)])

  ## variances remain stable when the child z's are not updated
  model2 <- model
  for(i in 1:10){
    z(model2) <- update_z(model2)
    ##z(model2) <- update_zchild(model2)
    model2@zfreq_parents <- tableZpar(model2)
    sigma2(model2) <- update_sigma2(model2)
    print(sigma2(model2))
    nu.0(model2) <- update_nu0(model2)
    sigma2.0(model2) <- update_sigma20(model2)
    model2@is_mendelian = update_mendelian(model2)
    model2@zfreq = tableZ(k(model2), z(model2))
    theta(model2) <- update_theta(model2)
    mu(model2) <- update_mu(model2)
    tau2(model2) <- update_tau2(model2)
    p(model2) <- update_p(model2)
    model2@pi_parents <- update_pp(model2)
    model2@u <- rchisq(3000, 100)
  }
  table(z(model), truez)
  table(z(model2), truez)
  table(z(model)[is_child(model)], truez[is_child(model)])
  table(z(model2)[is_child(model2)], truez[is_child(model2)])

  true.zo <- select(tdat, c(id, family_member, copy_number)) %>%
    spread(family_member, copy_number) %>%
    set_colnames(c("id", "true_m", "true_f", "true_o"))
  lr <- select(tdat, c(id, family_member, log_ratio)) %>%
    spread(family_member, log_ratio)
  truedat <- left_join(true.zo, lr, by="id")
  bayes1a <- bayesm <- select(tdat, c(id, family_member))
  is_off <- is_child(model2)
  model3 <- model2
  for(i in 1:10){
    z(model3) <- update_z(model3)
    tdat <- tdat %>%
      mutate(zz=factor(z(model3)),
             lr=log_ratio)
    ggplot(data=tdat, aes(zz, lr)) +
      geom_jitter(data=tdat[is_parent(base), ],
                  width=0.2, color="gray") +
      geom_jitter(data=tdat[!is_parent(base), ],
                  width=0.2, color="steelblue")
    z(model3) <- update_zchild(model3)
    ##
    ## Samples that were homozygously deleted are called hemizygous
    ## -- this increases the variance of the hemizygous component
    ##
    ## Suppose we have cf,cm = [2, 1]  (offspring is incorrect)
    ##   if M = 1, then p(c0 = 0 | data) is zero
    ##   if M = 0, then p(c0 = 0 | data) is near 1
    ##       the independence model is correct if assumption of ind. is true
    ##       the dependence model is always correct
    ##       but this does not give us the update we want
    ##       we kind of want to do the opposite:
    ##           - the independence model is correct
    ##           - the dependence model is incorrect
    ##   solutions:
    ##     - block update for cf, cm, co, M
    ##     - only accept cf, cm, co, M if co = 0
    ##
    ## - the probability of the dependent model is high
    ## - M is updated as 1
    ## - at some point we see [1, 1]
    ## - offspring will almost certainly be called 0
    ## - mendelian model may be favored a priori
    ## - parents get [2, 1]
    ## - offspring forced to zero
    ##
    ## intuition:
    ##
    ##  if offspring is clearly homozygous and the dependence model does not
    ##  allow this, the indicator for dependency should be updated to have value zero and the offspring should be updated to have zero copies
    ##  - dependence indicator should integrate over the possible copy number states
    ##
    ##
    ##     - the probability of the dependent model is 0
    ##       =>  M will be updated to have value 0
    ##       =>  offspring will most certainly be called 0 at next iteration
    ##       =>  M will be updated to have value 0
    ##     - at some point we could observe [1, 1, 0] (offspring correctly called)
    ##     - dependent model would be favored; M = 1
    ##     - parents could be updated to [2, 1]
    ##     - offspring forced to be 1  (incorrect)
    ##
    ##
    ## -- wouldn't the mendelian probability be small in these scenarios?
    ## p(M=1|co,cf,cm) = p(co,cf,cm|1)p(M=1)/p(co, cf, cm)
    ##                 = p(co|cf,cm,1)p(cf,cm|1)p(1)/p(co, cf, cm)
    ##            propto p(co|cf,cm,1)p(1)/p(co, cf, cm)
    tdat <- tdat %>%
      mutate(zz=factor(z(model3)),
             lr=log_ratio)
    ggplot(data=tdat, aes(zz, lr)) +
      geom_jitter(data=tdat[is_parent(base), ],
                  width=0.2, color="gray") +
      geom_jitter(data=tdat[!is_parent(base), ],
                  width=0.2, color="steelblue")
    id <- filter(tdat, log_ratio < -2.5 & zz == 2) %$%
      id %>%
      as.character
    ptrio=update_trioPr(model)
    pi <- model@pi_parents
    p <- prob_mendelian(model3)
    p [ unique(tdat$id) %in% id ]
    model3@is_mendelian = update_mendelian(model3)
    ismendel <- setNames(isMendelian(model3),
                         unique(tdat$id))
    ismendel[ names(ismendel) %in% id]

    model3@zfreq_parents <- tableZpar(model3)
    sigma2(model3) <- update_sigma2(model3)
    print(sigma2(model3))
    nu.0(model3) <- update_nu0(model3)
    sigma2.0(model3) <- update_sigma20(model3)
    model3@zfreq = tableZ(k(model3), z(model3))
    theta(model3) <- update_theta(model3)
    mu(model3) <- update_mu(model3)
    tau2(model3) <- update_tau2(model3)
    p(model3) <- update_p(model3)
    model3@pi_parents <- update_pp(model3)
    model3@u <- rchisq(3000, 100)
  }


  ##
  ## Why would the posterior predictive distribution look OK if sigma2 is exploding?
  ##
  dat1 <- filter(triodata(model), family_member=="o",
                copy_number==1) %>%
    mutate(batch=paste0("Batch ", batch))
  dat2 <- filter(triodata(model), family_member=="o",
                copy_number==2) %>%
    mutate(batch=paste0("Batch ", batch))
  ## The rug plots suggest that the offspring calls are reasonable
  ## -- nothing to indicate why sigma2 would explode
  ggMixture(model2) +
    geom_rug(data=dat1,
             aes(x=log_ratio), inherit.aes=FALSE)
  ggMixture(model2) +
    geom_rug(data=dat2,
             aes(x=log_ratio), inherit.aes=FALSE)


  datp1 <- filter(triodata(model), family_member=="f" | family_member=="m",
                  copy_number==1) %>%
    mutate(batch=paste0("Batch ", batch))
  datp2 <- filter(triodata(model), family_member=="f" | family_member=="m",
                  copy_number==2) %>%
    mutate(batch=paste0("Batch ", batch))
  ggMixture(model2) +
    geom_rug(data=datp1,
             aes(x=log_ratio), inherit.aes=FALSE)
  ggMixture(model2) +
    geom_rug(data=datp2,
             aes(x=log_ratio), inherit.aes=FALSE)


  if(FALSE) ggMixture(m2)
  ##
  ## When isMendelian is 1, the normal density is weighted by the mendelian transmission prob
  ## When isMendelian is 0, the normal density is weighted by the usual mixture probabilities
  ## The mendelian indicator is updated according to the likelihood of c_o under the mendelian model multiplied by the prior and the likelihood of c_o under the population mixture proportions.
  ##
  ## The full conditional for the probability Mi=0 is beta(n1 + 1, n1+n2 + 1), where
  ##
  m1 <- posteriorSimulation(model[[1]])
  if(FALSE){
    ggMixture(m1)
    prob.mendelian <- isMendelian(chains(m1))/(iter(m1) + iter(model[[1]]))
    dat <- tibble(prob=prob.mendelian, iter=seq_along(prob.mendelian))
    ggplot(dat, aes(iter, prob)) +
      geom_point() +
      xlab("Index for trio") +
      ylab("Prob transmission is Mendelian")
  }
  mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
                              batches=batch(model[[1]]),
                              mp=mp,
                              k_range=c(3, 3),
               max_burnin=100)
  mb2 <- mb2[[1]]
  theta(mb2)
  sigma(mb2)
  expect_equivalent(hyperParams(mb2), hyperParams(m1))
  expect_equal(sigma2(mb2)[1, ], sigma2, tolerance=0.1)
  expect_equal(p(mb2), p, tolerance=0.05)

  expect_equal(theta(mb2)[1, ], theta, tolerance=0.05)

  dat <- filter(truth$data, family_member != "o")
  mb <- gibbs(model="MB", dat=dat$log_ratio,
              batches=dat$batches,
              mp=mp,
              k_range=c(3, 3),
              max_burnin=100)
  mb <- mb[[1]]

  tmod <- model[[1]]
  theta(tmod) <- theta(mb)
  sigma2(tmod) <- sigma2(mb)
  p(tmod) <- p(mb)
  mu(tmod) <- mu(mb)
  tau2(tmod) <- tau2(mb)
  z(tmod) <- z(mb)
  sigma2.0(tmod) <- sigma2.0(mb)
  mp <- McmcParams(burnin=150, iter=200, nStarts=1)
  hp <- hyperParams(tmod)
  tmod2 <- TBM(triodata=truth$data,
               hp=hp,
               mp=mp,
               mprob=mprob,
               maplabel=maplabel)
  tmod2 <- posteriorSimulation(tmod2)
  ggMixture(tmod2)
  figs <- ggChains(tmod2)



  ggMixture(model[[1]])
  ggChains(model[[1]])
  ggMixture(mb2[[1]])
  ggChains(mb2[[1]])

  trio.results <- sum.fit(model, truth, maplabel)
  mb.results <- sum.fit(mb2, truth, maplabel)

  # other random stats
  true.stats <- component_stats(truth$data)
  truth.parent.pi <- table(trio.results$cn.truth.parent)/length(trio.results$cn.truth.parent)
  truth.child.pi <- table(trio.results$cn.truth.child)/length(trio.results$cn.truth.child)
  child.pi <- table(trio.results$results.child)/length(trio.results$results.child)
  parent.pi <- table(trio.results$results.parent)/length(trio.results$results.parent)
  logrr <- truth$data$log_ratio

  # input all the info
  summaryResults <- list(params = params,
                         sim.params = truth$params,
                         N = length(logrr)/3,
                         SimLogRR = logrr,
                         SimTruth = true.stats,
                         SimPi = true.stats$p,
                         SimParPi = truth.parent.pi,
                         SimChildPi = truth.child.pi,
                         SimThetas = true.stats$mean,
                         SimVar = (true.stats$sd)^2,
                         TruthCNcall = truth$data$copy_number,
                         TruthCNpar = trio.results$cn.truth.parent,
                         TruthCNoff = trio.results$cn.truth.child,
                         TrioCNcall = trio.results$results@z,
                         TrioCNpar = trio.results$results.parent,
                         TrioCNoff = trio.results$results.child,
                         MBCNcall = mb.results$results@z,
                         MBCNpar = mb.results$results.parent,
                         MBCNoff = mb.results$results.child,
                         Accuracy = trio.results$prop.true.overall,
                         AccuracyParents = trio.results$prop.true.parent,
                         AccuracyChild = trio.results$prop.true.child,
                         AccuracyMB = mb.results$prop.true.overall,
                         AccuracyMBParents = mb.results$prop.true.parent,
                         AccuracyMBChild = mb.results$prop.true.child,
                         ModelPi = model[[1]]@modes$mixprob,
                         ModelParentPi = parent.pi,
                         ModelChildPi = child.pi,
                         ModelTheta = model[[1]]@modes$theta,
                         ModelSigma2 = model[[1]]@modes$sigma2
  )

  data.read <-  simdata.read(summaryResults)
  stables <- summarise_results(data.read)
  simresults2 <- overall_results2(stables)
  simresults2
})
