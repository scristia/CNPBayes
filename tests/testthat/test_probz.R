context("Check z's")

test_that("test_probz", {
    set.seed(1)
    truth <- simulateData(N = 2500, p = rep(1/3, 3), theta = c(-1, 
        0, 1), sds = rep(0.1, 3))
    mp <- McmcParams(iter = 500, burnin = 500, nStarts=0)
    set.seed(123)
    model <- SingleBatchModel2(dat = y(truth), mp = mp,
                               hp=hpList(k=3)[["SB"]])
    model <- startAtTrueValues(model, truth)
    pz <- probz(model)
    expect_true(all(pz == 0))
    true_z <- z(truth)
    expect_equal(z(model), z(truth))
    if(FALSE){
      ##model <- posteriorSimulation(model)
      ##pz <- probz(model)
      iter(model, force = TRUE) <- 2L
      burnin(model) <- 0L
      zz <- map_z(posteriorSimulation(model))
      expect_equal(true_z, zz)
      model2 <- modelOtherModes(model, maxperm = 2)[[2]]
      z2 <- zz <- z(model2)
      expect_true(sum(zz != true_z) > 500)
      mcmcParams(model2) <- mcmcParams(model)
      model2 <- posteriorSimulation(model2)
      zz <- map_z(model2)
      expect_equal(true_z, zz)
      model3 <- modelOtherModes(model, maxperm = 3)[[3]]
      z3 <- z(model3)
      expect_true(sum(z3 != true_z) > 500)
      expect_true(sum(z2 != z3) > 500)
      model3 <- posteriorSimulation(model3)
      mz3 <- map_z(model3)
      ##table(mz3, true_z)
      expect_equal(true_z, mz3)
    }
})

