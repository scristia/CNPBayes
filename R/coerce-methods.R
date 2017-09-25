#' @include Deprecated-classes.R
NULL

setAs("HyperparametersBatch", "HyperparametersMultiBatch", function(from, to){
  new("HyperparametersMultiBatch",
      k=k(from),
      mu.0=mu.0(from),
      tau2.0=tau2.0(from),
      eta.0=eta.0(from),
      m2.0=m2.0(from),
      alpha=alpha(from),
      beta=from@beta,
      a=a(from),
      b=b(from))
})

setAs("BatchModel", "MultiBatchModel", function(from, to){
  hypp <- as(hyperParams(from), "HyperparametersMultiBatch")
  new("MultiBatchModel",
      k=k(from),
      hyperparams=hypp,
      theta=theta(from),
      batch=batch(from),
      sigma2=sigma2(from),
      nu.0=nu.0(from),
      sigma2.0=sigma2.0(from),
      pi=p(from),
      mu=mu(from),
      tau2=tau2(from),
      data=y(from),
      data.mean=from@data.mean,
      data.prec=from@data.prec,
      z=z(from),
      zfreq=zFreq(from),
      probz=probz(from),
      logprior=logPrior(from),
      loglik=log_lik(from),
      mcmc.chains=chains(from),
      mcmc.params=mcmcParams(from),
      label_switch=label_switch(from),
      marginal_lik=marginal_lik(from),
      .internal.constraint=from@.internal.constraint,
      .internal.counter=from@.internal.counter)
})

setAs("MarginalModel", "SingleBatchModel", function(from, to){
  ##hypp <- as(hyperParams(from), "Hyperparameters")
  hypp <- hyperParams(from)
  new("SingleBatchModel",
      k=k(from),
      hyperparams=hypp,
      theta=theta(from),
      batch=batch(from),
      sigma2=sigma2(from),
      nu.0=nu.0(from),
      sigma2.0=sigma2.0(from),
      pi=p(from),
      mu=mu(from),
      tau2=tau2(from),
      data=y(from),
      data.mean=from@data.mean,
      data.prec=from@data.prec,
      z=z(from),
      zfreq=zFreq(from),
      probz=probz(from),
      logprior=logPrior(from),
      loglik=log_lik(from),
      mcmc.chains=chains(from),
      mcmc.params=mcmcParams(from),
      label_switch=label_switch(from),
      marginal_lik=marginal_lik(from),
      .internal.constraint=from@.internal.constraint,
      .internal.counter=from@.internal.counter)
})


setAs("SingleBatchPooled", "SingleBatchCopyNumber", function(from, to){
  ##hypp <- as(hyperParams(from), "Hyperparameters")
  hypp <- hyperParams(from)
  new("SingleBatchCopyNumber",
      k=k(from),
      hyperparams=hypp,
      theta=theta(from),
      batch=batch(from),
      sigma2=sigma2(from),
      nu.0=nu.0(from),
      sigma2.0=sigma2.0(from),
      pi=p(from),
      mu=mu(from),
      tau2=tau2(from),
      data=y(from),
      data.mean=from@data.mean,
      data.prec=from@data.prec,
      z=z(from),
      zfreq=zFreq(from),
      probz=probz(from),
      logprior=logPrior(from),
      loglik=log_lik(from),
      mcmc.chains=chains(from),
      mcmc.params=mcmcParams(from),
      label_switch=label_switch(from),
      marginal_lik=marginal_lik(from),
      .internal.constraint=from@.internal.constraint,
      .internal.counter=from@.internal.counter,
      mapping=numeric(k(from)))
})

setAs("MultiBatchPooled", "MultiBatchCopyNumberPooled", function(from, to){
  ##hypp <- as(hyperParams(from), "Hyperparameters")
  hypp <- hyperParams(from)
  ub <- unique(batch(from))
  nbatch <- setNames(as.integer(table(batch(from))), ub)
  new("MultiBatchCopyNumberPooled",
      k=k(from),
      hyperparams=hypp,
      theta=theta(from),
      batch=batch(from),
      batchElements=nbatch,
      sigma2=sigma2(from),
      nu.0=nu.0(from),
      sigma2.0=sigma2.0(from),
      pi=p(from),
      mu=mu(from),
      tau2=tau2(from),
      data=y(from),
      data.mean=from@data.mean,
      data.prec=from@data.prec,
      z=z(from),
      zfreq=zFreq(from),
      probz=probz(from),
      logprior=logPrior(from),
      loglik=log_lik(from),
      mcmc.chains=chains(from),
      mcmc.params=mcmcParams(from),
      label_switch=label_switch(from),
      marginal_lik=marginal_lik(from),
      .internal.constraint=from@.internal.constraint,
      .internal.counter=from@.internal.counter,
      mapping=numeric(k(from)))
})
