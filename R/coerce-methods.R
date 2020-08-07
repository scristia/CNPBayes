##setAs("HyperparametersBatch", "HyperparametersMultiBatch", function(from, to){
##  new("HyperparametersMultiBatch",
##      k=k(from),
##      mu.0=mu.0(from),
##      tau2.0=tau2.0(from),
##      eta.0=eta.0(from),
##      m2.0=m2.0(from),
##      alpha=alpha(from),
##      beta=from@beta,
##      a=a(from),
##      b=b(from))
##})
