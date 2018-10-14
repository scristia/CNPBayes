load("../data/MarginalModelExample.RData")
SingleBatchModelExample <- as(MarginalModelExample, "SingleBatchModel")
save(SingleBatchModelExample, file="../data/SingleBatchModelExample.rda")

setMethod("coerce", "MultiBatchModel", "MultiBatchModel2", function(from, to){
  dat <- tibble(id=as.character(seq_along(y(from))),
                oned=y(from),
                batch=batch(from))
  
})




load("../data/BatchModelExample.RData")
MultiBatchModelExample <- as(BatchModelExample, "MultiBatchModel")
save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")

## updating these examples after migration to heavy-tails
model <- updateObject(MultiBatchPooledExample)
hp <- hyperParams(model)
hp@dfr <- 100
hyperParams(model) <- hp
model@u <- rchisq(length(y(model)), hp@dfr)
MultiBatchPooledExample <- model
save(MultiBatchPooledExample, file="../data/MultiBatchPooledExample.rda")

model <- updateObject(MultiBatchModelExample)
hp <- hyperParams(model)
hp@dfr <- 100
hyperParams(model) <- hp
model@u <- rchisq(length(y(model)), hp@dfr)
MultiBatchModelExample <- model
save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")


model <- updateObject(SingleBatchModelExample)
hp <- hyperParams(model)
hp@dfr <- 100
hyperParams(model) <- hp
model@u <- rchisq(length(y(model)), hp@dfr)
SingleBatchModelExample <- model
save(SingleBatchModelExample, file="../data/SingleBatchModelExample.rda")



data(SingleBatchModelExample)
SingleBatchModelExample <- updateObject(SingleBatchModelExample)

set.seed(123)
sb <- simulateBatchData(N=2500, p=rep(1/3, 3),
                        theta=matrix(c(-1, 0, 1), nrow=1),
                        sds=matrix(rep(0.1, 3), nrow=1),
                        batch=rep(1L, 2500),
                        df=100)
SingleBatchModelExample <- sb
save(SingleBatchModelExample, file="../data/SingleBatchModelExample.rda")

data(MultiBatchModelExample)
MultiBatchModelExample <- updateObject(MultiBatchModelExample)
save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")

data(MultiBatchPooledExample)
MultiBatchPooledExample <- updateObject(MultiBatchPooledExample)
save(MultiBatchPooledExample, file="../data/MultiBatchPooledExample.rda")
