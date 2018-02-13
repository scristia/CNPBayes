load("../data/MarginalModelExample.RData")
SingleBatchModelExample <- as(MarginalModelExample, "SingleBatchModel")
save(SingleBatchModelExample, file="../data/SingleBatchModelExample.rda")

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
