##load("../data/MarginalModelExample.RData")
##SingleBatchModelExample <- as(MarginalModelExample, "SingleBatchModel")
##save(SingleBatchModelExample, file="../data/SingleBatchModelExample.rda")

##load("../data/BatchModelExample.RData")
##MultiBatchModelExample <- as(BatchModelExample, "MultiBatchModel")
##save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")

## updating these examples after migration to heavy-tails
##model <- updateObject(MultiBatchPooledExample)
##hp <- hyperParams(model)
##hp@dfr <- 100
##hyperParams(model) <- hp
##model@u <- rchisq(length(y(model)), hp@dfr)
##MultiBatchPooledExample <- model
##save(MultiBatchPooledExample, file="../data/MultiBatchPooledExample.rda")

#model <- updateObject(MultiBatchModelExample)
#hp <- hyperParams(model)
#hp@dfr <- 100
#hyperParams(model) <- hp
#model@u <- rchisq(length(y(model)), hp@dfr)
#MultiBatchModelExample <- model
                                        #save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")

load("../data/MultiBatchModelExample.rda")
MultiBatchModelExample <- updateObject(MultiBatchModelExample)
MultiBatchModelExample <- posteriorSimulation(MultiBatchModelExample)
save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda",
     compression_level=9)

MultiBatchPooledExample <- updateObject(MultiBatchPooledExample)
MultiBatchPooledExample <- posteriorSimulation(MultiBatchPooledExample)
save(MultiBatchPooledExample, file="../data/MultiBatchPooledExample.rda",
     compression_level=9)
