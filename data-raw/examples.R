load("../data/MarginalModelExample.RData")
SingleBatchModelExample <- as(MarginalModelExample, "SingleBatchModel")
save(SingleBatchModelExample, file="../data/SingleBatchModelExample.rda")

load("../data/BatchModelExample.RData")
MultiBatchModelExample <- as(BatchModelExample, "MultiBatchModel")
save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")
