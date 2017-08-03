load("../data/MarginalModelExample.RData")
load("../data/BatchModelExample.RData")

from <- hyperParams(BatchModelExample)


})

from <- BatchModelExample



MultiBatchModelExample <- as(BatchModelExample, "MultiBatchModel")
save(MultiBatchModelExample, file="../data/MultiBatchModelExample.rda")
