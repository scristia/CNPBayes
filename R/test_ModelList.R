test_ModelList <- function(){
  m <- MarginalModelList()
  checkTrue(validObject(m))
  checkTrue(is.numeric(maxperm(m)))


  mode_indices <- CNPBayes:::trackModes(1:3, maxperm=3)
  index <- mode_indices[[1]]
  checkEquals(index[[3]],
              matrix(c(3, 1, 2, 2, 3, 1, 1, 3, 2), nrow=3, ncol=3, byrow=TRUE))
}
