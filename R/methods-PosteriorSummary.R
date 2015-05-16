PosteriorSummary <- function(p_theta=matrix(), chib=numeric(), berkhof=numeric(),
                             marginal=numeric(),
                             delta_marginal=numeric()){
  new("PosteriorSummary", p_theta=p_theta,
      chib=chib,
      berkhof=berkhof, marginal=marginal,
      delta_marginal=delta_marginal)
}

p_theta <- function(object) object@p_theta

chib <- function(object) object@chib

berkhof <- function(object) object@berkhof

setMethod("marginal", "PosteriorSummary", function(object) object@marginal)

deltaMarginal <- function(object) object@delta_marginal

setMethod("show", "PosteriorSummary", function(object){
  cat("Object of class 'PosteriorSummary'\n")
  cat("   Estimates of p(theta* | y ):\n")
  cat("     - Chib's :", round(chib(object), 2),    "\n")
  cat("     - Berkhof:", round(berkhof(object), 2), "\n")
  cat("   Marginal likelihood:", round(marginal(object), 2), "\n")
  cat("   see chib(), berkhof(), p_theta(), marginal()\n")
})



PosteriorSummaryList <- function(data=list(), names=character()){
  new("PosteriorSummaryList", data=data, names=names)
}

setMethod("[", "PosteriorSummaryList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@data <- x@data[i]
    x@names <- x@names[i]
  }
  x
})

setMethod("[[", "PosteriorSummaryList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x <- x@data[[i]]
  }
  x
})

setMethod("show", "PosteriorSummaryList", function(object){
  cat("Object of class 'PosteriorSummaryList'\n")

})

setMethod("summary", "PosteriorSummaryList", function(object, ...){
  chibs <- round(sapply(object, chib), 2)
  berk <- round(sapply(object, berkhof), 2)
  my <- round(sapply(object, marginal), 2)
  d <- round(sapply(object, deltaMarginal), 2)
  xx <- cbind(chibs, berk, my, d)
  rownames(xx) <- names(object)
  colnames(xx) <- c("chib", "berkhof", "marginal", "spread(marginal)")
  cat(xx, "\n")
  xx
})

setMethod("length", "PosteriorSummaryList", function(x) length(x@data))

setMethod("lapply", "PosteriorSummaryList", function(X, FUN, ...){
  results <- list()
  for(i in seq_along(X)){
    results[[i]] <- FUN(X[[i]], ...)
  }
  results
})

setMethod("sapply", "PosteriorSummaryList", function(X, FUN, ...){
  results <- lapply(X, FUN, ...)
  as.numeric(results)
})
