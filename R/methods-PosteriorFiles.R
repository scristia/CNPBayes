#' @include AllClasses.R

PosteriorFiles <- function(model=character(), post1=character(), ##post2=character(), post3=character(),
                           isMarginalModel=logical()){
  new("PosteriorFiles", model=model, post1=post1,##, post2=post2, post3=post3,
      isMarginalModel=isMarginalModel)
}

setMethod("isMarginalModel", "PosteriorFiles", function(object) object@isMarginalModel)
setMethod("model", "PosteriorFiles", function(object) object@model)
setMethod("post1", "PosteriorFiles", function(object) object@post1)
##setMethod("post2", "PosteriorFiles", function(object) object@post2)
##setMethod("post3", "PosteriorFiles", function(object) object@post3)



setMethod("postFiles", "PosteriorFiles", function(object) post1(object))


setMethod("[", "PosteriorFiles", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@model <- model(x)[i]
    x@post1 <- post1(x)[i]
##    x@post2 <- post2(x)[i]
##    x@post3 <- post3(x)[i]
  }
  x
})

setMethod("length", "PosteriorFiles", function(x) length(model(x)))

setMethod("show", "PosteriorFiles", function(object){
  cat("An object of class 'PosteriorFiles'\n")
  cat("   number of models (CNP loci):", length(object), "\n")
  type <- ifelse(isMarginalModel(object), "marginal", "batch")
  cat("   type of model (marginal|batch):", type, "\n")
  cat("   See model(), postFiles()\n")
})

#' @export
getFiles <- function(outdir, cnpids, model){
  model.files <- modelFiles(outdir, cnpids)
  post.files <- postSummaryFiles(outdir, cnpids)
  if(model=="batch"){
    batch <- PosteriorFiles(model=model.files[["batch"]],
                            post1=post.files$batch[, 1],
                            ##post2=post.files$batch[, 2],
                            ##post3=post.files$batch[, 3],
                            isMarginalModel=FALSE)
    return(batch)
  }
  if(model=="marginal"){
    marginal <- PosteriorFiles(model=model.files[["marginal"]],
                               post1=post.files$marginal[, 1],
                               ##post2=post.files$marginal[, 2],
                               ##post3=post.files$marginal[, 3],
                               isMarginalModel=TRUE)
    return(marginal)
  }
}

.getPosteriorStatsBatch <- function(xlist){
  k1 <- xlist[[1]][, "B1", drop=FALSE]
  k2 <- cbind(xlist[[1]][, "B2"], xlist[[2]][, "B2"])
  k3 <- do.call(cbind, lapply(xlist, "[", , "B3"))
  k4 <- do.call(cbind, lapply(xlist, "[", , "B4"))
  list(k1=k1,
       k2=k2,
       k3=k3,
       k4=k4)
}

.getPosteriorStatsMarginal <- function(xlist){
  k1 <- xlist[[1]][, "M1", drop=FALSE]
  k2 <- cbind(xlist[[1]][, "M2"], xlist[[2]][, "M2"])
  k3 <- do.call(cbind, lapply(xlist, "[", , "M3"))
  k4 <- do.call(cbind, lapply(xlist, "[", , "M4"))
  list(k1=k1,
       k2=k2,
       k3=k3,
       k4=k4)
}

#' @export
topTwo <- function(m.y){
  sort(m.y[is.finite(m.y)], decreasing=TRUE)[1:2]
}

#' @export
bayesFactor <- function(m.y){
  m.y <- topTwo(m.y)
  setNames(abs(diff(m.y)), paste(names(m.y), collapse="-"))
}

modelAbbrv <- function(x) ifelse(isMarginalModel(x), "M", "B")

#' @export
getPosteriorStats <- function(object){
  if(length(object) > 1) {
    message("Only extracting posterior summaries for first CNP")
    object <- object[1]
  }
  x <- readRDS(post1(object))
  xx <- do.call(cbind, lapply(x, rowMeans))
  post.range <- unlist(lapply(x, posteriorRange))
  marginals <- setNames(computeMarginal(xx),
                         paste0(modelAbbrv(object), 1:4))
  marginals
}

summarizePost <- function(x, thr=5){
  results <- foreach(i=seq_along(x), .combine="cbind") %do% {
    ptheta <- colSums(x[[i]][c("theta", "sigma2", "pmix"), , drop=FALSE])
    d <- diff(range(ptheta, na.rm=TRUE))
    if(is.nan(d)) d <- Inf
    result <- rowMeans(x[[i]], na.rm=TRUE)
    if(d > thr){
      result[c("theta", "sigma2", "pmix")] <- Inf
    }
    result
  }
  colnames(results) <- names(x)
  results
}

#' @export
computeMarginal <- function(stats){
  colSums(stats[c("prior", "loglik"), , drop=FALSE]) - colSums(stats[c("theta", "sigma2", "pmix"), , drop=FALSE])
}
