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
