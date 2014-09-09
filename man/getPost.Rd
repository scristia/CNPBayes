\name{getPost}
\alias{getPost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Posterior
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPost(r, kmin = 2, kmax = 5, delta = 0.15, S = 200, plot = F, burnin = 1:100, main = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{kmin}{
%%     ~~Describe \code{kmin} here~~
}
  \item{kmax}{
%%     ~~Describe \code{kmax} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (r, kmin = 2, kmax = 5, delta = 0.15, S = 200, plot = F, 
    burnin = 1:100, main = "") 
{
    r <- r[!is.na(r)]
    n <- length(r)
    loglik <- c()
    bic <- c()
    kpost <- list()
    for (k in kmin:kmax) {
        inits <- inits(r, k)
        mu0 <- inits$mu0
        sigma20 <- inits$sigma20
        nn <- inits$nn
        print(nn)
        alpha <- rep(1, k)
        print(mu0)
        post <- gibbs(r = r, nn = nn, tau20 = 0.1, K = k, alpha = alpha, 
            sigma20 = sigma20, mu0 = mu0, S = S, delta = delta)
        loglik[k - kmin + 1] <- post$loglik
        bic[k - kmin + 1] <- post$bic
        kpost[[k - kmin + 1]] <- post
    }
    results <- data.frame(K = seq(kmin, kmax), `log-likelihood` = loglik, 
        BIC = bic)
    print(results)
    if (plot) {
        plotPosts(r = r, posts = kpost, burnin = burnin, main = main)
    }
    return(kpost)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
