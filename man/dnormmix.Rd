\name{dnormmix}
\alias{dnormmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find density for normal mixture
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dnormmix(r, mixture, K, burnin = 1:50, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{mixture}{
%%     ~~Describe \code{mixture} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
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
function (r, mixture, K, burnin = 1:50, log = FALSE) 
{
    if (K == 1) {
        mix.mean <- mean(mixture[["means"]][-burnin, ])
        mix.prec <- mean(mixture[["precs"]][-burnin, ])
        p <- mean(mixture[["P"]][-burnin, ])
    }
    else {
        mix.mean <- apply(mixture[["means"]][-burnin, ], 2, mean)
        mix.prec <- apply(mixture[["precs"]][-burnin, ], 2, mean)
        p <- apply(mixture[["P"]][-burnin, ], 2, mean)
    }
    lik.comp <- function(r, comp) {
        p[comp] * dnorm(r, mean = mix.mean[comp], sd = 1/sqrt(mix.prec[comp]))
    }
    liks <- sapply(1:K, lik.comp, r = r)
    d <- rowSums(liks, na.rm = TRUE)
    if (log) {
        d <- log(d)
    }
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
