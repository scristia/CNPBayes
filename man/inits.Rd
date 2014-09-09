\name{inits}
\alias{inits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find initial values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inits(r, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (r, K) 
{
    if (K == 1) 
        return(list(mu0 = mean(r), sigma20 = var(r), nn = length(r)))
    hdel <- min(r) < -0.9 & (length(r[r < -0.75]) > 1)
    if (!hdel | length(r[r < -0.75]) > 1000) {
        pars <- kmeans(r, centers = K, nstart = 50)
        mu0 <- sort(pars$centers)
        s20 <- pars$withinss/(pars$size - 1)[order(pars$centers)]
        nn <- pars$size[order(pars$centers)]
        return(list(mu0 = mu0, sigma20 = s20, nn = nn))
    }
    else {
        mu1 <- median(r[r < -0.75])
        s201 <- 0.25^2
        pars <- kmeans(r[r > -0.75], centers = K - 1, nstart = 20)
        nn2 <- pars$size[order(pars$centers)]
        mu0 <- c(mu1, sort(pars$centers))
        s20 <- c(s201, pars$withinss/(pars$size - 1)[order(pars$centers)])
        nn <- c(length(r) - sum(nn2), nn2)
        return(list(mu0 = mu0, sigma20 = s20, nn = nn))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
