\name{dppgibbs}
\alias{dppgibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gibbs Sample with Dirichlet process prior
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dppgibbs(r, H, alpha = 1, mu0 = 0, tau20 = 0.1, a = 0.1, b = 0.1, S = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{H}{
%%     ~~Describe \code{H} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{mu0}{
%%     ~~Describe \code{mu0} here~~
}
  \item{tau20}{
%%     ~~Describe \code{tau20} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
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
function (r, H, alpha = 1, mu0 = 0, tau20 = 0.1, a = 0.1, b = 0.1, 
    S = 100) 
{
    if (missing(H)) 
        stop("H missing")
    n <- length(r)
    pi <- ns <- rep(0, H)
    v <- rep(1/H, H)
    v[H] <- 1
    mu <- rep(0, H)
    tau <- sigma2 <- rep(1, H)
    p <- tmp2 <- matrix(0, n, H)
    V <- Mu <- Sigma <- N <- Pi <- matrix(0, S, H)
    C <- matrix(0, S, n)
    grid <- seq(min(y), max(y), length = 500)
    Y <- array(0, dim = c(S, length(grid), H))
    for (i in 1:S) {
        cumv <- cumprod(1 - v)
        pi[1] <- v[1]
        for (h in 2:H) pi[h] <- v[h] * cumv[h - 1]
        for (h in 1:H) tmp2[, h] <- pi[h] * dnorm(y, mu[h], sqrt(sigma2[h]))
        p <- tmp2/apply(tmp2, 1, sum)
        C[i, ] <- c <- rMultinom(p, 1)
        Pi[i, ] <- pi
        for (h in 1:H) ns[h] <- length(c[c == h])
        for (h in 1:(H - 1)) v[h] <- rbeta(1, 1 + ns[h], alpha + 
            sum(ns[(h + 1):H]))
        V[i, ] <- v
        for (h in 1:H) {
            var <- 1/(tau20 + tau[h] * ns[h])
            m <- var * (tau20 * mu0 + tau[h] * sum(y[c == h]))
            Mu[i, h] <- mu[h] <- rnorm(1, m, sqrt(var))
            tau[h] <- rgamma(1, a + ns[h]/2, b + sum((y[c == 
                h] - mu[h])^2)/2)
            Sigma[i, h] <- sigma2[h] <- 1/tau[h]
            Y[i, , h] <- pi[h] * dnorm(grid, mu[h], sqrt(sigma2[h]))
        }
        N[i, ] <- ns
        if (i\%\%100 == 0) 
            print(i)
    }
    list(P = Pi, means = Mu, precs = 1/Sigma, Z = C, N = N)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
