\name{gibbs}
\alias{gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gibbs Sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gibbs(r, nn, alpha = rep(1, 3), mu0 = c(-4, -0.5, 0), tau20 = 0.1, nu0 = 1, kappa0 = 1, sigma20 = 0.1, K = 3, S = 100, delta = 0.15, burnin = 1:100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{nn}{
%%     ~~Describe \code{nn} here~~
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
  \item{nu0}{
%%     ~~Describe \code{nu0} here~~
}
  \item{kappa0}{
%%     ~~Describe \code{kappa0} here~~
}
  \item{sigma20}{
%%     ~~Describe \code{sigma20} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
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
function (r, nn, alpha = rep(1, 3), mu0 = c(-4, -0.5, 0), tau20 = 0.1, 
    nu0 = 1, kappa0 = 1, sigma20 = 0.1, K = 3, S = 100, delta = 0.15, 
    burnin = 1:100) 
{
    is.1ornull <- function(x) {
        (length(x) == 1) | is.null(x)
    }
    if (missing(mu0)) {
        pars <- inits(r, K)
        mu0 <- pars$mu0
        sigma20 <- pars$sigma20
        nn <- pars$nn
        print(mu0)
    }
    if (sum(nn) != length(r)) 
        stop("nn must sum to length(r)")
    z <- rep(NA, length(r))
    s2n <- tau2n <- mun <- numeric(K)
    p <- matrix(NA, sum(nn), K)
    nun <- nu0 + nn
    a0 <- min(r)
    b0 <- max(r)
    s2 <- var(r, na.rm = TRUE)
    precs <- means <- matrix(NA, S, K)
    rbar <- means[1, ] <- rnorm(K, mu0, tau20)
    Z <- matrix(NA, length(r), S - 1)
    precs[1, ] <- 1/sigma20
    PI <- matrix(NA, S, K)
    theta <- c()
    for (s in 2:S) {
        pi <- rdirichlet(1, alpha + nn)
        tau2n <- 1/(1/tau20 + nn * precs[s - 1, ])
        nun <- nu0 + nn
        s2n <- 1/nun * (nu0 * sigma20 + (nn - 1) * s2 + kappa0 * 
            nn/(kappa0 + nn) * (rbar - mu0)^2)
        mun <- (1/tau20)/(1/tau20 + nn * precs[s - 1, ]) * mu0 + 
            nn * precs[s - 1, ]/(1/tau20 + nn * precs[s - 1, 
                ]) * rbar
        theta <- means[s - 1, ]
        fint <- which(is.finite(theta))
        endpoints <- c(a0, theta[fint], b0)
        for (k in 1:length(fint)) {
            if (is.finite(mun[fint[k]])) {
                a <- endpoints[k] + delta
                b <- endpoints[k + 2] - delta
                if (mun[fint[k]] < a) 
                  mun[fint[k]] <- a
                if (mun[fint[k]] > b) 
                  mun[fint[k]] <- b
                theta[fint[k]] <- constr.draw(mun[fint[k]], tau2n[fint[k]], 
                  a, b)
            }
        }
        prec <- rgamma(K, nun/2, nun/2 * s2n)
        means[s, ] <- theta
        precs[s, ] <- prec
        d <- matrix(NA, nrow = length(r), ncol = K)
        for (i in 1:K) d[, i] <- pi[i] * dnorm(r, theta[i], sqrt(1/prec[i]))
        p <- d/apply(d, 1, sum)
        u <- runif(length(r))
        tmp <- p[, 1]
        z[u < tmp] <- 0
        if (K > 1) {
            for (i in 2:K) {
                z[tmp < u & u < tmp + p[, i]] <- i - 1
                tmp <- tmp + p[, i]
            }
        }
        for (i in 1:K) nn[i] <- sum(z == (i - 1))
        if (all(nn > 1)) {
            rbar <- sapply(split(r, z), mean, na.rm = TRUE)
            s2 <- sapply(split(r, z), var, na.rm = TRUE)
        }
        else {
            ww <- split(r, z)
            ww <- ww[!sapply(ww, is.1ornull)]
            s2[nn > 0] <- sapply(ww, var, na.rm = TRUE)
            if (any(nn == 1)) {
                s2[nn == 1] <- NA
                if (all(nn > 0)) 
                  rbar <- sapply(split(r, z), mean, na.rm = TRUE)
            }
            if (any(nn == 0)) {
                s2[nn == 0] <- NA
                rbar[nn != 0] <- sapply(split(r, z)[which(nn != 
                  0)], mean, na.rm = TRUE)
                rbar[nn == 0] <- NA
            }
        }
        Z[, s - 1] <- z
        PI[s, ] <- pi
    }
    post <- list(P = PI, means = means, precs = precs, Z = Z, 
        n = nn)
    loglik <- loglik.normmix(r, post, K = k, burnin = burnin)
    bic <- -2 * loglik + (3 * K - 1) * log(length(r))
    c(post, loglik = loglik, bic = bic, K = K)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
