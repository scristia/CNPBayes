\name{plotPosts}
\alias{plotPosts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot posterior density
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotPosts(r, posts, burnin = 1:100, main = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{posts}{
%%     ~~Describe \code{posts} here~~
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
function (r, posts, burnin = 1:100, main = "") 
{
    if (length(posts) > 1) {
        logliks <- sapply(posts, function(x) x$loglik)
        K <- sapply(posts, function(x) x$K)
        bics <- sapply(posts, function(x) x$bic)
        dens.comp <- function(y, comp, p, mix.mean, mix.prec) {
            p[comp] * dnorm(y, mean = mix.mean[comp], sd = 1/sqrt(mix.prec[comp]))
        }
        y <- seq(min(r), max(r), len = 1000)
        best.bic <- bics[order(bics)[1:2]]
        K2 <- K[order(bics)[1:2]]
        subs <- K2 - min(K) + 1
        if (K2[1] == 1) {
            p.1 <- mean(posts[[subs[1]]][["P"]][-burnin, ])
            mean.1 <- mean(posts[[subs[1]]][["means"]][-burnin, 
                ])
            prec.1 <- mean(posts[[subs[1]]][["precs"]][-burnin, 
                ])
        }
        else {
            p.1 <- apply(posts[[subs[1]]][["P"]][-burnin, ], 
                2, mean)
            mean.1 <- apply(posts[[subs[1]]][["means"]][-burnin, 
                ], 2, mean)
            prec.1 <- apply(posts[[subs[1]]][["precs"]][-burnin, 
                ], 2, mean)
        }
        if (K2[2] == 1) {
            p.2 <- mean(posts[[subs[2]]][["P"]][-burnin, ])
            mean.2 <- mean(posts[[subs[2]]][["means"]][-burnin, 
                ])
            prec.2 <- mean(posts[[subs[2]]][["precs"]][-burnin, 
                ])
        }
        else {
            p.2 <- apply(posts[[subs[2]]][["P"]][-burnin, ], 
                2, mean)
            mean.2 <- apply(posts[[subs[2]]][["means"]][-burnin, 
                ], 2, mean)
            prec.2 <- apply(posts[[subs[2]]][["precs"]][-burnin, 
                ], 2, mean)
        }
        d.1 <- rowSums(sapply(1:length(p.1), dens.comp, y = y, 
            p = p.1, mix.mean = mean.1, mix.prec = prec.1), na.rm = TRUE)
        d.2 <- rowSums(sapply(1:length(p.2), dens.comp, y = y, 
            p = p.2, mix.mean = mean.2, mix.prec = prec.2), na.rm = TRUE)
        hist(r, breaks = 200, col = "lightgray", border = "lightgray", 
            freq = FALSE, main = main)
        lines(y, d.1, col = "darkgreen", lwd = 2)
        lines(y, d.2, col = rgb(24, 167, 181, max = 255), lty = 2, 
            lwd = 2)
        legend("topleft", legend = paste(paste("K = ", K2, " BIC = ", 
            round(best.bic))), col = c("darkgreen", rgb(24, 167, 
            181, max = 255)), lty = 1:2, bty = "n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
