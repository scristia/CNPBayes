## hack of code found on stackexchange, which I can not find again
min.f1f2 <- function(x, mu1, mu2, sd1, sd2, pi1, pi2) {
    f1 <- pi1*dnorm(x, mean=mu1, sd=sd1)
    f2 <- pi2*dnorm(x, mean=mu2, sd=sd2)
    pmin(f1, f2)
}

avgRs <- readRDS("~/Software/CNPBayes/CNPBayes/data/avgRs_wc_ea.rds")
r <- avgRs[,301]
x <- readRDS("posteriors_EA.RDS")
xr <- x[[301]]

mu1 <- xr$MU[2]
mu2 <- xr$MU[3]
sd1 <- sqrt(xr$SIGMA[2])
sd2 <- sqrt(xr$SIGMA[3])
pi1 <- xr$PI[1]
pi2 <- xr$PI[2]
xs <- seq(min(mu1 - 3*sd1, mu2 - 3*sd2), max(mu1 + 3*sd1, mu2 + 3*sd2), .01)
f1 <- pi1*dnorm(xs, mean=mu1, sd=sd1)
f2 <- pi2*dnorm(xs, mean=mu2, sd=sd2)

ys <- min.f1f2(xs, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2, pi1=pi1, pi2=pi2)
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])

hist(r, breaks=200, col="lightgray", border="lightgray", freq=FALSE)
K <- xr$K
y <- seq(min(r), max(r), len=1000)
for(k in 1:K) lines(y, xr$PI[k]*dnorm(y, xr$MU[k], sqrt(xr$SIGMA[k])),
                    col='gray40', lwd=2)
polygon(xs, ys, col="gray40", density=50, angle=45)

############ ignore below for now
plot(xs, f1, type="l", ylim=c(0, max(f1,f2)), ylab="density")
lines(xs, f2, lty="dotted")
ys <- min.f1f2(xs, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2, pi1=pi1, pi2=pi2)
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])
polygon(xs, ys, col="gray")

### only works for sd1 = sd2
#SMD <- (mu1-mu2)/sd1
#2 * pnorm(-abs(SMD)/2)

### this works in general
integrate(min.f1f2, -Inf, Inf, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2, pi1=pi1, pi2=pi2)

