library(lattice)
library(MixtureModel)
avgRs <- readRDS("avgRs.rds")
y1 <- avgRs[1,][which(!is.na(avgRs[1,]))]
y2 <- avgRs[215,][which(!is.na(avgRs[215,]))]
y3 <- avgRs[17,][which(!is.na(avgRs[17,]))]
y4 <- avgRs[94,][which(!is.na(avgRs[94,]))]
y5 <- avgRs[95,][which(!is.na(avgRs[95,]))]
len <- length(y1)
df <- data.frame(y=c(y1, y2, y3, y4, y5), levels=rep(1:5, each=len))
df$levels <- factor(df$levels)


## this needs to be automated
post1 <- gibbs.mix(y1, k=2, delta=0.15, tau20=0.1, nu0=1, kappa0=0.1)
mean1 <- colMeans(post1$mean)
prec1 <- colMeans(post1$prec)
pi1 <- colMeans(post1$P)

post2 <- gibbs.mix(y2, k=4, delta=0.15, tau20=0.1, nu0=1, kappa0=0.1)
mean2 <- colMeans(post2$mean)
prec2 <- colMeans(post2$prec)
pi2 <- colMeans(post2$P)

post3 <- gibbs.mix(y3, k=3, delta=0.15, tau20=0.1, nu0=1, kappa0=0.1)
mean3 <- colMeans(post3$mean)
prec3 <- colMeans(post3$prec)
pi3 <- colMeans(post3$P)

post4 <- gibbs.mix(y4, k=3, delta=0.15, tau20=0.1, nu0=1, kappa0=0.1)
mean4 <- colMeans(post4$mean)
prec4 <- colMeans(post4$prec)
pi4 <- colMeans(post4$P)
## make list of density functions
xsim1 <- density(y1)$x
xsim2 <- density(y2)$x
xsim3 <- density(y3)$x
xsim4 <- density(y4)$x
dens <- list("y1" = pi1[1]*dnorm(xsim1, mean1[1], sqrt(1/prec1[1])) +
             pi1[2]*dnorm(xsim1, mean1[2], sqrt(1/prec1[2])),

             "y2" = pi2[1]*dnorm(xsim2, mean2[1], sqrt(1/prec2[1])) +
             pi2[2]*dnorm(xsim2, mean2[2], sqrt(1/prec2[2])) +
             pi2[3]*dnorm(xsim2, mean2[3], sqrt(1/prec2[3])) +
             pi2[4]*dnorm(xsim2, mean2[4], sqrt(1/prec2[4])),

             "y3" = pi3[1]*dnorm(xsim3, mean3[1], sqrt(1/prec3[1])) +
             pi3[2]*dnorm(xsim3, mean3[2], sqrt(1/prec3[2])) +
             pi3[3]*dnorm(xsim3, mean3[3], sqrt(1/prec3[3])),

             "y4" = pi4[1]*dnorm(xsim4, mean4[1], sqrt(1/prec4[1])) +
             pi4[2]*dnorm(xsim4, mean4[2], sqrt(1/prec4[2])) +
             pi4[3]*dnorm(xsim4, mean4[3], sqrt(1/prec4[3])))

layeredPlots <- function(df, col=rgb(0,0.5,1, alpha=0.6), shift=2, ylim) {
    panelfun <- function(x, levels, shift, ...) {
        nshift <- 0
        for(j in seq_along(levels(levels))) {
            d <- density(x[levels==j])
            lpolygon(c(d$x, rev(d$x)),
                     c(rep(nshift, length(d$x)), rev(d$y)+nshift),
                     col=rgb(0.7,0.7,0.7,alpha=0.6), border=rgb(0.7,0.7,0.7))
            #xsim = d$x
            #panel.lines(xsim, dens[[j]] + nshift, lwd=1)
            nshift <- 0.5*max(rev(d$y)) + nshift
        }
        panel.abline(v=0, col="gray50", lty=2)
        panel.abline(h=0, col="black")
    }

    y <- range(density(df$y)$y)
    levels <- df$levels

    histogram(~y, df, levels=levels, shift=shift, layout=c(1,1), breaks=1000,
              border="lightblue", col="lightblue", strip=FALSE, panel=panelfun,
              ylim=ylim, ylab="", par.settings=list(axis.line=list(col = 0)),
              scales=list(y=list(at=NULL)), xlab="Average log Rs")
}
pdf("layeredhists.pdf")
layeredPlots(df, shift=2, ylim=c(0,7))
dev.off()
