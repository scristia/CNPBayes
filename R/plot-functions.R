meltMultiBatchChains <- function(model){
  ch <- chains(model)
  th <- as.data.frame(theta(ch))
  th$param <- "theta"
  th$iter <- factor(1:nrow(th))
  th.m <- melt(th)
  ##
  ## 
  ##
  K <- k(model)
  B <- nBatch(model)
  btch <- matrix(uniqueBatch(model), B, K, byrow=FALSE)
  btch <- rep(as.character(btch), each=iter(model))
  comp <- matrix(1:K, B, K, byrow=TRUE)
  comp <- rep(as.numeric(comp), each=iter(model))
  th.m$batch <- factor(btch)
  th.m$comp <- factor(paste0("k=", comp))
  th.m$iter <- as.integer(th.m$iter)

  s <- as.data.frame(sigma(ch))
  s$iter <- th$iter
  s$param <- "sigma"
  s.m <- melt(s)
  s.m$iter <- th.m$iter
  s.m$batch <- th.m$batch
  s.m$comp <- th.m$comp

  nu0 <- as.data.frame(nu.0(ch))
  nu0$iter <- th$iter
  nu0$param <- "nu0"
  nu0.m <- melt(nu0)

  s20 <- as.data.frame(sigma2.0(ch))
  s20$iter <- th$iter
  s20$param <- "s20"
  s20.m <- melt(s20)

  mus <- as.data.frame(mu(ch))
  mus$iter <- th$iter
  mus$param <- "mu"
  mus.m <- melt(mus)
  mu.comp <- rep(seq_len(k(model)), each=iter(model))
  mus.m$comp <- factor(paste0("k=", mu.comp))

  prob <- as.data.frame(p(ch))
  prob$iter <- th$iter
  prob$param <- "p"
  prob.m <- melt(prob)
  prob.m$comp <- mus.m$comp

  taus <- as.data.frame(tau(ch))
  taus$iter <- th$iter
  taus$param <- "tau"
  taus.m <- melt(taus)
  ##taus.m$comp <- factor(rep(seq_len(k(model)), each=iter(model)))
  taus.m$comp <- mus.m$comp
  prob.m$iter <- taus.m$iter <- mus.m$iter <- as.integer(mus.m$iter)

  dat.batch <- rbind(th.m,
                     s.m)
  dat.comp <- rbind(mus.m, taus.m, prob.m)
  dat <- rbind(nu0.m,
               s20.m)
  dat$iter <- as.integer(dat$iter)
  list(batch=dat.batch,
       comp=dat.comp,
       single=dat)
}


meltSingleBatchChains <- function(model){
  ch <- chains(model)
  th <- as.data.frame(theta(ch))
  th$param <- "theta"
  th$iter <- factor(1:nrow(th))
  th.m <- melt(th)
  ##
  ## 
  ##
  K <- k(model)
  comp <- rep(seq_len(K), each=nrow(th))
  th.m$comp <- factor(paste0("k=", comp))
  th.m$iter <- as.integer(th.m$iter)

  s <- as.data.frame(sigma(ch))
  s$iter <- th$iter
  s$param <- "sigma"
  s.m <- melt(s)
  s.m$iter <- th.m$iter
  s.m$comp <- th.m$comp

  nu0 <- as.data.frame(nu.0(ch))
  nu0$iter <- th$iter
  nu0$param <- "nu0"
  nu0.m <- melt(nu0)

  s20 <- as.data.frame(sigma2.0(ch))
  s20$iter <- th$iter
  s20$param <- "s20"
  s20.m <- melt(s20)

  mus <- as.data.frame(mu(ch))
  mus$iter <- th$iter
  mus$param <- "mu"
  mus.m <- melt(mus)
  ##mus.m$comp <- factor(rep(seq_len(k(model)), each=iter(model)))

  prob <- as.data.frame(p(ch))
  prob$iter <- th$iter
  prob$param <- "p"
  prob.m <- melt(prob)
  p.comp <- rep(seq_len(k(model)), each=iter(model))
  prob.m$comp <- factor(paste0("k=", p.comp))

  taus <- as.data.frame(tau(ch))
  taus$iter <- th$iter
  taus$param <- "tau"
  taus.m <- melt(taus)
  ##taus.m$comp <- factor(rep(seq_len(k(model)), each=iter(model)))
  ##prob.m$iter <- taus.m$iter <- mus.m$iter <- as.integer(mus.m$iter)
  taus.m$iter <- as.integer(mus.m$iter)
  dat.comp <- rbind(th.m,
                    s.m,
                    prob.m)
  dat.comp$iter <- as.integer(dat.comp$iter)
  dat <- rbind(nu0.m,
               s20.m,
               mu=mus.m,
               tau=taus.m)
  dat$iter <- as.integer(dat$iter)
  list(comp=dat.comp,
       single=dat)
}

ggSingleBatchChains <- function(model){
  melt.ch <- meltSingleBatchChains(model)
  dat.comp <- melt.ch[["comp"]]     ## component-specific
  dat.single <- melt.ch[["single"]] ## single-parameter
  iter <- value <- comp <- param <- NULL
  p.comp <- ggplot(dat.comp, aes(iter, value, group=comp)) +
    geom_point(size=0.3, aes(color=comp)) +
    geom_line(aes(color=comp)) +
    facet_wrap(~param, scales="free_y")
  p.single <- ggplot(dat.single, aes(iter, value, group="param")) +
    geom_point(size=0.3, color="gray") +
    geom_line(color="gray") +
    facet_wrap(~param, scales="free_y")
   list(comp=p.comp, single=p.single)
}

ggMultiBatchChains <- function(model){
  melt.ch <- meltMultiBatchChains(model)
  dat.batch <- melt.ch$batch
  dat.comp <- melt.ch$comp
  dat.single <- melt.ch$single
  iter <- value <- batch <- param <- comp <- NULL
  p.batch <- ggplot(dat.batch, aes(iter, value, group=batch)) +
    geom_point(size=0.3, aes(color=batch)) +
    geom_line(aes(color=batch)) +
    facet_grid(param ~ comp, scales="free_y")

  p.comp <- ggplot(dat.comp, aes(iter, value, group=comp)) +
    geom_point(size=0.3, aes(color=comp)) +
    geom_line(aes(color=comp)) +
    facet_wrap(~param, scales="free_y")

  p.single <- ggplot(dat.single, aes(iter, value, group="param")) +
    geom_point(size=0.3, color="gray") +
    geom_line(color="gray") +
    facet_wrap(~param, scales="free_y")
  list(batch=p.batch, comp=p.comp, single=p.single)
}

ggSingleBatch <- function(model){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  df <- singleBatchDensities(model)
  df.observed <- data.frame(y=observed(model))
  ..density.. <- name <- x <- d <- NULL 
  ggplot(df, aes(x, d)) +
    geom_histogram(data=df.observed,
                   aes(y, ..density..),
                   bins=300, inherit.aes=FALSE) +
    geom_area(stat="identity", aes(color=name, fill=name),
              alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(""), color=guide_legend("")) 
}

ggMultiBatch <- function(model){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  df <- multiBatchDensities(model)
  nb <- nBatch(model)
  df.observed <- data.frame(y=observed(model),
                            batch=batch(model))
  y <- ..density.. <- name <- x <- d <- NULL
  ggplot(df, aes(x, d)) +
    geom_histogram(data=df.observed,
                   aes(y, ..density..),
                   bins=300, inherit.aes=FALSE) +
    geom_area(stat="identity", aes(color=name, fill=name),
              alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(""), color=guide_legend("")) +
    facet_wrap(~batch, nrow=nb)
}

#' Create a data.frame of the component densities for each batch
#'
#' @param object an object of class \code{BatchModel}
#' @return a \code{{data.frame}}
#' @export
#' @examples
#'    nbatch <- 3
#'    k <- 3
#'    means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
#'        0, 0.05), nbatch, k, byrow = FALSE)
#'    sds <- matrix(0.15, nbatch, k)
#'    sds[, 1] <- 0.3
#'    N <- 1000
#'    truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
#'                                                  length.out = N),
#'                               p = c(1/10, 1/5, 1 - 0.1 - 0.2), theta = means,
#'                               sds = sds)
#'    mcmcp <- McmcParams(iter = 1000, burnin = 500, thin = 1,
#'                        nStarts = 10)
#'
#'    ## this parameter setting for m2.0 allows a lot of varation of the thetas
#'    ## between batch
#'    hypp <- CNPBayes:::HyperparametersBatch(m2.0 = 1/60, eta.0 = 1800,
#'                                            k = 3, a = 1/6, b = 180)
#'    model <- BatchModel(data = y(truth), batch = batch(truth),
#'                        k = 3, mcmc.params = mcmcp, hypp = hypp)
#'    model <- posteriorSimulation(model)
#'    df <- multiBatchDensities(model)
#'    df.observed <- data.frame(y=y(model), batch=batch(model))
#'    library(ggplot2)
#'    ggplot(df, aes(x, d)) +
#'    geom_histogram(data=df.observed,
#'                   aes(y, ..density..),
#'                   bins=300, inherit.aes=FALSE) +
#'    geom_area(stat="identity", aes(color=name, fill=name),
#'              alpha=0.4) +
#'    xlab("quantiles") + ylab("density") +
#'    guides(fill=guide_legend(""), color=guide_legend("")) +
#'    facet_wrap(~batch, nrow=2)
multiBatchDensities <- function(object){
  probs <- p(object)
  thetas <- theta(object)
  sigmas <- sigma(object)
  P <- matrix(probs, nrow(thetas), ncol(thetas), byrow=TRUE)
  rownames(P) <- uniqueBatch(object)
  avglrrs <- observed(object)
  quantiles <- seq(min(avglrrs), max(avglrrs), length.out=500)
  batchPr <- table(batch(object))/length(y(object))
  dens.list <- batchDensities(quantiles, uniqueBatch(object),
                              thetas, sigmas, P, batchPr)
  ##component <- lapply(dens.list, rowSums)
  ##overall <- rowSums(do.call(cbind, component))
  ix <- order(thetas[1, ])
  dens.list <- dens.list[ix]
  d <- do.call(rbind, dens.list)
  K <- ncol(thetas)
  NB <- nBatch(object)
  over <- Reduce("+", dens.list)
  batches.overall <- rep(uniqueBatch(object), each=nrow(over))
  quantile.overall <- rep(quantiles, nBatch(object))
  overall <- as.numeric(over)

  d.vec <- as.numeric(d)
  d.vec <- c(d.vec, overall)
  batches <- c(rep(uniqueBatch(object), each=nrow(d)),
               batches.overall)
  K <- seq_len(ncol(thetas))
  name <- paste0("cn", K-1)
  name <- rep(rep(name, elementNROWS(dens.list)), nBatch(object))
  name <- c(name, rep("overall", length(overall)))
  x <- rep(rep(quantiles, length(dens.list)), nBatch(object))
  x <- c(x, quantile.overall)
  df <- data.frame(x=x, d=d.vec, name=name, batch=batches)
  df$batch <- factor(df$batch, uniqueBatch(object))
  df$name <- factor(df$name, levels=c("overall", paste0("cn", K-1)))
  df
}
