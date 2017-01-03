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

#' ggplot-style functions for diagnosing convergence
#'
#' @return a list of \code{ggplot} objects. Chains are grouped by the length of
#'   the parameter vector. For example, in the single-batch model, the means
#'   (theta) and variances (sigma2) are component-specific (length k, where k is
#'   number of components) and are plotted together in a single \code{ggplot}
#'   object.
#'
#' @examples
#' plist.sb <- ggSingleBatchChains(MarginalModelExample)
#' ## chains for parameter vectors of length k
#' plist.sb[["comp"]]
#' ## chains for parameters vectors of length 1
#' plist.sb[["single"]]
#' plist.mb <- ggMultiBatchChains(BatchModelExample)
#' ## chains for parameters that are batch- and component-specific
#' plist.mb[["batch"]]
#' ## chains for parameters vectors of length k
#' plist.mb[["comp"]]
#' ## chains for parameter vectors of length 1
#' plist.mb[["single"]]
#' @export
#' @rdname ggplot-functions
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


#' @export
#' @rdname ggplot-functions
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

#' ggplot wrapper for plotting the data at a single CNP and the model-based densities
#'
#' @param model a \code{BatchModel} or  \code{MarginalModel} object
#' @examples
#' ggMultiBatch(BatchModelExample)
#' ggSingleBatch(MarginalModelExample)
#' @export
#' @return a \code{ggplot} object
#' @rdname ggplot-functions
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


#' @export
#' @rdname ggplot-functions
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


#' @export
#' @examples
#' df <- multiBatchDensities(BatchModelExample)
#' head(df)
#' @rdname ggplot-functions
multiBatchDensities <- function(model){
  probs <- p(model)
  thetas <- theta(model)
  sigmas <- sigma(model)
  P <- matrix(probs, nrow(thetas), ncol(thetas), byrow=TRUE)
  rownames(P) <- uniqueBatch(model)
  avglrrs <- observed(model)
  quantiles <- seq(min(avglrrs), max(avglrrs), length.out=500)
  batchPr <- table(batch(model))/length(y(model))
  dens.list <- batchDensities(quantiles, uniqueBatch(model),
                              thetas, sigmas, P, batchPr)
  ##component <- lapply(dens.list, rowSums)
  ##overall <- rowSums(do.call(cbind, component))
  ix <- order(thetas[1, ])
  dens.list <- dens.list[ix]
  d <- do.call(rbind, dens.list)
  K <- ncol(thetas)
  NB <- nBatch(model)
  over <- Reduce("+", dens.list)
  batches.overall <- rep(uniqueBatch(model), each=nrow(over))
  quantile.overall <- rep(quantiles, nBatch(model))
  overall <- as.numeric(over)

  d.vec <- as.numeric(d)
  d.vec <- c(d.vec, overall)
  batches <- c(rep(uniqueBatch(model), each=nrow(d)),
               batches.overall)
  K <- seq_len(ncol(thetas))
  name <- paste0("cn", K-1)
  name <- rep(rep(name, elementNROWS(dens.list)), nBatch(model))
  name <- c(name, rep("overall", length(overall)))
  x <- rep(rep(quantiles, length(dens.list)), nBatch(model))
  x <- c(x, quantile.overall)
  df <- data.frame(x=x, d=d.vec, name=name, batch=batches)
  df$batch <- factor(df$batch, uniqueBatch(model))
  df$name <- factor(df$name, levels=c("overall", paste0("cn", K-1)))
  df
}
