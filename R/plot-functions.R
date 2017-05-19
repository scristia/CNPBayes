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

singleBatchDensities <- function(object){
  dnorm_poly(object)
}

.range_quantiles <- function(mean, sd){
  qnorm(c(0.001, 0.999), mean=mean, sd=sd)
}

dnorm_quantiles <- function(y, mean, sd){
  x <- list()
  for(i in seq_along(mean)){
    x[[i]] <- .range_quantiles(mean[i], sd[i])
  }
  xx <- unlist(x)
  limits <- c(min(y, min(xx)),
              max(y, max(xx)))
  seq(limits[1], limits[2], length.out=500)
}

.dnorm_poly <- function(x, p, mean, sd){
  y <- p*dnorm(x, mean=mean, sd=sd)
  yy <- c(y, rep(0, length(y)))
  xx <- c(x, rev(x))
  tmp <- data.frame(y=yy,
                    x=xx)
}

dnorm_singlebatch <- function(qtiles, p, mean, sd){
  df.list <- list()
  for(i in seq_along(mean)){
    dat <- .dnorm_poly(qtiles, p[i], mean[i], sd[i])
     if(i == 1){
      overall <- dat$y
    } else{
      overall <- overall + dat$y
    }
    df.list[[i]] <- dat
  }
  df <- do.call(rbind, df.list)
  L <- sapply(df.list, nrow)
  df$component <- factor(rep(seq_along(mean), L))
  df.overall <- data.frame(y=overall, x=dat$x)
  df.overall$component <- "marginal"
  df <- rbind(df, df.overall)
  df$component <- factor(df$component, levels=c("marginal", seq_along(mean)))
  if(length(mean) == 1){
    df <- df[df$component != "overall", ]
    df$component <- factor(df$component)
  }
  df
}

dnorm_poly <- function(model){
  mixprob <- p(model)
  means <- theta(model)
  sds <- sigma(model)
  df.list <- list()
  qtiles <- dnorm_quantiles(y(model), means, sds)
  df <- dnorm_singlebatch(qtiles, mixprob, means, sds)
  df
}

.gg_singlebatch <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  df.observed <- data.frame(y=observed(model))
  ## see stat_function
  ##  ggplot(df, aes(x, d, group=name)) +
  if(missing(bins))
    bins <- nrow(df.observed)/2
  dat <- dnorm_poly(model)
  component <- x <- y <- ..density.. <- NULL
  ggplot(dat, aes(x, y, group=component)) +
    geom_histogram(data=df.observed, aes(y, ..density..),
                   bins=bins,
                   inherit.aes=FALSE) +
    geom_polygon(aes(fill=component, color=component), alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_y_sqrt() +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(""), color=guide_legend(""))
}

#' ggplot wrapper for plotting the data at a single CNP and the model-based densities
#'
#' @param model a \code{BatchModel} or  \code{MarginalModel} object
#' @param bins length-one integer vector indicating the number of bins for the histograms (passed to \code{geom_histogram})
#' @examples
#' ggMultiBatch(BatchModelExample)
#' ggSingleBatch(MarginalModelExample)
#' @export
#' @return a \code{ggplot} object
#' @rdname ggplot-functions
setGeneric("ggSingleBatch", function(model, bins) standardGeneric("ggSingleBatch"))

#' @export
#' @rdname ggplot-functions
setGeneric("ggMultiBatch", function(model, bins) standardGeneric("ggMultiBatch"))

#' @aliases ggSingleBatch,MarginalModel-method
#' @rdname ggplot-functions
setMethod("ggSingleBatch", "MarginalModel", function(model, bins){
  .gg_singlebatch(model, bins)
})

.relabel_component <- function(dat, model){
  comp <- as.character(dat$component)
  index <- comp %in% as.character(seq_len(k(model)))
  comp2 <- comp[index]
  ##cn <- .remap(comp2, mapping(model))
  cn <- mapping(model)[as.integer(comp2)]
  cn <- as.character(cn)
  levs <- unique(cn)
  ##cn <- mapping(model)[z(model)]
  ##comp <- as.character(comp)
  comp[index] <- cn
  comp <- factor(comp, levels=c("marginal", levs))
  comp
}

.gg_singlebatch_copynumber <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  df.observed <- data.frame(y=observed(model))
  ## see stat_function
  ##  ggplot(df, aes(x, d, group=name)) +
  if(missing(bins))
    bins <- nrow(df.observed)/2
  dat <- dnorm_poly(model)
  dat$component <- .relabel_component(dat, model)
  if(FALSE)
    with(dat, table(component, component2))
  component <- x <- y <- ..density.. <- NULL
  ggplot(dat, aes(x, y, group=component)) +
    geom_histogram(data=df.observed, aes(y, ..density..),
                   bins=bins,
                   inherit.aes=FALSE) +
    geom_polygon(aes(fill=component, color=component), alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(""), color=guide_legend(""))
}


#' @aliases ggSingleBatch,SingleBatchCopyNumber-method
#' @rdname ggplot-functions
setMethod("ggSingleBatch", "SingleBatchCopyNumber", function(model, bins){
  .gg_singlebatch_copynumber(model, bins)
})

setMethod("ggMultiBatch", "MultiBatchCopyNumber", function(model, bins){
  .gg_multibatch_copynumber(model, bins)
})



dnorm_poly_multibatch <- function(model){
  mixprob <- p(model)
  means <- theta(model)
  sds <- sigma(model)
  df.list <- list()
  yy <- y(model)
  qtiles <- dnorm_quantiles(y(model),
                            as.numeric(means),
                            as.numeric(sds))
  df.list <- list()
  batch.labels <- unique(batch(model))
  w <- table(batch(model))/length(y(model))
  for(i in seq_len(nrow(means))){
    dat <- dnorm_singlebatch(qtiles, mixprob, means[i, ], sds[i, ])
    ## marginal in 'dat' is the marginal for that batch
    dat$y <- dat$y * w[i]
    dat$batch <- paste0("batch: ", batch.labels[i])
    if(i == 1){
      overall <- dat$y
    } else{
      overall <- overall + dat$y
    }
    df.list[[i]] <- dat
  }
  df <- do.call(rbind, df.list)
  L <- sapply(df.list, nrow)
  ##df$component <- factor(rep(seq_along(means), L))
  df$batch <- factor(df$batch)
  df.overall <- data.frame(y=overall, x=dat$x,
                           component=dat$component,
                           batch="overall")
  df <- rbind(df, df.overall)
  df$batch <- factor(df$batch,
                     levels=c("overall", paste0("batch: ", batch.labels)))
  df
}

batchDensities <- function(x, batches, thetas, sds, P, batchPr){
  K <- length(batches)
  mlist <- vector("list", K)
  for(j in seq_len(ncol(thetas))){
    mlist[[j]] <- .batchdens(x, batches, thetas[, j], sds[, j], P[, j], batchPr)
  }
  names(mlist) <- paste0("component", seq_len(K))
  mlist
}

.batchdens <- function(x, batches, thetas, sds, p1, p2){
  marginal <- matrix(NA, length(x), length(batches))
  for(b in seq_along(batches)){
    marginal[, b] <- dens(x, thetas[b], sds[b], p1[b], p2[b])
  }
  colnames(marginal) <- batches
  marginal
}

#' @export
#' @examples
#' df <- multiBatchDensities(BatchModelExample)
#' head(df)
#' @rdname ggplot-functions
multiBatchDensities <- function(model){
  dnorm_poly_multibatch(model)
}

.gg_multibatch <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  ##df <- multiBatchDensities(model)
  dat <- dnorm_poly_multibatch(model)
  nb <- nBatch(model)
  df.observed <- data.frame(y=observed(model),
                            batch=paste0("batch: ", batch(model)))
  if(missing(bins))
    bins <- nrow(df.observed)/2
  component <- y <- ..density.. <- x <- y <- NULL
  ggplot(dat, aes(x, y, group=component)) +
    geom_histogram(data=df.observed, aes(y, ..density..),
                   bins=bins,
                   inherit.aes=FALSE) +
    geom_polygon(aes(fill=component, color=component), alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_y_sqrt() +
    guides(fill=guide_legend(""), color=guide_legend("")) +
    facet_wrap(~batch, nrow=nb)
}

.gg_multibatch_copynumber <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")

  ##df <- multiBatchDensities(model)
  dat <- dnorm_poly_multibatch(model)
  dat2 <- dat
  index <- dat$component %in% seq_len(k(model))
  comp <- dat$component
  comp2 <- comp[index]
  comp2 <- as.integer(as.character(comp2))
  cn <- .remap(comp2, mapping(model))
  comp[index] <- cn
  comp <- factor(comp)
  dat$component <- comp
  nb <- nBatch(model)
  df.observed <- data.frame(y=observed(model),
                            batch=paste0("batch: ", batch(model)))
  if(missing(bins))
    bins <- nrow(df.observed)/2
  component <- ..density.. <- x <- y <- NULL
  ## how to make the density and polygon on same y-scale
  ggplot(dat, aes(x, y, group=component)) +
    geom_histogram(data=df.observed, aes(y, ..density..),
                   ##bins=bins,
                   inherit.aes=FALSE, binwidth=0.05) +
    geom_polygon(aes(fill=component, color=component),
                 alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
      scale_fill_manual(values=colors) +
      scale_y_sqrt() +
    guides(fill=guide_legend(""), color=guide_legend("")) +
    facet_wrap(~batch, nrow=nb)
}

#' @rdname ggplot-functions
setMethod("ggMultiBatch", "BatchModel", function(model, bins){
  .gg_multibatch(model, bins)
})

#' @rdname ggplot-functions
setMethod("ggMultiBatch", "MultiBatchCopyNumber", function(model, bins){
  .gg_multibatch_copynumber(model, bins)
})
