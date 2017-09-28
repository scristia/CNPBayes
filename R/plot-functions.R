.meltMultiBatchChains <- function(model){
  ch <- chains(model)
  th <- as.data.frame(theta(ch))
  th$param <- "theta"
  th$iter <- factor(1:nrow(th))
  ##
  ## to suppress annoying notes
  ##
  param <- iter <- NULL
  th.m <- gather(th, key="variable", values=-c(param, iter)) %>%
    as.tibble
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
  . <- value <- NULL
  s.m <- gather(s, key="variable", values=-c(param, iter)) %>%
    as.tibble
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

  dat.batch <- rbind(th.m, s.m)
  dat.comp <- rbind(mus.m, taus.m, prob.m)
  dat <- rbind(nu0.m, s20.m)
  dat$iter <- as.integer(dat$iter)
  list(batch=dat.batch, comp=dat.comp, single=dat)
}

.meltMultiBatchPooledChains <- function(model){
  ch <- chains(model)
  th <- as.data.frame(theta(ch))
  th$param <- "theta"
  th$iter <- factor(1:nrow(th))
  . <- param <- iter <- NULL
  th.m <- gather(th, key="variable", values=-c(param, iter)) %>%
    as.tibble
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

  value <- NULL
  s <- as.data.frame(sigma(ch))
  s.m <- s %>% as.tibble %>%
    mutate(iter=seq_len(iter(model))) %>%
    gather(param, value, -iter) %>%
    mutate(param="sigma",
           iter=factor(iter),
           batch=factor(rep(uniqueBatch(model), each=iter(model))))

  nu0.m <- as.tibble(nu.0(ch)) %>%
    mutate(iter=th$iter,
           param="nu0")
  s20.m <- as.tibble(sigma2.0(ch)) %>%
    mutate(iter=th$iter,
           param="s20")

  mus.m <- as.tibble(mu(ch)) %>% 
    mutate(iter=th$iter) %>%
    gather(param, value, -iter) %>%
    mutate(param="mu",
           comp=rep(seq_len(k(model)), each=iter(model))) %>%
    mutate(comp=factor(paste0("k=", .$comp)))

  prob.m <- as.tibble(p(ch)) %>%
    mutate(iter=th$iter) %>%
    gather(param, value, -iter) %>%
    mutate(param="p",
           comp=mus.m$comp)

  taus.m <- as.tibble(tau(ch)) %>%
    mutate(iter=th$iter) %>%
    gather(param, value, -iter) %>%
    mutate(param="tau",
           comp=mus.m$comp)

  ##dat.batch <- bind_rows(th.m, s.m)
  dat.comp <- rbind(mus.m, taus.m, prob.m)
  dat <- bind_rows(nu0.m, s20.m)
  dat$iter <- as.integer(dat$iter)
  list(theta=th.m,
       sigma=s.m,
       comp=dat.comp,
       single=dat)
}


meltSingleBatchChains <- function(model){
  parameter <- NULL
  ch <- chains(model)
  th <- as.tibble(theta(ch)) %>%
    set_colnames(paste0("theta", seq_len(k(model))))
  th$iter <- 1:nrow(th)
  th.m <- gather(th, key="parameter", value="value", -iter) %>%
    mutate(comp=gsub("theta", "", parameter))
  sig <- as.tibble(sigma(ch)) %>%
    set_colnames(paste0("sigma", seq_len(k(model))))
  if(class(model) == "SingleBatchPooled"){
    sig$iter <- 1:nrow(sig)
  } else sig$iter <- th$iter
  s.m <- gather(sig, key="parameter", value="value", -iter) %>%
    mutate(comp=gsub("sigma", "", parameter))
  nu0 <- as.tibble(nu.0(ch)) %>%
    set_colnames("nu0")
  nu0$iter <- th$iter
  nu0.m <- gather(nu0, key="parameter", value="value", -iter)

  s20 <- as.tibble(sigma2.0(ch)) %>%
    set_colnames("s20")
  s20$iter <- th$iter
  s20.m <- gather(s20, key="parameter", value="value", -iter)

  mus <- as.tibble(mu(ch))%>%
    set_colnames("mu")
  mus$iter <- th$iter
  mus.m <- gather(mus, key="parameter", value="value", -iter)

  prob <- as.tibble(p(ch)) %>%
    set_param_names("p")
  prob$iter <- th$iter
  prob.m <- gather(prob, key="parameter", value="value", -iter) %>%
    mutate(comp=gsub("p", "", parameter))

  taus <- as.tibble(tau(ch)) %>%
    set_colnames("tau")
  taus$iter <- th$iter
  taus.m <- gather(taus, key="parameter", value="value", -iter)

  ll <- as.tibble(log_lik(ch)) %>%
    set_colnames("log lik") %>%
    mutate(iter=th$iter) %>%
    gather(key="parameter", value="value", -iter)
  dat.comp <- bind_rows(th.m,
                        s.m,
                        prob.m)

  dat <- bind_rows(nu0.m,
                   s20.m,
                   mu=mus.m,
                   tau=taus.m,
                   loglik=ll)
  list(comp=dat.comp,
       single=dat)
}

meltSingleBatchPooledChains <- function(model){
  ch <- chains(model)
  th <- as.tibble(theta(ch)) %>%
    set_colnames(paste0("theta", seq_len(k(model))))
  th$iter <- 1:nrow(th)
  parameter <- NULL
  th.m <- gather(th, key="parameter", value="value", -iter) %>%
    mutate(comp=gsub("theta", "", parameter))
  sig <- as.tibble(sigma(ch)) %>%
    set_colnames("value") %>%
    mutate(iter=seq_len(iter(model))) %>%
    mutate(parameter="sigma")
##  if(class(model) == "SingleBatchPooled"){
##    sig$iter <- 1:nrow(sig)
##  } else sig$iter <- th$iter
##  s.m <- gather(sig, key="parameter", value="value", -iter) %>%
##    mutate(comp=gsub("sigma", "", parameter))
  nu0 <- as.tibble(nu.0(ch)) %>%
    set_colnames("nu0")
  nu0$iter <- th$iter
  nu0.m <- gather(nu0, key="parameter", value="value", -iter)
  nu0$iter <- th$iter
  nu0.m <- gather(nu0, key="parameter", value="value", -iter)

  prob <- as.tibble(p(ch)) %>%
    set_param_names("p")
  prob$iter <- th$iter
  prob.m <- gather(prob, key="parameter", value="value", -iter) %>%
    mutate(comp=gsub("p", "", parameter))

  s20 <- as.tibble(sigma2.0(ch)) %>%
    set_colnames("s20")
  s20$iter <- th$iter
  s20.m <- gather(s20, key="parameter", value="value", -iter)

  mus <- as.tibble(mu(ch))%>%
    set_colnames("mu")
  mus$iter <- th$iter
  mus.m <- gather(mus, key="parameter", value="value", -iter)

  taus <- as.tibble(tau(ch)) %>%
    set_colnames("tau")
  taus$iter <- th$iter
  taus.m <- gather(taus, key="parameter", value="value", -iter)

  ll <- as.tibble(log_lik(ch)) %>%
    set_colnames("log lik") %>%
    mutate(iter=th$iter) %>%
    gather(key="parameter", value="value", -iter)
  dat.comp <- bind_rows(th.m,
                        ##sig,
                        prob.m)
  sig <- sig[ , colnames(nu0.m) ]
  dat <- bind_rows(sig,
                   nu0.m,
                   s20.m,
                   mu=mus.m,
                   tau=taus.m,
                   loglik=ll)
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
#' sb <- SingleBatchModelExample
#' plist.sb <- ggChains(sb)
#' \dontrun{
#' ## chains for parameter vectors of length k
#' plist.sb[["comp"]]
#' ## chains for parameters vectors of length 1
#' plist.sb[["single"]]
#' }
#'
#' mb <- MultiBatchModelExample
#' plist.mb <- ggChains(mb)
#' \dontrun{
#' ## chains for parameters that are batch- and component-specific
#' plist.mb[["batch"]]
#' ## chains for parameters vectors of length k
#' plist.mb[["comp"]]
#' ## chains for parameter vectors of length 1
#' plist.mb[["single"]]
#' }
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
    facet_wrap(~parameter, scales="free_y")
  p.single <- ggplot(dat.single, aes(iter, value, group="parameter")) +
    geom_point(size=0.3, color="gray") +
    geom_line(color="gray") +
    facet_wrap(~parameter, scales="free_y")
  list(comp=p.comp, single=p.single)
}


ggSingleBatchPooledChains <- function(model){
  melt.ch <- gatherChains(model)
  dat.comp <- melt.ch[["comp"]]     ## component-specific
  dat.single <- melt.ch[["single"]] ## single-parameter
  iter <- value <- comp <- param <- NULL
  p.comp <- ggplot(dat.comp, aes(iter, value, group=comp)) +
    geom_point(size=0.3, aes(color=comp)) +
    geom_line(aes(color=comp)) +
    facet_wrap(~parameter, scales="free_y")
  p.single <- ggplot(dat.single, aes(iter, value, group="parameter")) +
    geom_point(size=0.3, color="gray") +
    geom_line(color="gray") +
    facet_wrap(~parameter, scales="free_y")
   list(comp=p.comp, single=p.single)
}

setMethod("gatherChains", "MultiBatchModel", function(object){
  .meltMultiBatchChains(object)
})

setMethod("gatherChains", "SingleBatchPooled", function(object){
  meltSingleBatchPooledChains(object)
})

setMethod("gatherChains", "MultiBatchPooled", function(object){
  .meltMultiBatchPooledChains(object)
})

.ggMultiBatchChains <- function(model){
  melt.ch <- gatherChains(model)
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

.ggMultiBatchPooledChains <- function(model){
  melt.ch <- gatherChains(model)
  iter <- value <- batch <- param <- comp <- NULL
  p.theta <- ggplot(melt.ch$theta, aes(iter, value, group=batch)) +
    geom_point(size=0.3, aes(color=batch)) +
    geom_line(aes(color=batch)) +
    facet_grid(param ~ comp, scales="free_y")

  p.sigma <- ggplot(melt.ch$sigma, aes(iter, value, group=batch)) +
    geom_point(size=0.3, aes(color=batch)) +
    geom_line(aes(color=batch)) 
    ##facet_grid(param ~ comp, scales="free_y")

  p.comp <- ggplot(melt.ch$comp, aes(iter, value, group=comp)) +
    geom_point(size=0.3, aes(color=comp)) +
    geom_line(aes(color=comp)) +
    facet_wrap(~param, scales="free_y")

  p.single <- ggplot(melt.ch$single, aes(iter, value, group="param")) +
    geom_point(size=0.3, color="gray") +
    geom_line(color="gray") +
    facet_wrap(~param, scales="free_y")
  list(theta=p.theta, sigma=p.sigma, comp=p.comp, single=p.single)
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
  df.overall$component <- "SB"
  df <- rbind(df, df.overall)
  df$component <- factor(df$component, levels=c("SB", seq_along(mean)))
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
  ##if(class(model) == "SingleBatchPooled"){
  if(length(sds) != k(model)){
    sds <- rep(sds, k(model))
  }
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
  comp <- factor(comp, levels=c("SB", levs))
  comp
}

.gg_singlebatch_copynumber <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  df.observed <- tibble(y=observed(model))
  ## see stat_function
  ##  ggplot(df, aes(x, d, group=name)) +
  if(missing(bins))
    bins <- nrow(df.observed)/2
  dat <- dnorm_poly(model) %>%
    as.tibble
  dat$component <- .relabel_component(dat, model) 
  component <- x <- y <- ..density.. <- NULL
  ggplot(dat, aes(x, y, group=component)) +
    geom_histogram(data=df.observed, aes(y, ..density..),
                   bins=bins,
                   inherit.aes=FALSE) +
    geom_polygon(aes(fill=component, color=component), alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_y_sqrt() +
    guides(fill=guide_legend(""), color=guide_legend(""))
}


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

dnorm_poly_multibatch_pooled <- function(model){
  mixprob <- p(model)
  means <- theta(model)
  nc <- k(model)
  sds <- matrix(rep(sigma(model), nc),
                nrow(means), nc)
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

.gg_multibatch <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  ##df <- multiBatchDensities(model)
  dat <- dnorm_poly_multibatch(model)
  nb <- nBatch(model)
  df.observed <- tibble(y=observed(model),
                        batch=paste0("batch: ", batch(model)))
  if(missing(bins))
    bins <- nrow(df.observed)/2
  component <- y <- ..density.. <- x <- y <- NULL
  dat <- as.tibble(dat)
  ##
  ## Replace histogram with geom_polygon.  Need x and y
  ##  -
  ggplot(dat, aes(x, y, group=component)) +
    ##geom_density() +
    ##facet_wrap(~batch, nrow=nb) +
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

mb_predictive <- function(model, predict, adjust=1/3){
  batches <- factor(batch(model), labels=paste("batch", unique(batch(model))))
  dat <- tibble(y=y(model),
                batch=batches,
                predictive="empirical")
  predict$batch <- factor(predict$batch, labels=paste("batch", unique(batch(model))))
  colnames(predict)[3] <- "predictive"
  predict$predictive <- "posterior\npredictive"
  predictive <- NULL
  dat2 <- rbind(dat, predict) %>%
    mutate(predictive=factor(predictive,
                             levels=c("empirical", "posterior\npredictive")))
  fig <- ggplot(dat2, aes(y, fill=predictive)) +
    geom_density(alpha=0.4, adjust=adjust) +
    facet_wrap(~batch, ncol=1) +
    guides(fill=guide_legend(title="")) +
    theme(panel.background=element_rect(fill="white"))
  fig
}

sb_predictive <- function(model, predict, adjust=1/3){
  dat <- tibble(y=y(model), predictive="empirical")
  colnames(predict)[2] <- "predictive"
  predict$predictive <- "posterior\npredictive"
  predictive <- NULL
  dat2 <- rbind(dat, predict) %>%
    mutate(predictive=factor(predictive,
                             levels=c("empirical", "posterior\npredictive")))
  fig <- ggplot(dat2, aes(y, fill=predictive)) +
    geom_density(alpha=0.4, adjust=adjust) +
    guides(fill=guide_legend(title="")) +
    theme(panel.background=element_rect(fill="white"))
  fig
}

#' Compare the posterior predictive distribution to the empirical data
#'
#' @param model a SB, MB, SBP, or MBP model
#' @param predict a \code{tibble} of the posterior predictive values, batch (only for MB and MBP models), and mixture component assignments
#' @param adjust a length-one numeric vector passed to \code{geom_density} -- controls the smoothness of the kernal density
#' @examples
#'   bmodel <- MultiBatchModelExample
#'   mp <- McmcParams(iter=500, burnin=150, nStarts=4)
#'   mcmcParams(bmodel) <- mp
#'   \dontrun{
#'      ## this is preferred to posteriorSimulation, but takes longer
#'      bmodel <- gibbs(model="MB", dat=y(bmodel), mp=mp, hp.list=hpList()[["MB"]],
#'                      batches=batch(bmodel))
#'   }
#'   bmodel <- posteriorSimulation(bmodel)
#'   tab <- posteriorPredictive(bmodel)
#'   ggPredictive(bmodel, tab)
#' @export
#' @return a `gg` object
ggPredictive <- function(model, predict, adjust=1/3){
  if(!isSB(model)){
    fig <- mb_predictive(model, predict, adjust)
  } else {
    fig <- sb_predictive(model, predict, adjust)
  }
  fig
}

.gg_multibatch_pooled <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  ##df <- multiBatchDensities(model)
  dat <- dnorm_poly_multibatch_pooled(model)
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
  if(class(model) == "MultiBatchCopyNumber"){
    dat <- dnorm_poly_multibatch(model)
  } else {
    dat <- dnorm_poly_multibatch_pooled(model)
  }
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

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchCopyNumber-method
setMethod("ggMixture", "MultiBatchCopyNumber", function(model, bins){
  .gg_multibatch_copynumber(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchCopyNumberPooled-method
setMethod("ggMixture", "MultiBatchCopyNumberPooled", function(model, bins){
  .gg_multibatch_copynumber(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,SingleBatchModel-method
setMethod("ggMixture", "SingleBatchModel", function(model, bins){
  .gg_singlebatch(model, bins)
})


#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchModel-method
setMethod("ggMixture", "MultiBatchModel", function(model, bins){
  .gg_multibatch(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchPooled-method
setMethod("ggMixture", "MultiBatchPooled", function(model, bins){
  .gg_multibatch_pooled(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,SingleBatchCopyNumber-method
setMethod("ggMixture", "SingleBatchCopyNumber", function(model, bins){
  .gg_singlebatch_copynumber(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggChains,MultiBatchModel-method
setMethod("ggChains", "MultiBatchModel", function(model){
  .ggMultiBatchChains(model)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggChains,MultiBatchPooled-method
setMethod("ggChains", "MultiBatchPooled", function(model){
  .ggMultiBatchPooledChains(model)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggChains,SingleBatchPooled-method
setMethod("ggChains", "SingleBatchPooled", function(model){
  ggSingleBatchPooledChains(model)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggChains,SingleBatchModel-method
setMethod("ggChains", "SingleBatchModel", function(model){
  ggSingleBatchChains(model)
})
