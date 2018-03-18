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

setMethod("gatherChains", "MultiBatchModel", function(object){
  .meltMultiBatchChains(object)
})

setMethod("gatherChains", "MultiBatchPooled", function(object){
  .meltMultiBatchPooledChains(object)
})

.ggMultiBatchChains <- function(model){
  melt.ch <- gatherChains(model)
  dat.batch <- melt.ch$batch %>%
    mutate(iter=as.integer(iter))
  dat.comp <- melt.ch$comp %>%
    as.tibble
  dat.single <- melt.ch$single %>%
    as.tibble
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
  for(i in seq_along(melt.ch)){
    melt.ch[[i]]$iter <- as.integer(melt.ch[[i]]$iter)
  }
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

.gg_multibatch <- function(model, bins){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  predictive <- posteriorPredictive(model) %>%
    mutate(component=factor(component))
  predictive.summary <- predictive %>%
    group_by(model, batch) %>%
    summarize(n=n())
  predictive <- left_join(predictive, predictive.summary,
                          by=c("model", "batch")) %>%
    mutate(batch=paste("Batch", batch))
  colors <- colors[seq_len(k(model))]
  ##df <- multiBatchDensities(model)
  full.data <- tibble(y=y(model),
                      batch=batch(model)) %>%
    mutate(batch=paste("Batch", batch)) %>%
    mutate(model=modelName(model))
  ##xlimit <- c(-5, 1)
  batch.data <- full.data %>%
    group_by(batch, model) %>%
    summarize(n=n()) %>%
    mutate(x=-Inf, y=Inf)
  fig <- ggplot(predictive, aes(x=y, n_facet=n,
                                y=..count../n_facet,
                                fill=component)) +
    geom_histogram(data=full.data, aes(y, ..density..),
                   bins=bins,
                   inherit.aes=FALSE,
                   color="gray70",
                   fill="gray70",
                   alpha=0.1) +
    geom_density(adjust=1, alpha=0.4, size=0.75, color="gray30") +
    geom_text(data=batch.data, aes(x=x, y=y, label=paste0("  n=", n)),
              hjust="inward", vjust="inward",
              inherit.aes=FALSE,
              size=3) +
    ## show marginal density
    theme(panel.background=element_rect(fill="white"),
          axis.line=element_line(color="black"),
          legend.position="bottom",
          legend.direction="horizontal") +
    ##facet_grid(batch~model,
    facet_wrap(~batch, ncol=1, strip.position="right") +
    ##labeller=labeller(model=ml)) +
    scale_y_sqrt() +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    xlab("average copy number") +
    ylab("density") +
    ##coord_cartesian(xlim=xlimit) +
    guides(color=FALSE) ##+
  return(fig)
##    ggtitle(cnp.id)
  ##df <- multiBatchDensities(model)
##  dat <- dnorm_poly_multibatch(model)
##  nb <- nBatch(model)
##  df.observed <- tibble(y=observed(model),
##                        batch=paste0("batch: ", batch(model)))
##  if(missing(bins))
##    bins <- nrow(df.observed)/2
##  component <- y <- ..density.. <- x <- y <- NULL
##  dat <- as.tibble(dat)
##  ##
##  ## Replace histogram with geom_polygon.  Need x and y
##  ##  -
##  ggplot(dat, aes(x, y, group=component)) +
##    ##geom_density() +
##    ##facet_wrap(~batch, nrow=nb) +
##    geom_histogram(data=df.observed, aes(y, ..density..),
##                   bins=bins,
##                   inherit.aes=FALSE) +
##    geom_polygon(aes(fill=component, color=component), alpha=0.4) +
##    xlab("quantiles") + ylab("density") +
##    scale_color_manual(values=colors) +
##    scale_fill_manual(values=colors) +
##    scale_y_sqrt() +
##    guides(fill=guide_legend(""), color=guide_legend("")) +
##    facet_wrap(~batch, nrow=nb)
}

mb_pred_data <- function(model, predict){
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
  dat2 <- dat2[, c("y", "predictive", "batch")]
  dat2
}

mb_predictive <- function(model, predict, adjust=1/3){
  dat <- mb_pred_data(model, predict)
  predictive <- NULL
  fig <- ggplot(dat, aes(y, fill=predictive)) +
    geom_density(alpha=0.4, adjust=adjust) +
    facet_wrap(~batch, ncol=1) +
    guides(fill=guide_legend(title="")) +
    theme(panel.background=element_rect(fill="white"))
  fig
}

.predictiveDataTable <- function(model, predict){
  model.name <- modelName(model)
  if(class(model) %in% c("SingleBatchModel", "SingleBatchPooled")){
    dat <- sb_pred_data(model, predict)
    dat$model <- model.name
  } else {
    dat <- mb_pred_data(model, predict)
    dat$model <- model.name
  }
  dat
}

#' Create a tibble of posterior predictive distributions for a list of models
#'
#' @param models list of models
#' @examples
#'   sb <- SingleBatchModelExample
#'   mcmcParams(sb) <- McmcParams(iter=500, burnin=50)
#'   sb <- posteriorSimulation(sb)
#'   models <- list(MultiBatchModelExample, sb)
#'   tab <- predictiveDataTable(models)
#'   \dontrun{
#'     library(ggplot2)
#'     ggplot(tab, aes(y, fill=predictive)) +
#'       geom_density(alpha=0.4, adjust=1/2) +
#'       facet_grid(model~batch) +
#'       guides(fill=guide_legend(title="")) +
#'       theme(panel.background=element_rect(fill="white"))
#'   }
#' @export
predictiveDataTable <- function(models){
  tab.list <- vector("list", length(models))
  for(i in seq_along(models)){
    pred <- posteriorPredictive(models[[i]])
    tab.list[[i]] <- .predictiveDataTable(models[[i]], pred)
  }
  tab <- do.call(rbind, tab.list)
  ## assume models were ordered by marginal likelihood
  mod.names <- sapply(models, modelName)
  tab$model <- factor(tab$model, levels=mod.names)
  tab
}

sb_pred_data <- function(model, predict){
  dat <- tibble(y=y(model), predictive="empirical")
  ##browser()
  ##predict.bak=predict
  colnames(predict)[2] <- "predictive"
  predict$predictive <- "posterior\npredictive"
  predict$batch <- factor("batch 1")
  dat$batch <- factor("batch 1")
  predictive <- NULL
  dat2 <- rbind(dat, predict) %>%
    mutate(predictive=factor(predictive,
                             levels=c("empirical", "posterior\npredictive")))
  dat2
}

sb_predictive <- function(model, predict, adjust=1/3){
  dat2 <- sb_pred_data(model, predict)
  predictive <- NULL
  fig <- ggplot(dat2, aes(y, fill=predictive, color=predictive)) +
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

multibatch_figure <- function(theoretical, empirical, model){
  nb <- nBatch(model)
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  scale_col <-  scale_color_manual(values=colors)
  scale_fill <- scale_fill_manual(values=colors)
  scale_y <- scale_y_sqrt()
  lrr <- NULL
  ..count.. <- NULL
  ghist <- geom_histogram(data=empirical, aes(lrr, ..count..), binwidth=0.01, inherit.aes=FALSE)
  gobj <- ggplot() + ghist + facet_wrap(~batch)
  gb <- ggplot_build(gobj)
  ylimit <- gb$layout$panel_ranges[[1]][["y.range"]]
  theoretical.sum <- group_by(theoretical, batch, component) %>%
    summarize(maxy=max(y))
  theoretical$y <- rescale(theoretical$y, c(0, ylimit[2]))
  component <- x <- y <- NULL
  gpolygon <-  geom_polygon(aes(x, y, fill=component, color=component), alpha=0.4)
  ggplot(theoretical) +
    ghist +
    gpolygon +
    scale_col +
    scale_fill +
    xlab("quantiles") + ylab("count") +
    guides(fill=guide_legend(""), color=guide_legend("")) +
    facet_wrap(~batch, nrow=nb, as.table=TRUE, scales="free_y")
}

.gg_multibatch_pooled <- function(model, bins){
  fig <- .gg_multibatch(model, bins)
  fig
}



.gg_multibatch_copynumber <- function(model, bins=400){
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  predictive <- posteriorPredictive(model) %>%
    mutate(component=factor(component))
  predictive.summary <- predictive %>%
    group_by(model, batch) %>%
    summarize(n=n())
  predictive <- left_join(predictive, predictive.summary,
                          by=c("model", "batch")) %>%
    mutate(batch=paste("Batch", batch))
  zz <- map_z(model)
  comp_labels <- mapping(model)
  predictive$copynumber <- comp_labels[predictive$component]
  colors <- colors[seq_along(comp_labels)]
  ##df <- multiBatchDensities(model)
  full.data <- tibble(y=y(model),
                      batch=batch(model)) %>%
    mutate(batch=paste("Batch", batch)) %>%
    mutate(model=modelName(model))
  ##xlimit <- c(-5, 1)
  batch.data <- full.data %>%
    group_by(batch, model) %>%
    summarize(n=n()) %>%
    mutate(x=-Inf, y=Inf)
  fig <- ggplot(predictive, aes(x=y, n_facet=n,
                                y=..count../n_facet,
                                fill=copynumber)) +
    geom_histogram(data=full.data, aes(y, ..density..),
                   bins=bins,
                   inherit.aes=FALSE,
                   color="gray70",
                   fill="gray70",
                   alpha=0.1) +
    geom_density(adjust=1, alpha=0.4, size=0.75, color="gray30") +
    geom_text(data=batch.data, aes(x=x, y=y, label=paste0("  n=", n)),
              hjust="inward", vjust="inward",
              inherit.aes=FALSE,
              size=3) +
    ## show marginal density
    theme(panel.background=element_rect(fill="white"),
          axis.line=element_line(color="black"),
          legend.position="bottom",
          legend.direction="horizontal") +
    ##facet_grid(batch~model,
    facet_wrap(~batch, ncol=1, strip.position="right") +
    ##labeller=labeller(model=ml)) +
    scale_y_sqrt() +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    xlab("average copy number") +
    ylab("density") +
    ##coord_cartesian(xlim=xlimit) +
    guides(color=FALSE) ##+
  return(fig)
##    ggtitle(cnp.id)
##
##  empirical <- tibble(lrr=observed(model),
##                      batch=paste0("batch: ", batch(model)),
##                      z=map_z(model))
##  empirical$cn <- .remap(empirical$z, mapping(model))
##  empirical2 <- empirical %>%
##    mutate(batch="overall")
##  empirical3 <- bind_rows(empirical, empirical2) %>%
##    mutate(batch=factor(batch, levels=c(unique(empirical$batch),
##                                        "overall")))
##
##  scale_col <-  scale_color_manual(values=colors)
##  scale_fill <- scale_fill_manual(values=colors)
##  scale_y <- scale_y_sqrt()
##  ..count.. <- lrr <- NULL
##  ghist <- geom_histogram(data=empirical3, aes(lrr, ..count..), binwidth=0.01, inherit.aes=FALSE)
##  gobj <- ggplot() + ghist + facet_wrap(~batch)
##  gb <- ggplot_build(gobj)
##  ylimit <- gb$layout$panel_ranges[[1]][["y.range"]]
##
##  theoretical$y <- scales::rescale(theoretical$y, c(0, ylimit[2]))
##  x <- y <- cn <- NULL
##  gpolygon <-  geom_polygon(aes(x, y, fill=cn, color=cn), alpha=0.4)
##  nb <- nBatch(model)
##  ggplot(theoretical) +
##    ghist +
##    gpolygon +
##    scale_col +
##    scale_fill +
##    xlab("quantiles") + ylab("count") +
##    guides(fill=guide_legend(""), color=guide_legend("")) +
##    facet_wrap(~batch, nrow=nb, as.table=TRUE, scales="free_y")
}

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchCopyNumber-method
setMethod("ggMixture", "MultiBatchCopyNumber", function(model, bins=100){
  .gg_multibatch_copynumber(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchCopyNumberPooled-method
setMethod("ggMixture", "MultiBatchCopyNumberPooled", function(model, bins=100){
  .gg_multibatch_copynumber(model, bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchModel-method
setMethod("ggMixture", "MultiBatchModel", function(model, bins=100){
  .gg_multibatch(model, bins=bins)
})

#' @export
#' @rdname ggplot-functions
#' @aliases ggMixture,MultiBatchPooled-method
setMethod("ggMixture", "MultiBatchPooled", function(model, bins=100){
  .gg_multibatch_pooled(model, bins)
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
