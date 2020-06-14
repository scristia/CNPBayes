##
## problem: batch not clustered with other batches and does not have hom del comp, while other batches clearly do
##
## Solution 1: augment data for batches without homdel by simulating 5 observations
## - simulate(rnorm(median(min), 0.3))
## - fit model
## - provide posterior probs only for non-augmented data
## potential problems : cutoff of 1 is  very arbitrary
##   - would need to filter augmented data in ggMixture; predictive distribution may look funny
##
##
augmentData <- function(full.data){
  full.data$augmented <- FALSE
  dat <- group_by(full.data, batch) %>%
    summarize(nhom=sum(medians < -1, na.rm=TRUE),
              min_homdel=min(medians, na.rm=TRUE),
              nambiguous=sum(medians < -0.9, na.rm=TRUE))
  nzero <- sum(dat$nhom == 0, na.rm=TRUE)
  if(nzero == 0 || nzero == nrow(dat))
    return(full.data)
  dat2 <- dat %>%
    filter(nhom == 0)
  if(!"id" %in% colnames(full.data)){
    full.data$id <- seq_len(nrow(full.data)) %>%
      as.character
  }
  current <- full.data
  for(i in seq_len(nrow(dat2))){
    augment <- filter(full.data, batch == dat2$batch[i])  %>%
      "["(1:5, ) %>%
      mutate(medians = rnorm(5, mean(dat$min_homdel), 0.3),
             augmented=TRUE,
             id=paste0(id, "*"))
    latest <- bind_rows(current, augment)
    current <- latest %>%
      arrange(batch_index)
  }
  current
}

useSingleBatchValues <- function(mb, sb){
  cv.mb <- current_values(mb)
  cv.sb <- current_values(sb)

  cv.mb$theta <- replicate(numBatch(mb), cv.sb$theta) %>%
    "["(1,,) %>%
    t
  cv.mb$sigma2 <- replicate(numBatch(mb), cv.sb$sigma2) %>%
    "["(1,,) %>%
    t
  cv.mb$p <- replicate(numBatch(mb), cv.sb$p) %>%
    "["(1,,) %>%
    t
  current_values(mb) <- cv.mb
  modes(mb) <- cv.mb
  mb
}


setMethod("augmentData2", "MultiBatchList", function(object){
  ##
  ## - only makes sense to do this if multibatch models with 3 or 4 components are included in the list
  ##
  sp <- specs(object) %>%
    filter(k == 3 & substr(model, 1, 2) == "MB")
  if(nrow(sp) == 0) return(object)
  sp <- sp[1, ]
  object2 <- object[ specs(object)$model %in% sp$model ]
  ix <- order(specs(object2)$k, decreasing=FALSE)
  ## order models by number of components.
  object2 <- object2[ix]
  SB <- toSingleBatch(object2[[1]])
  iter(SB) <- max(iter(object), 150L)
  SB <- posteriorSimulation(SB)
  ##
  ## Running MCMC for SingleBatch model to find posterior predictive distribution
  ##
  message("Checking whether possible homozygous deletions occur in only a subset of batches...")
  modes(SB) <- computeModes(SB)
  SB <- setModes(SB)
  mn.sd <- c(theta(SB)[1], sigma(SB)[1])
  limits <- mn.sd[1] + c(-1, 1)*2*mn.sd[2]
  ## record number of observations in each batch that are within 2sds of the mean
  freq.del <- assays(object) %>%
    group_by(batch) %>%
    summarize(n = sum(oned < limits[[2]]))
  fewobs <- freq.del$n <= 2
  iszero <- freq.del$n == 0
  if( all(iszero) ){
    ## no homozygous deletions
    assays(object)$is_simulated <- FALSE
    return(object)
  }
  if( !any(fewobs) ){
    ## many homozygous deletions in each batch
    assays(object)$is_simulated <- FALSE
    return(object)
  }
  ## Some of the batches have 2 or fewer homozygous deletions
  ## Augment data with 10 observations to allow fitting this component
  ##
  ## - sample a minimum of 10 observations (with replacement) from the posterior predictive distribution of the other batches
  ##
  batches <- freq.del$batch [ fewobs ]
  ##zerobatch <- freq.del$batch[ freq.del$n == 0 ]
  dat <- assays(object2[[1]])
  expected_homdel <- modes(SB)[["p"]][1] * table(dat$batch)
  expected_homdel <- ceiling(expected_homdel [ unique(dat$batch) %in% batches ])
  nsample <- pmax(10L, expected_homdel)
  pred <- predictiveTibble(SB) %>%
    ##filter(!(batch %in% batches)  & component == 0) %>%
    filter(component == 0) %>%
    "["(sample(seq_len(nrow(.)), sum(nsample), replace=TRUE), ) %>%
    select(oned) %>%
    mutate(batch=rep(batches, nsample),
           id=paste0("augment_", seq_len(nrow(.))),
           oned=oned+rnorm(nrow(.), 0, mn.sd[2])) %>%
    select(c(id, oned, batch))
  newdat <- bind_rows(assays(object), pred) %>%
    arrange(batch) %>%
    mutate(is_simulated = seq_len(nrow(.)) %in% grep("augment_", id))
  mbl <- MultiBatchList(data=newdat,
                        parameters=parameters(object))
  mbl
})


augmentTest <- function(object){
  ##
  ## - only makes sense to do this if multibatch models
  ##   with 3 or 4 components are included in the list
  ##
  ##
  ## Idea:
  ##  1. Run SingleBatch independently for each batch
  ##  2. Assess if there is a potential problem
  ##     - components with high standard deviation, or models with small standard deviation of thetas
  ##     - components with very few observations
  ##  3. If no problems detected, return object
  ##
  ##   Unusually high standard deviations
  ##     - is homozygous deletion component missing?
  ##       assign well estimated components to theta1 theta2 of theta matrix
  ##       set theta0 to NA for these batches
  ##
  ##   Small standard deviations
  ##     - set components with most observations to theta 1 theta2 of theta matrix
  ##     - set theta0 to NA
  ##  4. Impute missing theta 10 times assuming MVN
  ##  5. Augment data with the imputed thetas
  ##
  mb <- object[["MB3"]]
  iter(mb) <- max(iter(object), 150L)
  burnin(mb) <- 0L
  mbm <- as(mb, "MultiBatchModel")
  zfreq <- tableBatchZ(mbm)
  if(any(zfreq == 0)){

  }
  mb <- posteriorSimulation(mb)
  ub <- unique(batch(mb))
  mb.list <- vector("list", length(ub))
  for(i in seq_along(ub)){
    B <- ub[i]
    mb2 <- mb[ batch(mb) == B ]
    mb.list[[i]] <- posteriorSimulation(mb2)
  }
  th <- lapply(mb.list, theta) %>%
    do.call(rbind, .) %>%
    round(2)
  sds <- lapply(mb.list, sigma) %>%
    do.call(rbind, .) %>%
    round(2)
  zz <- lapply(mb.list, function(x) table(z(x))) %>%
    do.call(rbind, .)

  r1 <- order(rowSds(th), decreasing=TRUE)

  batchfreq <- table(batch(object))
  B <- which(batchfreq==max(batchfreq))[[1]]
  mb <- object[["MB3"]]
  sb <- mb[ batch(mb) == B ]
  iter(mb) <- iter(sb) <- max(iter(object), 150L)
  sb <- posteriorSimulation(sb)
  modes(sb) <- computeModes(sb)
  sb <- setModes(sb)
  MB <- useSingleBatchValues(mb=mb, sb=sb)
  mbm <- as(MB, "MultiBatchModel")
  z(mbm) <- update_z(mbm)
  zfreq <- tableBatchZ(mbm)
  if(any(zfreq)==0){
    
  }
  MB <- posteriorSimulation(MB)
  modes(MB) <- computeModes(MB)
  MB <- setModes(MB)
  sds <- sigma(MB)
  th <- theta(MB)
  fc <- sds[, 1]/min(sds[, 1])
  if(!any(fc > 2.5 )) {
    assays(object)$is_simulated <- FALSE
    return(object)
  }
  if(any(fc > 2.5)){
    th1 <- th[, 1]
    th1[ fc > 2.5 ] <- NA
    th[, 1] <- NA

    th.imputed <- replicate(10, Impute, simplify=FALSE)
    th.imputed2 <- lapply(th.imputed, function(x) x$yimp[, 1]) %>%
      do.call(cbind, .)
    th.imputed <- th.imputed2[which(is.na(th1)), drop=FALSE]
  }
  dat <- assays(object[[1]])
  newdat <- bind_rows(assays(object), pred) %>%
    arrange(batch) %>%
    mutate(is_simulated = seq_len(nrow(.)) %in% grep("augment_", id))
  mbl <- MultiBatchList(data=newdat,
                        parameters=parameters(object))
  mbl
}

impute <- function(model, loc.scale, start.index){
  x <- assays(model) %>%
    group_by(batch) %>%
    summarize(N=n()) %>%
    ##left_join(tab, by="batch") %>%
    left_join(loc.scale, by="batch") %>%
    mutate(##n = ifelse(is.na(n), 0, n),
      ##phat=n/N,
      expected=2*ceiling(N*max(phat, na.rm=TRUE))) %>%
    filter(!is.na(theta)) %>%
    mutate(expected=pmax(expected, 5))
    ##filter(n < 3)
  imp.list <- vector("list", nrow(x))
  for(i in seq_along(imp.list)){
    imp.list[[i]] <- rnorm(x$expected[i], mean=x$theta[i],
                           sd=sqrt(x$sigma2[i]))
  }
  imp <- unlist(imp.list)
  index <- seq_along(imp) + start.index - 1
  impdat <- tibble(id=paste0("augment_", index),
                   oned=imp,
                   provisional_batch=NA,
                   ##likely_hd=TRUE,
                   likely_deletion=TRUE,
                   batch=rep(x$batch, x$expected),
                   is_simulated=TRUE)
  batch_mapping <- assays(model) %>%
    group_by(batch) %>%
    summarize(batch_labels=unique(batch_labels))
  impdat2 <- left_join(impdat, batch_mapping, by="batch")
  impdat2
}
