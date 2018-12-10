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
  ##
  ## Running MCMC for SingleBatch model to find posterior predictive distribution
  ##
  message("Checking whether possible homozygous deletions occur in only a subset of batches...")
  iter(SB) <- max(iter(object), 150L)
  SB <- posteriorSimulation(SB)
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
