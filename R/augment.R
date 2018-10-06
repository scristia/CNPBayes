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
    summarize(nhom=sum(medians < -1),
              min_homdel=min(medians),
              nambiguous=sum(medians < -0.9))
  nzero <- sum(dat$nhom <= 2)
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
