trioData <- function(id=character(),
                     oned=numeric(),
                     batch=integer(),
                     member=character()){
  tibble(id=id,
         oned=oned,
         batch=batch,
         member=member)
}

trio_spec <- function(model, data) {
  models <- c("TSB", "TSBP")
  K <- 1:5
  avail.models <- lapply(models, paste0, K) %>%
    unlist
  if(missing(model)) model <- "TSB3" else model <- model[1]
  if(!model %in% avail.models) stop("Model not recognized")
  number_batches <- length(unique(data$batch))
  k <- substr(model, nchar(model), nchar(model)) %>%
    as.integer
  is_SB <- substr(model, 1, 2) == "SB"
  if(nrow(data) > 0){
    number_batches <- ifelse(is_SB, 1, number_batches)
  }
  number_obs <- nrow(data)
  tab <- tibble(model=model,
                k=k,
                number_batches=as.integer(number_batches),
                number_obs=as.integer(number_obs))
  tab
}


Trios <- function(model="TBP3",
                  data=modelData(),
                  specs=trio_spec(model, data),
                  iter=1000L,
                  burnin=500L,
                  thin=1L,
                  nStarts=4L,
                  hp=Hyperparameters(k=specs$k),
                  mp=McmcParams(iter=iter, thin=thin,
                                burnin=burnin,
                                nStarts=nStarts),
                  parameters=modelParameters(mp=mp, hp=hp),
                  chains=mcmc_chainsP(specs, parameters),
                  current_values=modelValuesP(specs, data, hp),
                  summaries=modelSummariesP(specs),
                  flags=modelFlags()){
  if(nrow(data) > 0){
    data$batch <- 1L
  }
  model <- new("Trios",
               data=data,
               specs=specs,
               parameters=parameters,
               chains=chains,
               current_values=current_values,
               summaries=summaries,
               flags=flags)
  model
}
