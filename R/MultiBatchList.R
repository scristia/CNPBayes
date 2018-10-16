#' @include MultiBatch.R
NULL

setClass("MultiBatchList", representation(data="tbl_df",
                                          model_names="character",
                                          parameters="list",
                                          chains="list",
                                          current_values="list",
                                          summaries="list",
                                          flags="list"))

modelSpecs <- function(models, K, data) {
  if(missing(models)) models <- c("SB", "SBP", "MB", "MBP")
  if(missing(K)) K <- 1:4
  models <- lapply(models, paste0, 1:4) %>%
    unlist
  number_batches <- max(1, length(unique(data$batch)))
  k <- substr(models, nchar(models), nchar(models)) %>%
    as.integer
  is_SB <- substr(models, 1, 2) == "SB"
  number_batches <- ifelse(is_SB, 1, number_batches)
  number_obs <- sum(data$is_sampled)
  tab <- tibble(model=models, k=k,
                number_batches=number_batches,
                number_obs=number_obs)
  tab
}

listChains2 <- function(model_specs, parameters){
  S <- iter(parameters[["mp"]])
  B <- model_specs$number_batches
  K <- model_specs$k
  mc.list <- vector("list", nrow(model_specs))
  for(i in seq_len(nrow(model_names))){
    mc.list[[i]] <- initialize_mcmc(K[i], S, B[i])
  }
  mc.list
}

listValues <- function(model_specs){
  
}

MultiBatchList <- function(data=modelData(),
                           K,
                           models,
                           model_specs=modelSpecs(models, K, data),
                           num_models=nrow(model_specs),
                           parameters=modelParameters(),
                           chains=listChains2(num_models, parameters),
                           current_values=listValues(model_specs),
                           summaries=listSummaries(model_specs),
                           flags=listFlags(model_specs)){
  new("MultiBatchList",
      data=data,
      model_names=model_names$model,
      parameters=parameters,
      chains=chains,
      current_values=current_values,
      summaries=summaries,
      flags=flags)
}
