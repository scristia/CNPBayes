set_param_names <- function(x, nm){
  K <- seq_len(ncol(x))
  set_colnames(x, paste0(nm, K))
}

mcmcList <- function(model.list){
  ch.list <- purrr::map(model.list, chains)
  theta.list <- purrr::map(ch.list, theta)
  theta.list <- purrr::map(theta.list, set_param_names, "theta")
  sigma.list <- purrr::map(ch.list, "sigma")
  sigma.list <- map(sigma.list, set_param_names, "sigma")
  p.list <- map(ch.list, "p")
  p.list <- map(p.list, set_param_names, "p")
  half <- floor(nrow(theta.list[[1]])/2)
  first_half <- function(x, half){
    x[seq_len(half), ]
  }
  last_half <- function(x, half){
    i <- (half + 1):(half*2)
    x[i, ]
  }
  theta.list <- c(map(theta.list, first_half, half),
                  map(theta.list, last_half, half))
  sigma.list <- c(map(sigma.list, first_half, half),
                  map(sigma.list, last_half, half))
  p.list <- c(map(p.list, first_half, half),
              map(p.list, last_half, half))
  vars.list <- vector("list", length(p.list))
  for(i in seq_along(vars.list)){
    vars.list[[i]] <- cbind(theta.list[[i]], sigma.list[[i]], p.list[[i]])
  }
  vars.list <- map(vars.list, mcmc)
  mlist <- mcmc.list(vars.list)
  mlist      
 }
