#' Create simulated batch data for testing.
#'
#' @param N number of observations
#' @param p a vector indicating probability of membership to each component
#' @param theta a matrix of means.  Columns are components and rows are batches.
#' @param sds a matrix of standard deviations.  Columns are components and rows are batches.
#' @param batch a vector of labels indication from which batch each simulation should come from
#' @param zz a vector indicating latent variable membership. Can be omitted.
#' @param df length-1 numeric vector for the t-distribution degrees of freedom
#' @return An object of class 'MultiBatchModel'
#' @examples
#' k <- 3
#' nbatch <- 3
#' means <- matrix(c(-1.2, -1.0, -0.8,
#'                   -0.2, 0, 0.2,
#'                   0.8, 1, 1.2), nbatch, k, byrow=FALSE)
#' sds <- matrix(0.1, nbatch, k)
#' N <- 1500
#' truth <- simulateBatchData(N=N,
#'                            batch=rep(letters[1:3], length.out=N),
#'                            theta=means,
#'                            sds=sds,
#'                            p=c(1/5, 1/3, 1-1/3-1/5))
#' @export
simulateBatchData <- function(N=2500, p, theta, sds, batch, zz, df=10){
  ## order ys by batch
  if(!is.matrix(p)) p <- matrix(p, nrow(theta), ncol(theta), byrow=TRUE)
  if(ncol(p) != ncol(theta)) stop("length of p must be same as ncol(theta)")
  if(!all(rowSums(p)==1)) stop("elements of p must sum to 1")
  if(missing(batch)) {
    batch <- rep(1L, N)
  } else {
    batch <- as.integer(factor(batch))
  }
  batch <- sort(batch)
  if(missing(zz)) {
    zz <- simulateZ(N, p)
  }
  yy <- rep(NA, N)
  ub <- unique(batch)
  rownames(theta) <- rownames(sds) <- ub
  for(b in ub){
    index <- which(batch==b)
    nn <- length(index)
    cn <- zz[index]
    mu <- theta[b, ]
    s <- sds[b, ]
    yy[index] <- rnorm(nn, mu[cn], s[cn])/sqrt(rchisq(nn, df)/df)
  }
  ##ix <- order(batch)
  ##object <- MultiBatchModel(yy, batch=batch, k=ncol(theta))
  object <- MultiBatchModel2(dat=yy, batches=batch,
                             hpList(k=ncol(theta))[["MB"]])
  z(object) <- as.integer(factor(zz))
  ##
  ## Must initialize the slots independently (must point to different
  ## locations in memory)
  ##
  theta(object) <- computeMeans(object)
  dataMean(object) <- computeMeans(object)
  mu(object) <- colMeans(dataMean(object))
  sigma2(object) <- computeVars(object)
  dataPrec(object) <- 1/computeVars(object)
  counts <- tableBatchZ(object)
  P <- counts/rowSums(counts)
  p(object) <- P
  log_lik(object) <- computeLoglik(object)
  object
}

#' Create simulated data for testing.
#'
#' @param N number of observations
#' @param p a vector indicating probability of membership to each component
#' @param theta a vector of means, one per component
#' @param sds a vector of standard deviations, one per component
#' @param df length-1 numeric vector for the t-distribution degrees of freedom
#' @return An object of class 'SingleBatchModel'
#' @examples
#' truth <- simulateData(N=2500, p=rep(1/3, 3),
#'                       theta=c(-1, 0, 1),
#'                       sds=rep(0.1, 3))
#' @export
simulateData <- function(N, p, theta, sds, df=10){
  theta <- matrix(theta, nrow=1)
  sds <- matrix(sds, nrow=1)
  p <- matrix(p, nrow=1)
  object <- simulateBatchData(N=N,
                              p=p,
                              theta=theta,
                              sds=sds,
                              df=df,
                              batch=rep(1L, N))
  return(object)
##  zz <- simulateZ(N, p)
##  y <- rnorm(N, theta[zz], sds[zz]/sqrt(rchisq(N, df)/df))
##  ## y <- rnorm(N, theta[zz], sds[zz])
##  ##object <- SingleBatchModel(data=y, k=length(theta))
##  object <- SingleBatchModel2(dat=y, hp=hpList(k=length(theta))[["SB"]])
##  z(object) <- as.integer(factor(zz, levels=unique(sort(zz))))
##  p(object) <- p
##  theta(object) <- matrix(as.numeric(sapply(split(y(object), z(object)), mean)),
##                          nrow=1)
##  sigma2(object) <- matrix(as.numeric(sapply(split(y(object), z(object)), var)),
##                           nrow=1)
##  p(object) <- as.numeric(sapply(split(y(object),
##                                       z(object)),
##                                 length)/length(z(object)))
##  mu(object) <- mean(theta(object))
##  tau2(object) <- var(theta(object))
##  log_lik(object) <- computeLoglik(object)
##  logPrior(object) <- computePrior(object)
##  object
}

simulateZ <- function(N, p){
  P <- rdirichlet(N, p[1, ])
  cumP <- t(apply(P, 1, cumsum))
  u <- runif(N)
  zz <- rep(NA, N)
  zz[u < cumP[, 1]] <- 1
  k <- 2
  while(k <= ncol(P)){
    zz[u < cumP[, k] & u >= cumP[, k-1]] <- k
    k <- k+1
  }
  zz
}

cumProbs <- function(p, k){
  pcum <- list()
  cols <- 2:(k-1)
  for(j in seq_along(cols)){
    g <- cols[j]
    pcum[[j]] <- rowSums(p[, 1:g, drop=FALSE])
  }
  pcum2 <- cbind(p[, 1], do.call(cbind, pcum), 1)
}

simulateProbeLevel <- function(samples=2000,
                               cnvs=200, K=4,
                               probes=10,
                               arguments, qual="easy") {
  l <- probes
  sl.good <- arguments$sl.good
  sl.bad <- arguments$sl.bad
  prbias <- arguments$prbias
  n <- arguments$n
  prvar <- arguments$prvar

  Z <- array(data=NA, dim=c(cnvs,samples,K))
  Z[,,1] <- 1L
  ## Simulate copy numbers from multinomial distribution according to HWE
  theta = 30
  alpha = rbeta(200, 13.3, 13.3)
  ## probability of starting with homozygous,hemizygous,diploid
  ##
  for(j in 2:K) {
    k = j
    i = 0:(k-1)
    lambda = matrix(NA, cnvs, k)
    for(s in 1:cnvs) {
      w = choose(k-1, i) * alpha[s]^i * (1 - alpha[s])^((k-1) - i)
      lambda[s,] <- rdirichlet(1, w * theta)
    }
    ##        w = choose(k-1, i) * alpha^i * (1 - alpha)^((k-1) - i)
    ##        lambda <- rdirichlet(cnvs, w * theta)
    Z[,,j] <- rMultinom(lambda, samples)
    ##randomly shift multinomial sample to begin with 0, 1, or 2 copy number.
    ##Z[,,j] + sample(c(-1,0,1), 1)
  }

  A <- array(data=NA, dim=c(samples, l, cnvs, K))
  Pb  <- matrix(NA, cnvs, l)
  Pv <- matrix(NA, cnvs, l)
  for(i in 1:nrow(Pb)) Pb[i,] <- rnorm(l, 0, prbias)
  for(i in 1:nrow(Pb)) Pv[i,]   <- rgamma(l, shape = prvar[1], scale = prvar[2])

  corrupted <- sample(samples, ceiling(0.004*samples))
  for(k in 1:K) {
    for(j in 1:cnvs) {
      for(i in 1:samples) {
        p.mean <- c(rep(sl.good,3),rep(sl.bad,7)) * Z[j,i,k] + Pb[j,]
        p.sd <- n + Pv[j,]

        A[i,,j,k] <- rnorm(l,  p.mean, p.sd)
      }
      ## null measurements for medium and bad quality
      if(qual == "medium" || qual == "hard") {
        v <- range(A[,,j,k])
        for(ind in corrupted) A[ind,,j,k] <- runif(l, v[1], v[2])
      }
    }
  }
  return(list("measurements"=A, "assignments"=Z))
}


hapmapTruth <- function(avg.lrr, simulate_truth){
  ##hapmap.dir <- "PancCnvsData2/inst/extdata/hapmap"
  hapmap.dir <- system.file("extdata/hapmap", package="PancCnvsData2")
  if(!simulate_truth){
    true.z <- readRDS(file.path(hapmap.dir, "cnp57_truth.rds"))
    return(true.z)
  }
  mp <- McmcParams(iter = 1000, burnin = 1000, nStarts = 4, thin = 4)
  model <- SB(dat=avg.lrr, hp=hpList(k=3)[["SB"]], mp=mp)
  ##
  ## idea is that any method would get this right without batch effects,
  ## so we define the truth based on cnpbayes or cnvcall fit
  ##
  truth <- posteriorSimulation(model)
  true.z <- z(truth)
  true.z <- true.z - 1
  ##saveRDS(true.z, file=file.path(hapmap.dir, "hapmap_truth.rds"))
  saveRDS(true.z, file=file.path(hapmap.dir, "cnp57_truth.rds"))
  true.z
}

#' @export
simulateBatchEffect <- function(hap.list,
                                experiment){
  scale <- experiment$scale[[1]]
  shift <- experiment$batch.effect[[1]]
  r <- hap.list[["r"]] %>%
    group_by(plate, id) %>%
    summarize(oned=mean(lrr, na.rm=TRUE)) %>%
    left_join(hap.list[["truth"]], by="id")
  cn.stats <- ungroup(r) %>%
    group_by(cn) %>%
    summarize(mu=median(oned),
              sd=sd(oned))
  snr <- diff(cn.stats$mu)/cn.stats$sd[1:2]
  ## the first component is ~9 sds from componenent 2
  ## the second component is ~6 sds from component 3
  ##
  ## Batch effect
  ##
  nplates <- length(unique(r$plate))
  batches <- tibble(plate=unique(r$plate),
                    batch=sample(c(0, 1),
                                 nplates,
                                 replace=TRUE,
                                 prob=c(0.5, 0.5)))
  r2 <- left_join(r, batches, by="plate") %>%
    ungroup() %>%
    mutate(mu=cn.stats$mu[ as.integer(cn)] ,
           delta=ifelse(batch==0, 0, shift),
           centered=oned-mu) %>%
    mutate(rescaled=centered*scale + mu,
           raw=oned,
           oned=rescaled + rnorm(nrow(.), delta, 0.02)) %>%
    select(c(plate, id, oned, raw, cn, batch))
  b <- hap.list[["b"]] %>%
    select(-probe)
  g <- hap.list[["g"]] %>%
    select(-probe)
  r3 <- left_join(r2, b, by=c("id", "plate")) %>%
    left_join(g, by=c("id", "plate"))
  r3
}
  ##number.plates <- sum(batches == 0)
  ##  n.plate <- sapply(lrr.batch, ncol)
  ##  plates <- factor(rep(names(lrr.batch), n.plate),
  ##                   levels=names(lrr.batch))
  ##  z.list <- split(true.z, plates)
  ##medians2 <- medians[true.z + 1]
  ##medians.list <- split(medians2, plates)
  ##delta <- ifelse(batch.id == 0, shift, 0)
  ##lrr.sim <- lrr.batch
  ##
  ##  Below is the random part of the simulation -- adding shift and scale
  ##
##  for(j in seq_along(batch.id)){
##    r <- lrr.batch[[j]]
##    ## centering by z.meds will not remove copy number effect
##    z.meds <- matrix(medians.list[[j]],
##                     nrow(r), ncol(r), byrow=TRUE)
##    centered <- r - z.meds
##    rr <- (centered)*scale + z.meds
##    if(batch.id[j] == 0) {
##      ##
##      ## On average, each observation is shifted by amount [shift]
##      ##
##      ##rr <- rr + rnorm(nrow(rr) * ncol(rr), shift, 0.01)
##      rr <- rr + shift
##    }
##    lrr.sim[[j]] <- rr
##  }
##  lrr.sim



getHapmapIds <- function(hap.list){
  ids <- lapply(hap.list[[1]], colnames) %>%
    lapply("[", -c(1, 2)) %>%
    unlist
  ids
}

##hapmapData <- function(experiment, hap.list){
##  ids <- getHapmapIds(hap.list)
##  shift <- experiment$batch.effect
##  scale <- experiment$scale
##  seed <- experiment$seed
##  set.seed(seed)
##  simulate_truth <- experiment$id == "001"
##  ## plates will be in different batches even if the overall modes
##  ## look the same. This is because the representation of copy number in the different plates is different.  homdel are rare in plates 12 and 14
##  sim.data <- hapmapSimulation2(shift, scale,
##                                simulate_truth=simulate_truth,
##                                thr=0.0001)
##  return(sim.data)
##  full.data <- sim.data$full.data %>%
##    mutate(id=ids) %>%
##    arrange(batch_index)
##
##  cnp57 <- readRDS("PancCnvsData2/inst/extdata/hapmap/cnp_57.rds")
##  probe_id <- cnp57[["r"]][[1]]$probeset_id
##  sample_ids <- lapply(cnp57[["r"]], colnames) %>%
##    lapply("[", -c(1,2)) %>%
##    unlist
##  lrr.sim <- sim.data$lrr.sim
##  xx <- t(do.call(cbind, lrr.sim))
##  dimnames(xx) <- list(sample_ids, probe_id)
##  list(full.data=full.data, probe.level=xx)
##}

hapmapAvgLrr <- function(dat, min_plate=0){
  ##hapmap.dir <- "PancCnvsData2/inst/extdata/hapmap"
  ##lrr.batch <- readRDS(file.path(hapmap.dir, "lrr_roi_1.rds"))
  cnp_57 <- dat
  ##cnp_57 <- readRDS(file.path(hapmap.dir, "cnp_57.rds"))
  matFun <- function(x) as.matrix(x[, -c(1, 2)])
  plates <- rep(names(cnp_57[["r"]]), sapply(cnp_57[["r"]], ncol)-2)
  np <- table(plates)
  np <- np[ np > min_plate ]
  lrr.batch <- cnp_57[["r"]][names(np)] %>%
    lapply(matFun)
  ##lrr.batch <- lapply(lrr.batch, function(x) x/1000)
  ##
  ## Simulate batches
  ##
  ##batch.id <- c(rep(0,8), rep(1, 8))
  batch.id <- sample(c(0, 1), length(lrr.batch), replace=TRUE,
                     prob=c(0.5, 0.5))
  avg.lrr <- unlist(lapply(lrr.batch, colMeans, na.rm=TRUE))
  list(lrr.batch=lrr.batch, avg.lrr=avg.lrr, batch.id=batch.id)
}

hapmapSimulation2 <- function(hapmap.data, 
                              shift, scale,
                              simulate_truth=FALSE,
                              thr=0.01, min_plate=0){
  avglrrs <- unlist(lapply(lrr.sim, colMeans, na.rm=TRUE))
  ncols <- sapply(lrr.batch, ncol)
  plates <- rep(names(lrr.batch), ncols)
  plate.index <- plates %>%
    factor(levels=unique(.)) %>%
    as.integer
  true_batch <- rep(batch.id, sapply(lrr.batch, ncol)) %>%
    as.character
  dat <- tibble(id=seq_along(avglrrs),
                oned=avglrrs,
                provisional_batch=plate.index,
                batch=1L,
                true_batch=true_batch)
  dat
}

hapmapData <- function(cnp_57, true.z, experiment){
  hap.list <- hapmapAvgLrr(cnp_57)
  hap.list[["true.z"]] <- true.z
  lrr.list <- simulateBatchEffect(hap.list, experiment)
  avglrrs <- unlist(lapply(lrr.list, colMeans, na.rm=TRUE))
  ids <- strsplit(names(avglrrs), "\\.") %>%
    sapply("[", 2)
  lrr.batch <- hap.list[["lrr.batch"]]
  batch.id <- hap.list[["batch.id"]]
  ncols <- sapply(lrr.batch, ncol)
  plates <- rep(names(lrr.batch), ncols)
  plate.index <- plates %>%
    factor(levels=unique(.)) %>%
    as.integer
  true_batch <- rep(batch.id, sapply(lrr.batch, ncol)) %>%
    as.character
  dat <- tibble(id=ids,
                oned=avglrrs,
                provisional_batch=plate.index,
                batch=1L,
                true_batch=true_batch)
  dat
}

#' @export
hapmapSummarizedExperiment <- function(hapmap, gr){
  B <- hapmap[["b"]] 
  B2 <- select(B, c(id, baf)) %>%
    spread("id", "baf") %>%
    as.matrix
  rownames(B2) <- unique(B$probe)
  ##  B <- cnp57[["b"]] %>%
  ##    lapply(function(x) as.matrix(x[, -c(1, 2)])) %>%
  ##    do.call(cbind, .)
  ##rownames(B) <- cnp57$b[[1]]$probeset_id
  G <- hapmap[["g"]]
  G2 <- select(G, c(id, genotype)) %>%
    spread("id", "genotype") %>%
    as.matrix
  rownames(G2) <- unique(G$probe)
  ##  G <- cnp57[["g"]] %>%
  ##    lapply(function(x) as.matrix(x[, -c(1, 2)])) %>%{
  ##      do.call(cbind, .)
  ##    }
  ##  G[ !G %in% 1:3 ] <- NA
  ##  rownames(G) <- rownames(B)
  ##rr <- rowRanges(cnp_se)["CNP_057"]
  ##names(gr) <- cnp57$b[[1]]$probeset_id
  names(gr) <- rownames(B2)
  snpdat <- SummarizedExperiment(assays=SimpleList(baf=B2,
                                                   GT=G2),
                                 rowRanges=gr)
  snpdat
}

hapmapExperiment <- function(){
  batch.effect <- seq(0, 0.5, by=0.1)
  sep <- batch.effect
  scale <- seq(1, 2, by=0.25)
  paramgrid <- expand.grid(batch.effect, scale)
  nseeds <- 10
  experiment <- tibble(batch.effect=rep(paramgrid[, 1], each=nseeds),
                       scale=rep(paramgrid[, 2], each=nseeds)) %>%
    mutate(seed=as.integer(round(runif(nrow(.), 1, 100000), 0))) %>%
    mutate(id = str_pad(seq_len(nrow(.)), 3, pad="0"))
  experiment
}

#' @export
performanceStats <- function(truth, cn){
  k <- length(table(cn))
  tib <- tibble(truth=as.character(truth),
                cn=as.character(cn))
  stats <- tib %>%
    group_by(truth) %>%
    summarize(n=n(),
              correct=sum(cn == truth, na.rm=TRUE),
              incorrect=sum(cn != truth, na.rm=TRUE)) %>%
    mutate(number_components=k)
  TP <- sum(stats$correct[1:2])
  P <- sum(stats$n[1:2])
  N <- stats$n[3]
  TN <- stats$correct[3]
  sensitivity <- TP/P
  specificity <- TN/N
  list(stats, sensitivity, specificity)
}
