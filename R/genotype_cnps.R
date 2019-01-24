genotypes <- function(se){
  assays(se)[["GT"]]
}

bafs <- function(se){
  assays(se)[["baf"]]
}

shapeParams <- function(){
  m <- matrix(c(1, 1,
                10, 1,
                1, 10,
                30, 30,
                10, 20,
                20, 10,
                30, 10,
                10, 30), ncol=2,
              byrow=TRUE)
  rownames(m) <- c("null",
                   "allB",
                   "zeroB",
                   "balanced",
                   "one-thirdB",
                   "two-thirdsB",
                   "one-fourthB",
                   "three-fourthsB")
  m
}

d_beta <- function(params){
  function(x, ...) dbeta(x, params[1], params[2], ...)
}
null <- function(){
  params <- shapeParams()["null", ]
  d_beta(params)
}
aa <- function(){
  params <- shapeParams()["zeroB", ]
  d_beta(params)
}
bb <- function(){
  params <- shapeParams()["allB", ]
  d_beta(params)
}
ab <- function(){
  params <- shapeParams()["balanced", ]
  d_beta(params)
}
aab <- function(){
  params <- shapeParams()["one-thirdB", ]
  d_beta(params)
}
abb <- function(){
  params <- shapeParams()["two-thirdsB", ]
  d_beta(params)
}
aaab <- function(){
  params <- shapeParams()["one-fourthB", ]
  d_beta(params)
}
abbb <- function(){
  params <- shapeParams()["three-fourthsB", ]
  d_beta(params)
}

.p_b <- function(gt){
  gt <- factor(gt, levels=1:3)
  tab <- table(gt)
  p.bb <- tab[3]/length(gt)
  p.ab <- tab[2]/length(gt)
  p.b <- p.bb + 1/2*p.ab
  as.numeric(p.b)
}

p_b <- function(gt){
  as.numeric(apply(gt, 1, .p_b))
}

d0 <- null()

## p.b = population frequency of B allele
## p.b should be the same length as the x-vector
## x: vector of B allele frequencies
d1 <- function(x, p.b){
  daa <- aa()
  dbb <- bb()
  p.outlier <- 1e-4
  p.outlier + (1-p.outlier)*(
    (1-p.b)*daa(x) + p.b*dbb(x)
  )
}

d2 <- function(x, p.b){
  daa <- aa()
  dab <- ab()
  dbb <- bb()
  p.outlier <- 1e-4
  p.outlier + (1-p.outlier)*(
    (1-p.b)^2*daa(x) + p.b^2*dbb(x) +
    2*p.b*(1-p.b)*dab(x)
  )
}

d3 <- function(x, p.b){
  p.a <- 1-p.b
  daa <- aa()
  daab <- aab()
  dabb <- abb()
  dbb <- bb()
  p.outlier <- 1e-4
  p.outlier + (1-p.outlier)*(
    (1-p.b)^3*daa(x) + choose(3, 1)*p.a^2*p.b*daab(x) +
    choose(3, 1)*p.a*p.b^2*dabb(x) + p.b^3*dbb(x)
  )
}

d4 <- function(x, p.b){
  p.a <- 1-p.b
  daaaa <- aa()
  daaab <- aaab()
  daabb <- ab()
  dabbb <- abbb()
  dbbbb <- bb()
  p.outlier <- 1e-4
  p.outlier + (1-p.outlier)*(
    (1-p.b)^4*daaaa(x) + choose(4, 1)*p.a^3*p.b*daaab(x) +
    choose(4, 2)*p.a^2*p.b^2*daabb(x) +
    choose(4, 1)*p.a*p.b^3*dabbb(x) +
    p.b^4*dbbbb(x)
  )
}

#' @param snpdat SummarizedExperiment with gt and baf assays
#' @param cn.model a CopyNumberModel
mixtureProbs <- function(snpdat, cn.model){
  snpdat2 <- snpdat
  p.b <- p_b(genotypes(snpdat2))
  B <- bafs(snpdat2)
  keep <- rowMeans(is.na(B)) < 0.1
  B <- B[keep, , drop=FALSE]
  p.b <- p.b[ keep ]
  baf <- NULL
  B <- B %>%
    as.tibble %>%
    mutate(p.b=p.b) %>%
    gather("id", "baf", -p.b) %>%
    mutate(p0=d0(baf),
           p1=d1(baf, p.b),
           p2=d2(baf, p.b),
           p3=d3(baf, p.b),
           p4=d4(baf, p.b)) %>%
    filter(!is.na(baf))
  B
}

pBaf <- function(snpdat, cn.model){
  pz <- probCopyNumber(cn.model) %>%
    as.tibble %>%
    set_colnames(paste0("prob_cn", unique(mapping(cn.model)))) %>%
    mutate(id=colnames(snpdat))
  ##pwr <- 1/10
  pwr <- 1
  p0 <- p1 <- p2 <- p3 <- p4 <- baf <- id <- NULL
  B <- mixtureProbs(snpdat, cn.model) %>%
    group_by(id) %>%
    summarize(baf=mean(baf, na.rm=TRUE),
              p0=prod(p0, na.rm=TRUE)^pwr,
              p1=prod(p1, na.rm=TRUE)^pwr,
              p2=prod(p2, na.rm=TRUE)^pwr,
              p3=prod(p3, na.rm=TRUE)^pwr,
              p4=prod(p4, na.rm=TRUE)^pwr) %>%
    left_join(pz, by="id")
  B
}

.pBAF_012 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p0 <- p1 <- p2 <- NULL
  prob_cn0 <- prob_cn1 <- prob_cn2 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("baf", "p0", "p1", "p2", cols)) %>%
    mutate(prob=prob_cn0*p0 + prob_cn1*p1 + prob_cn2*p2)
}

pBAF_012 <- function(snpdat, cn.model){
  B <- .pBAF_012(snpdat, cn.model)
}

pBAF_0123 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p0 <- p1 <- p2 <- p3 <- NULL
  prob_cn0 <- prob_cn1 <- prob_cn2 <- prob_cn3 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("p0", "p1", "p2", "p3", cols)) %>%
    mutate(prob=prob_cn0*p0 + prob_cn1*p1 + prob_cn2*p2 + prob_cn3*p3)
}

pBAF_01234 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p0 <- p1 <- p2 <- p3 <- p4 <- NULL
  prob_cn0 <- prob_cn1 <- prob_cn2 <- prob_cn3 <- prob_cn4 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("p0", "p1", "p2", "p3", "p4", cols)) %>%
    mutate(prob=prob_cn0*p0 + prob_cn1*p1 + prob_cn2*p2 + prob_cn3*p3 + prob_cn4*p4)
}

pBAF_123 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p0 <- p1 <- p2 <- p3 <- NULL
  prob_cn0 <- prob_cn1 <- prob_cn2 <- prob_cn3 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("baf", "p1", "p2", "p3", cols)) %>%
    mutate(prob=prob_cn1*p1 + prob_cn2*p2 + prob_cn3*p3)
}

.pBAF_12 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p0 <- p1 <- p2 <- p3 <- NULL
  prob_cn0 <- prob_cn1 <- prob_cn2 <- prob_cn3 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("baf", "p1", "p2", cols)) %>%
    mutate(prob=prob_cn1*p1 + prob_cn2*p2)
}

pBAF_12 <- function(snpdat, cn.model){
  B <- .pBAF_12(snpdat, cn.model)
}

pBAF_2 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p2 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select("baf", "p2", cols) %>%
    mutate(prob=p2)
}

pBAF_234 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p2 <- p3 <- p4 <- NULL
  prob_cn2 <- prob_cn3 <- prob_cn4 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("baf", "p2", "p3", "p4", cols)) %>%
    mutate(prob=prob_cn2*p2 + prob_cn3*p3 + prob_cn4*p4)
}

pBAF_23 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p2 <- p3 <- p4 <- NULL
  prob_cn2 <- prob_cn3 <- prob_cn4 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("baf", "p2", "p3", cols)) %>%
    mutate(prob=prob_cn2*p2 + prob_cn3*p3)
}

pBAF_02 <- function(snpdat, cn.model){
  cols <- paste0("prob_cn", unique(mapping(cn.model)))
  p0 <- p2 <- p3 <- p4 <- NULL
  prob_cn0 <- prob_cn2 <- prob_cn3 <- prob_cn4 <- NULL
  B <- pBaf(snpdat, cn.model) %>%
    select(c("baf", "p0", "p2", cols)) %>%
    mutate(prob=prob_cn0*p0 + prob_cn2*p2)
}

cmap <- function(model){
  paste(unique(mapping(model)), collapse="")
}

.modelProb <- function(cn.model, snpdata){
  cn_map <- paste0("m_", cmap(cn.model))
##  cn.model <- specs(cn.model)$cn.model
##  cn.model <- strsplit(cn.model, ",") %>%
##    unlist
##  cn.model <- paste(unique(cn.model), collapse="")
##  ##cn_map <- paste0("m_", cmap(cn.model))
##  cn_map <- paste0("m_", cn.model)
  fun <- switch(cn_map,
                m_012=pBAF_012,
                m_021=pBAF_012,
                m_02=pBAF_02,
                m_0123=pBAF_0123,
                m_01234=pBAF_01234,
                m_123=pBAF_123,
                m_12=pBAF_12,
                m_234=pBAF_234,
                m_23=pBAF_23,
                m_2=pBAF_2)
  B <- fun(snpdata, cn.model)
  B
}

#' @export
modelProb <- function(cn.model, snpdata){
  B <- .modelProb(cn.model, snpdata)
  marginal_prob <- sum(log(B$prob))
  marginal_prob
}

candidateModels <- function(cn.model){
  if(k(cn.model) == 4){
    candidate_models <- list(c(0, 1, 2, 2),
                             c(0, 2, 2, 2),
                             c(2, 2, 2, 2),
                             c(0, 1, 1, 2),
                             c(0, 1, 2, 3),
                             c(1, 2, 3, 3)) %>%
      lapply(as.character)
  }
  if(k(cn.model) == 3){
    candidate_models <- list(c(0, 1, 2),
                             ## hemizygous component can not be
                             ## distinguished
                             c(0, 2, 2),
                             ##c(1, 1, 2),
                             c(1, 2, 2),
                             c(2, 2, 2),
                             c(1, 2, 3),
                             c(2, 3, 4),
                             c(2, 3, 3)) %>%
      lapply(as.character)
  }
  if(k(cn.model) == 2){
    candidate_models <- list(c(0, 2),
                             c(1, 2),
                             c(2, 3),
                             c(2, 2)) %>%
      lapply(as.character)
  }
  if(k(cn.model) == 1){
    candidate_models <- list("2")
  }
  model.list <- list()
  for(i in seq_along(candidate_models)){
    m <- cn.model
    mapping(m) <- candidate_models[[i]]
    model.list[[i]] <- m
  }
  model.list
}

ccmap <- function(model){
  paste(mapping(model), collapse="")
}

#' Calculate the likelihood of the observed B allele frequencies for a given copy number model
#'
#' @param cn.model a copy number model
#' @param snpdata a \code{SummarizedExperiment} with assay element `GT` containing an integer coding (1, 2, and 3) for the generic genotypes AA, AB, and BB, respectively.
#' @export
#' @return a list. The first element is the log likeihood and the second element is the copy number model that maximized the likelihood.
bafLikelihood <- function(cn.model, snpdata){
  g <- genotypes(snpdata) %>%
    as.integer
  g <- g[!is.na(g)]
  if(!all(g %in% 1:3)){
    stop("Genotypes are assumed to be integers with values 1, 2, or 3.")
  }
  model.list <- candidateModels(cn.model)
  if(length(model.list) > 1){
    logprobs <- sapply(model.list, modelProb,
                       snpdata=snpdata)
    m <- model.list[[which.max(logprobs)]]
  } else {
    logprobs <- setNames(0, "2")
  }
  model.names <- sapply(model.list, ccmap)
  loglik <- tibble(model=model.names,
                   loglik=logprobs)
  model <- model.list[[which.max(logprobs)]]
  list(loglik=loglik, model=model)
}
