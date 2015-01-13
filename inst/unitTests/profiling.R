##
## example (from test_batch.R)
##
set.seed(1234)
k <- 3
nbatch <- 3
means <- matrix(c(-1.9, -2, -1.85,
                  -0.45, -0.4, -0.35,
                  -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
sds <- matrix(0.3, nbatch, k)
truth <- simulateBatchData(2500,
                           means=means,
                           sds=sds,
                           .batch=rep(letters[1:3], length.out=2500),
                           .alpha=c(5, 500, 1000))

mcmcp <- McmcParams(iter=1000, burnin=1000)

## use the true k
k <- 3
## marginal
params <- ModelParams("marginal", y=y(truth), k=k,
                      batch=rep(letters[1:3], length.out=2500),
                      mcmc.params=mcmcp)
model <- initializeModel(params)
#model <- posteriorSimulation(model, mcmcp)

updateTheta(model, constrain=TRUE)
updateThetaCpp(model, constrain=TRUE)
## Need to profile posteriorSimulation.  
Rprof()
posteriorSimulation(model, mcmcp)



##### Rcpp S4 testing stuff
library(inline)
src <- '
S4 foo(x) ; foo.slot(".Data") = "bar" ; return(foo);
'
fun <- cxxfunction(signature(x="any"), src,
                   plugin="Rcpp")
setClass( "S4ex", contains = "character",
representation( x = "numeric" ) ) 
x <- new( "S4ex", "bla", x = 10 )
fun(x)
str(fun(x))

