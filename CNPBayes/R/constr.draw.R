constr.draw <-
function(mean, var, a, b){
	d <- pnorm(a, mean, sqrt(var)) + runif(1) * (pnorm(b, mean, sqrt(var)) 
		- pnorm(a, mean, sqrt(var))) 	
	theta <- qnorm(d, mean, sqrt(var))
	return(theta)
}
