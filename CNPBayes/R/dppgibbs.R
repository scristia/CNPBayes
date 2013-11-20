dppgibbs <-
function(r, ## data
		     H, ## max number of clusters in block DP
		     alpha=1, ## concentration parameter
		     mu0=0, ## prior mean of theta for all clusters
		     tau20=0.1, ## prior prec for theta, all clusters
		     a=0.1,
		     b=0.1,
		     S=100
		     ){
	if(missing(H)) stop("H missing")
	n <- length(r) ## number of subjects
	##
	#########
	# Inits #
	#########
	 pi<-ns<-rep(0,H) 	# Mixing weights and number of subjects per cluster
	 v<-rep(1/H,H)		# Conditional weights -- pr(c_i=h|c_i not in l<h)
	 v[H]<-1			# Apply DP truncation to H classes
	 mu<-rep(0,H)		# Cluster-specific means
	 tau<-sigma2<-rep(1,H)	
	 p<-tmp2<-matrix(0,n,H) # p[i,h] = posterior prob that subject i belongs to cluster h
	 
	#########
	# Store #
	#########
	 V<-Mu<-Sigma<-N<-Pi<-matrix(0,S,H)
	 C<-matrix(0,S,n)
	 grid<-seq(min(y),max(y),length=500)
	 Y<-array(0,dim=c(S,length(grid),H))

	#########
	# GIBBS #
	#########
	for (i in 1:S) {
	 # Update c, cluster indicator
	  cumv<-cumprod(1-v)
	  pi[1]<-v[1]
	  for (h in 2:H) pi[h]<-v[h]*cumv[h-1]
	  for (h in 1:H) tmp2[,h]<-pi[h]*dnorm(y,mu[h],sqrt(sigma2[h]))
	  p<-tmp2/apply(tmp2,1,sum)
	  C[i,]<-c<-rMultinom(p,1)
	  Pi[i,]<-pi
	  for (h in 1:H) ns[h]<-length(c[c==h])  # Must allow zeros for empty clusters

	 # Update v
	  for (h in 1:(H-1)) v[h]<-rbeta(1,1+ns[h],alpha+sum(ns[(h+1):H]))
	  V[i,]<-v
	   
	 # Update mu and sigma2 and Yhat (density estimate)
	 for (h in 1:H) {
	   var<-1/(tau20+tau[h]*ns[h])
	   m<-var*(tau20*mu0+tau[h]*sum(y[c==h])) 
	   Mu[i,h]<-mu[h]<-rnorm(1,m,sqrt(var))
	   tau[h]<-rgamma(1,a+ns[h]/2,b+sum((y[c==h]-mu[h])^2)/2)
	   Sigma[i,h]<-sigma2[h]<-1/tau[h]
	   Y[i,,h]<-pi[h]*dnorm(grid,mu[h],sqrt(sigma2[h]))
	 }
	 N[i,]<-ns  		# Number per cluster
	 if (i%%100==0) print(i)

	}
	list(P=Pi, means=Mu, precs=1/Sigma, Z=C, N=N)
}
