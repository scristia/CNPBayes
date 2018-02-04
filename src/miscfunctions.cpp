#ifndef _miscfunctions_H
#define _miscfunctions_H

#include <Rcpp.h>


using namespace Rcpp;
// FUNCTIONS FOR ACCESSING HYPERPARAMETERS
// [[Rcpp::export]]
int getK(Rcpp::S4 hyperparams) {
  int k = hyperparams.slot("k");
  return k;
}

// [[Rcpp::export]]
double getDf(Rcpp::S4 hyperparams) {
  double df = hyperparams.slot("dfr");
  return df;
}

// [[Rcpp::export]]
Rcpp::IntegerVector uniqueBatch(Rcpp::IntegerVector x) {
  IntegerVector tmp = unique(x) ;
  IntegerVector b = clone(tmp) ;
  std::sort(b.begin(), b.end()) ;
  return b ;
}

// [[Rcpp::export]]
Rcpp::IntegerVector tableZ(int K, Rcpp::IntegerVector z){
  Rcpp::IntegerVector nn(K) ;
  for(int k = 0; k < K; k++){
    nn[k] = sum(z == (k+1)) ;
  }
  return nn ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix tableBatchZ(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  int K = getK(model.slot("hyperparams")) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  IntegerVector z = model.slot("z") ;
  NumericMatrix nn(B, K) ;
  for(int j = 0; j < B; ++j){
    for(int k = 0; k < K; k++){
      nn(j, k) = sum((z == (k+1)) & (batch == ub[j]));
    }
  }
  return nn ;
}

// Accessors
Rcpp::IntegerVector getZ(Rcpp::S4 model) {
    IntegerVector z = model.slot("z");
    return z;
}

Rcpp::NumericVector getData(Rcpp::S4 model) {
    NumericVector y = model.slot("data");
    return y;
}



Rcpp::NumericVector getMu(Rcpp::S4 hyperparams) {
  NumericVector mu = hyperparams.slot("mu");
  return mu;
}

Rcpp::NumericVector getTau2(Rcpp::S4 hyperparams) {
    NumericVector tau2 = hyperparams.slot("tau2");
    return tau2;
}

Rcpp::IntegerVector getAlpha(Rcpp::S4 hyperparams) {
    IntegerVector alpha = hyperparams.slot("alpha");
    return alpha;
}

Rcpp::LogicalVector nonZeroCopynumber(Rcpp::IntegerVector z) {
//nonZeroCopynumber <- function(object) as.integer(as.integer(z(object)) > 1)
 LogicalVector nz = z > 1;
 return nz;
}



using namespace Rcpp;
// Function to simulate from dirichlet distribution
void rdirichlet(Rcpp::NumericVector a, Rcpp::NumericVector pr) {
  Rcpp::NumericVector sample(a.size());
  double sample_sum = 0;
  for(int i=0; i<a.size(); i++) {
    sample[i] = as<double>(rgamma(1, a[i], 1));
    sample_sum += sample[i];
  }
  for(int i = 0; i<a.size(); i++) {
    pr[i] = sample[i] / sample_sum ;
  }
}

// generate multinomial random variables with varying probabilities
// [[Rcpp::export]]
Rcpp::IntegerMatrix rMultinom(Rcpp::NumericMatrix probs, int m) {
    // set RNG
    Rcpp::RNGScope scope;

    // get dimensions of probs
    int n = probs.nrow();
    int k = probs.ncol();
    
    Rcpp::IntegerMatrix ran(n, m);

    Rcpp::NumericVector z(n);
    Rcpp::NumericMatrix U(k, n);
    
    for (int i = 0; i < n; i++) {
        z[i] = Rcpp::sum(probs(i, Rcpp::_));
        Rcpp::NumericVector cumsum_temp = Rcpp::cumsum(probs(i, Rcpp::_));
        U(Rcpp::_, i) = cumsum_temp;
    }

    for (int i = 0; i < m; i++) {
        Rcpp::NumericVector rand = Rcpp::runif(n);
        Rcpp::NumericVector un(k * n);
        int index = 0;

        Rcpp::IntegerMatrix compare(k, n);
        int ind = 0;

        // C++ equivalent of `un <- rep(rand, rep(k, n))`
        for (int a = 0; a < n; a++) {
            std::fill(un.begin() + index, un.begin() + index + k, rand[a]);
            index += k;

            for (int b = 0; b < k; b++) {
                compare(b, a) = un[ind] > U(b, a);
                ind++;
            }
            
            ran(a, i) = Rcpp::sum(compare(Rcpp::_, a)) + 1;
        }
    }

    return ran;
}

// Function for drawing from contrained normal distribution for theta
double cons_normal(double mean, double var, double a, double b) {
    double p = R::pnorm(a, mean, sqrt(var), 1, 0) + as<double>(runif(1)) *
        (R::pnorm(b, mean, sqrt(var), 1, 0) - R::pnorm(a, mean, sqrt(var), 1, 0));
    return R::qnorm(p, mean, sqrt(var), 1, 0);
}

// truncated normal using inverse probability transform.
// This is not stable when endpoint is far from mean.
double trunc_norm(double mean, double sd) {
    double p = R::pnorm(0, mean, sd, 1, 0) + as<double>(runif(1)) *
        ( 1 - R::pnorm(0, mean, sd, 1, 0));
    return R::qnorm(p, mean, sd, 1, 0);
}

// Distribution function for skew normal
Rcpp::NumericVector dsn(NumericVector r, double xi, double omega, double alpha) {
    NumericVector z, logN, logS, logPDF; z = (r - xi)/omega;
    logN = -log(sqrt(2.0 * PI)) - log(omega) - pow(z, 2) / 2.0;
    logS = log(pnorm(alpha * z));
    logPDF = logN + logS - R::pnorm(0.0, 0.0, 1.0, 1, 1);
    return exp(logPDF);
}

// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================
double norm_rs(double a, double b)
{
   double  x;
   x = Rf_rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = norm_rand();
   return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

double half_norm_rs(double a, double b)
{
   double   x;
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling.
// ======================================================================

double unif_rs(double a, double b)
{
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   if(a <= 0.0) xstar = 0.0;
   else xstar = a;
   logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

   x = R::runif(a, b);
   logu = log(R::runif(0.0, 1.0));
   while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
   {
      x = R::runif(a, b);
      logu = log(R::runif(0.0, 1.0));
   }
   return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

double exp_rs(double a, double b)
{
  double  z, u, rate;

//  Rprintf("in exp_rs");
  rate = 1/a;
//1/a

   // Generate a proposal on (0, b-a)
   z = R::rexp(rate);
   while(z > (b-a)) z = R::rexp(rate);
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
      z = R::rexp(rate);
      while(z > (b-a)) z = R::rexp(rate);
      u = R::runif(0.0,1.0);
   }
   return(z+a);
}




// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

double rnorm_trunc (double mu, double sigma, double lower, double upper)
{
int change;
 double a, b;
 double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
 double z, tmp, lograt;

 change = 0;
 a = (lower - mu)/sigma;
 b = (upper - mu)/sigma;

 // First scenario
 if( (a == R_NegInf) || (b == R_PosInf))
   {
     if(a == R_NegInf)
       {
     change = 1;
     a = -b;
     b = R_PosInf;
       }

     // The two possibilities for this scenario
     if(a <= 0.45) z = norm_rs(a, b);
     else z = exp_rs(a, b);
     if(change) z = -z;
   }
 // Second scenario
 else if((a * b) <= 0.0)
   {
     // The two possibilities for this scenario
     if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
       {
     z = norm_rs(a, b);
       }
     else z = unif_rs(a,b);
   }
 // Third scenario
 else
   {
     if(b < 0)
       {
     tmp = b; b = -a; a = -tmp; change = 1;
       }

     lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
     if(lograt <= logt2) z = unif_rs(a,b);
     else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
     else z = exp_rs(a,b);
     if(change) z = -z;
   }
   double output;
   output = sigma*z + mu;
 return (output);
}

// [[Rcpp::export]]
Rcpp::NumericVector dlocScale_t(NumericVector x, double df, double mu, double sigma) {
    double coef = tgamma((df + 1.0)/2.0)/(sigma*sqrt(df*PI)*tgamma(df/2.0));
    NumericVector d = coef*pow(1 + pow((x - mu)/sigma, 2.0)/df, -(df+1.0)/2.0);
    return d;
}


// [[Rcpp::export]]
Rcpp::NumericVector compute_u_sums(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  // IntegerVector nn = model.slot("zfreq") ;
  NumericVector sums( K ) ;
  NumericVector u = model.slot("u") ;
  for(int i = 0; i < n; i++){
    for(int k = 0; k < K; k++){
      if(z[i] == k+1){
        sums[k] += u[i] ;
      }
    }
  }
  return sums ;
}


// [[Rcpp::export]]
Rcpp::NumericVector compute_heavy_sums(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 xmod = clone(object) ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  // IntegerVector nn = model.slot("zfreq") ;
  
  NumericVector sums( K ) ;
  NumericVector u = model.slot("u") ;
  x = x * u ;
  for(int i = 0; i < n; i++){
    for(int k = 0; k < K; k++){
      if(z[i] == k+1){
        sums[k] += x[i] ;
      }
    }
  }
  return sums ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_heavy_means(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  IntegerVector z = model.slot("z") ;
  int K = getK(hypp) ;
  IntegerVector nn = tableZ(K, z) ;
  NumericVector means = compute_heavy_sums(xmod) ;
  for(int k = 0; k < K; k++){
    means[k] = means[k] / nn[k] ;
  }
  return means ;
}



// [[Rcpp::export]]
Rcpp::NumericVector compute_u_sums_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector u = model.slot("u") ;
  int n = u.size() ;

  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  NumericMatrix sums(B, K) ;
  for(int i = 0; i < n; i++){
      for(int b = B; b < B; b++) {
          for(int k = 0; k < K; k++){
              if(z[i] == k+1){
                  sums(b, k) += u[i] ;
              }
          }
      }
  }
  return sums ;
}

// [[Rcpp::export]] Rcpp::NumericVector
Rcpp::NumericMatrix compute_heavy_sums_batch(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 xmod = clone(object) ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;

  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  NumericVector u = model.slot("u") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;

  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  // IntegerVector nn = model.slot("zfreq") ;
  NumericMatrix sums(B, K) ;
  
  x = x * u ;
  for(int i = 0; i < n; i++){
      for(int b = B; b < B; b++) {
          for(int k = 0; k < K; k++){
              if(z[i] == k+1){
                  sums(b, k) += x[i] ;
              }
          }
      }
  }
  return sums ;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_heavy_means_batch(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  IntegerVector z = model.slot("z") ;
  int K = getK(hypp) ;

  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = uniqueBatch(batch) ;
  int B = ub.size() ;
  NumericMatrix nn = tableBatchZ(xmod) ;
  NumericMatrix means = compute_heavy_sums_batch(xmod) ;
  for(int b = B; b < B; b++) {
      for(int k = 0; k < K; k++){
          means(b, k) = means(b, k) / nn(b, k) ;
      }
  }
  return means ;
}

#endif
