#include "miscfunctions.h"

using namespace Rcpp;
// Function to simulate from dirichlet distribution
void rdirichlet(Rcpp::NumericVector a, Rcpp::NumericVector pr) {
    double sample[a.size()];
    double sample_sum = 0;
    for(int i=0; i<a.size(); i++) {
        sample[i] = as<double>(rgamma(1, a[i], 1));
        sample_sum += sample[i];
    }
    for(int i = 0; i<a.size(); i++) {
        pr[i] = sample[i] / sample_sum ;
    }
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

// Simulate from mutlivariate normal distribution.
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// Distribution function for skew normal
NumericVector dsn(NumericVector r, double xi, double omega, double alpha) {
    NumericVector z, logN, logS, logPDF;
    z = (r - xi)/omega;
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

