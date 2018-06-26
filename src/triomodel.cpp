// [[Rcpp::depends(RcppArmadillo)]]

//#include <RcppArmadillo.h>

#include "miscfunctions.h" // for rdirichlet
#include "multibatch.h" 
#include <Rmath.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <list>

using namespace Rcpp ;
//using namespace RcppArmadillo;


// [[Rcpp::export]]
Rcpp::CharacterVector family_member(Rcpp::S4 object){
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  Rcpp::DataFrame triodat(model.slot("triodata"));
  IntegerVector batch = model.slot("batch") ;
  int n=batch.size();
  Rcpp::CharacterVector family_member(n);
  Rcpp::NumericVector cn(n);
  family_member = triodat["family_member"] ;
  // cn = triodat["copy_number"] ;
  return family_member;
  // return cn;
}

// [[Rcpp::export]]
Rcpp::NumericVector lookup_mprobs(Rcpp::S4 model, int father, int mother){
  Rcpp::NumericMatrix mprob = model.slot("mprob");
  IntegerVector f = model.slot("father");
  IntegerVector m = model.slot("mother");
  IntegerVector map = model.slot("maplabel");
  int fat2 = father - 1;
  int mot2 = mother - 1;
  int f2 = map[fat2];
  int m2 = map[mot2];
  Rcpp::LogicalVector ind(m.size());
  Rcpp::LogicalVector ind2(f.size());
  ind = f == f2;
  ind2 = m == m2;
  Rcpp::LogicalVector is_parental_cn(ind.size());
  is_parental_cn = ind==TRUE & ind2==TRUE;
  int nr=f.size();
  int j = 0;
  for(int i = 0; i < nr; i++){
    if(is_parental_cn[i] == TRUE){
      j = i ;
      break;
    }
  }
  Rcpp::NumericVector result(mprob.ncol());
  result=mprob(j, _);
  return result;
}

//
// RS notes:Currently, component labels do not necessarily correspond to the integer copy
// number. Update of the component labels for the parents is independent of copy
// number. Updates of the component labels for the offspring depends on how
// components are mapped to copy number.
//
// update_trioPr is the Mendelian prob lookup module

// [[Rcpp::export]]
Rcpp::NumericMatrix update_trioPr(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::DataFrame triodat(model.slot("triodata"));
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector z = model.slot("z");
  CharacterVector fam = family_member(xmod);
  //std::string level = Rcpp::as<std::string>(level_of_species[0]);
  Rcpp::LogicalVector fat_ind(fam.size());
  Rcpp::LogicalVector mat_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    fat_ind[i] = (fam[i] == "f");
    mat_ind[i] = (fam[i] == "m");
  }
  Rcpp::IntegerVector zf = z[fat_ind];
  Rcpp::IntegerVector zm = z[mat_ind];
  int trio_size = zf.size();
  Rcpp::NumericMatrix zo_prob(trio_size, K);
  Rcpp::NumericVector trans_probs(K);
  for (int i = 0; i < trio_size; i++){
    trans_probs = lookup_mprobs(xmod, zf[i], zm[i]);
    zo_prob(i,_) = trans_probs ;
  }
  return zo_prob;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix tableBatchZpar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  int K = getK(model.slot("hyperparams")) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zp = z[!child_ind];
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  NumericMatrix nn(B, K) ;
  for(int j = 0; j < B; ++j){
    for(int k = 0; k < K; k++){
      nn(j, k) = sum((zp == (k+1)) & (batch == ub[j]));
    }
  }
  return nn ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix tableBatchZchd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  int K = getK(model.slot("hyperparams")) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zc = z[child_ind];
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  NumericMatrix nc(B, K) ;
  for(int j = 0; j < B; ++j){
    for(int k = 0; k < K; k++){
      nc(j, k) = sum((zc == (k+1)) & (batch == ub[j]));
    }
  }
  return nc ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_mupar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  double tau2_0 = hypp.slot("tau2.0") ;
  double tau2_0_tilde = 1/tau2_0 ;
  double mu_0 = hypp.slot("mu.0") ;
  
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector tau2_tilde = 1/tau2 ;
  IntegerVector z = model.slot("z") ;
  NumericMatrix theta = model.slot("theta") ;
  // below line not referenced anywhere in function
  // IntegerVector nn = model.slot("zfreq") ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  
  NumericVector tau2_B_tilde(K) ;;
  for(int k = 0; k < K; ++k) tau2_B_tilde[k] = tau2_0_tilde + B*tau2_tilde[k] ;
  
  NumericVector w1(K) ;
  NumericVector w2(K) ;
  for(int k = 0; k < K; ++k){
    w1[k] = tau2_0_tilde/(tau2_0_tilde + B*tau2_tilde[k]) ;
    w2[k] = B*tau2_tilde[k]/(tau2_0_tilde + B*tau2_tilde[k]) ;
  }
  NumericMatrix n_b = tableBatchZpar(model) ;
  NumericVector theta_bar(K) ;
  NumericVector th(K) ;
  
  for(int k = 0; k < K; ++k){
    double n_k = 0.0 ; // number of observations for component k
    double colsumtheta = 0.0;
    for(int i = 0; i < B; ++i){
      colsumtheta += n_b(i, k)*theta(i, k) ;
      n_k += n_b(i, k) ;
    }
    theta_bar[k] = colsumtheta/n_k ;
  }
  NumericVector mu_n(K) ;
  NumericVector mu_new(K) ;
  double post_prec ;
  for(int k=0; k<K; ++k){
    post_prec = sqrt(1.0/tau2_B_tilde[k]) ;
    mu_n[k] = w1[k]*mu_0 + w2[k]*theta_bar[k] ;
    mu_new[k] = as<double>(rnorm(1, mu_n[k], post_prec)) ;
  }
  // simulate from prior if NAs
  LogicalVector isnan = is_nan(mu_new) ;
  if(!is_true(any(isnan)))
    return mu_new ;
  
  for(int k = 0; k < K; ++k){
    if(isnan[k])
      mu_new[k] = as<double>(rnorm(1, mu_0, sqrt(tau2_0))) ;
  }
  return mu_new ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_muchd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model_(xmod);
  Rcpp::S4 model = clone(model_);
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  double tau2_0 = hypp.slot("tau2.0") ;
  double tau2_0_tilde = 1/tau2_0 ;
  double mu_0 = hypp.slot("mu.0") ;
  
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector tau2_tilde = 1/tau2 ;
  IntegerVector z = model.slot("z") ;
  NumericMatrix theta = model.slot("theta") ;
  // below line not referenced anywhere in function
  // IntegerVector nn = model.slot("zfreq") ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  
  NumericVector tau2_B_tilde(K) ;;
  for(int k = 0; k < K; ++k) tau2_B_tilde[k] = tau2_0_tilde + B*tau2_tilde[k] ;
  
  NumericVector w1(K) ;
  NumericVector w2(K) ;
  for(int k = 0; k < K; ++k){
    w1[k] = tau2_0_tilde/(tau2_0_tilde + B*tau2_tilde[k]) ;
    w2[k] = B*tau2_tilde[k]/(tau2_0_tilde + B*tau2_tilde[k]) ;
  }
  NumericMatrix n_b = tableBatchZchd(model) ;
  NumericVector theta_bar(K) ;
  NumericVector th(K) ;
  
  for(int k = 0; k < K; ++k){
    double n_k = 0.0 ; // number of observations for component k
    double colsumtheta = 0.0;
    for(int i = 0; i < B; ++i){
      colsumtheta += n_b(i, k)*theta(i, k) ;
      n_k += n_b(i, k) ;
    }
    theta_bar[k] = colsumtheta/n_k ;
  }
  NumericVector mu_n(K) ;
  NumericVector mu_new(K) ;
  double post_prec ;
  for(int k=0; k<K; ++k){
    post_prec = sqrt(1.0/tau2_B_tilde[k]) ;
    mu_n[k] = w1[k]*mu_0 + w2[k]*theta_bar[k] ;
    mu_new[k] = as<double>(rnorm(1, mu_n[k], post_prec)) ;
  }
  // simulate from prior if NAs
  LogicalVector isnan = is_nan(mu_new) ;
  if(!is_true(any(isnan)))
    return mu_new ;
  
  for(int k = 0; k < K; ++k){
    if(isnan[k])
      mu_new[k] = as<double>(rnorm(1, mu_0, sqrt(tau2_0))) ;
  }
  return mu_new ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_tau2chd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod));
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  double m2_0 = hypp.slot("m2.0") ;
  int K = getK(hypp) ;
  double eta_0 = hypp.slot("eta.0") ;
  
  NumericVector mu = model.slot("mu_chd") ;
  NumericMatrix theta = model.slot("theta_chd") ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  double eta_B = eta_0 + B ;
  
  NumericVector s2_k(K) ;
  for(int k = 0; k < K; ++k){
    for(int i = 0; i < B; ++i){
      s2_k[k] += pow(theta(i, k) - mu[k], 2) ;
    }
  }
  NumericVector tau2(K) ;
  for(int k = 0; k < K; ++k) {
    double m2_k = 0.0 ;
    m2_k = 1.0/eta_B*(eta_0*m2_0 + s2_k[k]) ;
    tau2[k] = 1.0/as<double>(rgamma(1, 0.5*eta_B, 2.0/(eta_B*m2_k))) ;
  }
  return tau2 ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_sigma20chd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod));
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  NumericVector a = hypp.slot("a") ;
  NumericVector b = hypp.slot("b") ;
  NumericVector nu_0 = model.slot("nu.0chd") ;  
  NumericMatrix sigma2 = model.slot("sigma2_chd") ;
  Rcpp::NumericVector sigma2_0_old = model.slot("sigma2.0_chd");
  NumericVector prec(1) ;
  
  for(int i = 0; i < B; ++i){
    for(int k = 0; k < K; ++k){
      prec[0] += 1.0/sigma2(i, k) ;
    }
  }
  
  NumericVector a_k(1) ;
  NumericVector b_k(1) ;
  a_k[0] = a[0] + 0.5*(K * B)*nu_0[0] ;
  b_k[0] = b[0] + 0.5*nu_0[0]*prec[0] ;
  double rate ;
  rate = 1.0/b_k[0] ;
  //return b_k ;
  NumericVector sigma2_0(1) ;
  sigma2_0[0] = as<double>(rgamma(1, a_k[0], rate)) ;
  double constraint = model.slot(".internal.constraint");
  if (constraint > 0) {
    if (sigma2_0[0] < constraint) {
      return sigma2_0_old;
    }
    else {
      return sigma2_0;
    }
  }
  else {
    return sigma2_0;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector update_nu0chd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericMatrix sigma2 = model.slot("sigma2_chd") ;
  int B = sigma2.nrow() ;  
  double sigma2_0 = model.slot("sigma2.0_chd") ;  
  double prec = 0.0;
  double lprec = 0.0 ;
  double betas = hypp.slot("beta") ;
  for(int i = 0; i < B; ++i){
    for(int k = 0; k < K; ++k){
      prec += 1.0/sigma2(i, k) ;
      lprec += log(1.0/sigma2(i, k)) ;
    }
  }
  NumericVector x(100) ;  
  for(int i = 0; i < 100; i++)  x[i] = i+1 ;
  NumericVector lpnu0(100);
  NumericVector y1(100) ;
  NumericVector y2(100) ;
  NumericVector y3(100) ;
  NumericVector prob(100) ;
  y1 = (B*K)*(0.5*x*log(sigma2_0*0.5*x) - lgamma(x*0.5)) ;
  y2 = (0.5*x - 1.0) * lprec ;
  y3 = x*(betas + 0.5*sigma2_0*prec) ;
  lpnu0 =  y1 + y2 - y3 ;
  prob = exp(lpnu0) ; // - maxprob) ;
  prob = prob/sum(prob) ;
  NumericVector nu0(1) ;
  //int u ;
  NumericVector u(1) ;
  double cumprob = 0.0;
  // sample x with probability prob
  for(int i = 0; i < 100; i++){
    cumprob += prob[i] ;
    u = runif(1) ;
    if (u[0] < cumprob){
      nu0[0] = x[i] ;
      break ;
    }
  }
  return nu0 ;  
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_multinomialPrPar(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  NumericVector pp = model.slot("pi") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;

  double df = getDf(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  
  Rcpp::NumericVector xp = x[!child_ind];
  int M = xp.size() ;
  NumericMatrix lik(M, K) ;
  NumericVector this_batch(M) ;
  NumericVector tmp(M) ;
  NumericVector rowtotal(M) ;
  
  // more feedback
  for(int k = 0; k < K; ++k){
    NumericVector dens(M) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      double sigma = sqrt(sigma2(b, k));
      //tmp = p[k] * pp[k] * dlocScale_t(xp, df, theta(b, k), sigma) * this_batch ;
      //tmp = ((p[k]+pp[k])/2) * dlocScale_t(xp, df, theta(b, k), sigma) * this_batch ;
      tmp = pp[k] * dlocScale_t(xp, df, theta(b, k), sigma) * this_batch ;
      dens += tmp ;
    }
    lik(_, k) = dens ;
    rowtotal += dens ;
  }
  
  NumericMatrix PP(M, K) ;
  for(int k=0; k<K; ++k){
    PP(_, k) = lik(_, k)/rowtotal ;
  }
  return PP ;
}

// [[Rcpp::export]]
Rcpp::IntegerVector update_parents(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericMatrix theta = model.slot("theta") ;
  IntegerVector batch = model.slot("batch") ;
  int B = theta.nrow() ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zp = z[!child_ind];
  int parents_size = zp.size();
  NumericMatrix p(parents_size, K);
  p = update_multinomialPrPar(xmod) ;  // number trios x K
  
  //NumericMatrix cumP(n, K) ;
  //  Make more efficient
  //return cumP ;
  NumericVector upar = runif(parents_size) ;
  IntegerVector zpar_(parents_size) ;
  IntegerVector zpar = clone(zpar_) ;
  IntegerMatrix freq(B, K) ;
  
  int b ;
  for(int i=0; i < parents_size; i++){
    //initialize accumulator ;
    double acc = 0 ;
    for(int k = 0; k < K; k++){
      acc += p(i, k) ;
      if( upar[i] < acc ) {
        zpar[i] = k + 1 ;
        b = batch[i] - 1 ;
        freq(b, k) += 1 ;
        break ;
      }
    }
  }
  if(is_true(all(freq > 1))){
    return zpar ;
  } else{
    return zp ;
  }

}  

// [[Rcpp::export]]
Rcpp::IntegerVector update_zparents(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::IntegerVector ztrio = model.slot("z");
  
  //Rcpp::CharacterVector family_member=triodat["family_member"];
  // Rcpp::LogicalVector is_offspring;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  // return child_ind;
  
  IntegerVector zz_parents;
  zz_parents = update_parents(model); 
  int n = ztrio.size() ;
  int j = 0;
  for(int i = 0; i < n; i++){
    if(child_ind[i] == FALSE){
      ztrio[i] = zz_parents[j];
      j++;
    }
  }
  return ztrio;
}

// [[Rcpp::export]]
Rcpp::IntegerVector tableZpar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zp = z[!child_ind];
  
  Rcpp::IntegerVector nn(K) ;
  for(int k = 0; k < K; k++){
    nn[k] = sum(zp == (k+1)) ;
  }
  return nn ;
}

// [[Rcpp::export]]
Rcpp::IntegerVector tableZchd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zo = z[child_ind];
  
  Rcpp::IntegerVector nn(K) ;
  for(int k = 0; k < K; k++){
    nn[k] = sum(zo == (k+1)) ;
  }
  return nn ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_pp(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  // IntegerVector z = model.slot("z") ;  
  IntegerVector nn = model.slot("zfreq_parents");
  IntegerVector alpha = hypp.slot("alpha") ;
  NumericVector alpha_n(K) ;  // really an integer vector, but rdirichlet expects numeric
  for(int k=0; k < K; k++) alpha_n[k] = alpha[k] + nn[k] ;
  NumericVector pp(K) ;
  // pass by reference
  rdirichlet(alpha_n, pp) ;
  return pp ;
  //  return alpha_n ;
}

// [[Rcpp::export]]
Rcpp::NumericVector update_pc(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  // IntegerVector z = model.slot("z") ;  
  IntegerVector nn = model.slot("zfreq_chd");
  IntegerVector alpha = hypp.slot("alpha") ;
  NumericVector alpha_n(K) ;  // really an integer vector, but rdirichlet expects numeric
  for(int k=0; k < K; k++) alpha_n[k] = alpha[k] + nn[k] ;
  NumericVector pc(K) ;
  // pass by reference
  rdirichlet(alpha_n, pc) ;
  return pc ;
  //  return alpha_n ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_multinomialPrChild(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  NumericMatrix ptrio = update_trioPr(xmod) ;
  NumericVector pc = model.slot("pi_chd") ;
  NumericMatrix sigma2 = model.slot("sigma2_chd") ;
  NumericMatrix theta = model.slot("theta_chd") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;

  double df = getDf(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  
  Rcpp::NumericVector xo = x[child_ind];
  int M = xo.size() ;
  NumericMatrix lik(M, K) ;
  NumericVector this_batch(M) ;
  NumericVector tmp(M) ;
  NumericVector rowtotal(M) ;
  
  for(int k = 0; k < K; ++k){
    NumericVector dens(M) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      double sigma = sqrt(sigma2(b, k));
      tmp = pc[k] * ptrio(_,k) * dlocScale_t(xo, df, theta(b, k), sigma) * this_batch ;
      dens += tmp ;
    }
    lik(_, k) = dens;
    rowtotal += dens ;
  }
  NumericMatrix PC(M, K) ;
  for(int k=0; k<K; ++k){
    PC(_, k) = lik(_, k)/rowtotal ;
  } 
  return PC;
}

// [[Rcpp::export]]
Rcpp::IntegerVector update_offspring(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericMatrix theta = model.slot("theta_chd") ;
  IntegerVector batch = model.slot("batch") ;
  int B = theta.nrow() ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zo = z[child_ind];
  int child_size = zo.size();
  NumericMatrix p(child_size, K);
  //p = update_multinomialPrOff(xmod) ;  
  p = update_multinomialPrChild(xmod) ;
  
  //NumericMatrix cumP(n, K) ;
  //  Make more efficient
  //return cumP ;
  NumericVector uc = runif(child_size) ;
  IntegerVector zc_(child_size) ;
  IntegerVector zc = clone(zc_) ;
  IntegerMatrix freq(B, K) ;
  
  int b ;
  for(int i=0; i < child_size; i++){
    //initialize accumulator ;
    double acc = 0 ;
    for(int k = 0; k < K; k++){
      acc += p(i, k) ;
      if( uc[i] < acc ) {
        zc[i] = k + 1 ;
        b = batch[i] - 1 ;
        freq(b, k) += 1 ;
        break ;
      }
    }
  }
  if(is_true(all(freq > 1))){
    return zc ;
  } else {
    return zo;
  }
  
}  

// this selectively updates z for offspring

// [[Rcpp::export]]
Rcpp::IntegerVector update_zchild(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::IntegerVector ztrio = model.slot("z");
  
  //return ztrio ;
  //Rcpp::CharacterVector family_member=triodat["family_member"];
  // Rcpp::LogicalVector is_offspring;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  // return child_ind;
  
  IntegerVector zz_offspring;
  zz_offspring = update_offspring(model); // length 300
  int n = ztrio.size() ;
  int j = 0;
  for(int i = 0; i < n; i++){
    if(child_ind[i] == TRUE){
      ztrio[i] = zz_offspring[j];
      j++;
    }
  }
return ztrio;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_u_sums_batchpar(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zp = z[!child_ind];
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector u = model.slot("u") ;
  NumericVector upar = u[!child_ind];
  int n = upar.size() ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  NumericMatrix sums(B, K) ;
  for(int i = 0; i < n; i++){
    for(int b = 0; b < B; b++) {
      for(int k = 0; k < K; k++){
        if(zp[i] == k+1 & batch[i] == b+1){
          sums(b, k) += upar[i] ;
        }
      }
    }
  }
  //Rcpp::Rcout << "u sums:" << std::endl << sums << std::endl;
  return sums ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_heavy_sums_batchpar(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 xmod = clone(object) ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  Rcpp::NumericVector xp = x[!child_ind];
  int n = xp.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::IntegerVector zp = z[!child_ind];
  NumericVector u = model.slot("u") ;
  NumericVector upar = u[!child_ind];
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  // IntegerVector nn = model.slot("zfreq") ;
  NumericMatrix sums(B, K) ;
  xp = xp * upar ;
  for(int i = 0; i < n; i++){
    for(int b = 0; b < B; b++) {
      for(int k = 0; k < K; k++){
        if(zp[i] == k+1 & batch[i] == b+1){
          sums(b, k) += xp[i] ;
        }
      }
    }
  }
  //Rcpp::Rcout << "heavy sums:" << std::endl << sums << std::endl;
  return sums ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_thetapar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix n_hb = tableBatchZpar(model) ;
  NumericVector mu = model.slot("mu") ;
  int B = n_hb.nrow() ;
  NumericMatrix theta_new(B, K) ;
  double w1 = 0.0 ;
  double w2 = 0.0 ;
  double mu_n = 0.0 ;
  double tau_n = 0.0 ;
  double post_prec = 0.0 ;
  NumericVector u = model.slot("u") ;
  double df = getDf(hypp) ;
  // find heavy means by batch
  NumericMatrix data_sum =  compute_heavy_sums_batchpar(model) ;
  NumericMatrix sumu = compute_u_sums_batchpar(model) ;
  double heavyn = 0.0;
  double heavy_mean = 0.0;
  for (int b = 0; b < B; ++b) {
    for(int k = 0; k < K; ++k){
      heavyn = sumu(b, k) / df;
      post_prec = 1.0/tau2[k] + heavyn*1.0/sigma2(b, k) ;
      if (post_prec == R_PosInf) {
        throw std::runtime_error("Bad simulation. Run again with different start.");
      }
      w1 = (1.0/tau2[k])/post_prec ;
      w2 = (heavyn * 1.0/sigma2(b, k))/post_prec ;
      heavy_mean = data_sum(b, k) / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean ;
      tau_n = sqrt(1.0/post_prec) ;
      theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
    }
  }
  return theta_new ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_u_sums_batchchd(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(xmod) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zo = z[child_ind];
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector u = model.slot("u") ;
  NumericVector uo = u[child_ind];
  int n = uo.size() ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  NumericMatrix sums(B, K) ;
  for(int i = 0; i < n; i++){
    for(int b = 0; b < B; b++) {
      for(int k = 0; k < K; k++){
        if(zo[i] == k+1 & batch[i] == b+1){
          sums(b, k) += uo[i] ;
        }
      }
    }
  }
  //Rcpp::Rcout << "u sums:" << std::endl << sums << std::endl;
  return sums ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_heavy_sums_batchchd(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 xmod = clone(object) ;
  Rcpp::S4 model(xmod) ;
  NumericVector x = model.slot("data") ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  Rcpp::NumericVector xo = x[child_ind];
  int n = xo.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::IntegerVector zo = z[child_ind];
  NumericVector u = model.slot("u") ;
  NumericVector uo = u[child_ind];
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  // IntegerVector nn = model.slot("zfreq") ;
  NumericMatrix sums(B, K) ;
  xo = xo * uo ;
  for(int i = 0; i < n; i++){
    for(int b = 0; b < B; b++) {
      for(int k = 0; k < K; k++){
        if(zo[i] == k+1 & batch[i] == b+1){
          sums(b, k) += xo[i] ;
        }
      }
    }
  }
  //Rcpp::Rcout << "heavy sums:" << std::endl << sums << std::endl;
  return sums ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_thetachd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector tau2 = model.slot("tau2_chd") ;
  NumericMatrix sigma2 = model.slot("sigma2_chd") ;
  NumericMatrix n_hb = tableBatchZchd(model) ;
  NumericVector mu = model.slot("mu_chd") ;
  int B = n_hb.nrow() ;
  NumericMatrix theta_new(B, K) ;
  double w1 = 0.0 ;
  double w2 = 0.0 ;
  double mu_n = 0.0 ;
  double tau_n = 0.0 ;
  double post_prec = 0.0 ;
  NumericVector u = model.slot("u") ;
  double df = getDf(hypp) ;
  // find heavy means by batch
  NumericMatrix data_sum =  compute_heavy_sums_batchchd(model) ;
  NumericMatrix sumu = compute_u_sums_batchchd(model) ;
  double heavyn = 0.0;
  double heavy_mean = 0.0;
  for (int b = 0; b < B; ++b) {
    for(int k = 0; k < K; ++k){
      heavyn = sumu(b, k) / df;
      post_prec = 1.0/tau2[k] + heavyn*1.0/sigma2(b, k) ;
      if (post_prec == R_PosInf) {
        throw std::runtime_error("Bad simulation. Run again with different start.");
      }
      w1 = (1.0/tau2[k])/post_prec ;
      w2 = (heavyn * 1.0/sigma2(b, k))/post_prec ;
      heavy_mean = data_sum(b, k) / heavyn / df;
      mu_n = w1*mu[k] + w2*heavy_mean ;
      tau_n = sqrt(1.0/post_prec) ;
      theta_new(b, k) = as<double>(rnorm(1, mu_n, tau_n)) ;
    }
  }
  return theta_new ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_sigma2par(Rcpp::S4 xmod){
  Rcpp::RNGScope scope;
  Rcpp::S4 model(clone(xmod)) ;
  // get model
  // get parameter estimates
  Rcpp::NumericMatrix theta = model.slot("theta");
  double nu_0 = model.slot("nu.0");
  double sigma2_0 = model.slot("sigma2.0");
  
  // get data and size attributes
  Rcpp::NumericVector x = model.slot("data");
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  Rcpp::NumericVector xp = x[!child_ind];
  Rcpp::IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zp = z[!child_ind];
  int n = xp.size();
  int K = theta.ncol();
  int B = theta.nrow();
  NumericVector u = model.slot("u") ;
  double df = getDf(model.slot("hyperparams")) ;
  
  //IntegerVector nn = model.slot("zfreq");
  // get batch info
  Rcpp::NumericMatrix tabz = tableBatchZchd(model);
  Rcpp::IntegerVector batch = model.slot("batch");
  Rcpp::IntegerVector ub = unique_batch(batch);
  Rcpp::NumericMatrix ss(B, K);
  
  for (int i = 0; i < n; ++i) {
    for (int b = 0; b < B; ++b) {
      if (batch[i] != ub[b]) {
        continue;
      }
      
      for (int k = 0; k < K; ++k){
        if(zp[i] == k+1 & batch[i] == b+1){
          ss(b, k) += u[i] * pow(x[i] - theta(b, k), 2);
          //ss(b, k) += pow(x[i] - theta(b, k), 2);
        }
      }
    }
  }
  
  //NumericMatrix sigma2_nh(B, K);
  double shape;
  double rate;
  double sigma2_nh;
  double nu_n;
  Rcpp::NumericMatrix sigma2_tilde(B, K);
  Rcpp::NumericMatrix sigma2_(B, K);
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) {
      nu_n = nu_0 + tabz(b, k);
      sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k)/df);
      // sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k));
      shape = 0.5 * nu_n;
      rate = shape * sigma2_nh;
      sigma2_tilde(b, k) = Rcpp::as<double>(rgamma(1, shape, 1.0/rate));
      sigma2_(b, k) = 1.0 / sigma2_tilde(b, k);
    }
  }
  
  return sigma2_;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_sigma2chd(Rcpp::S4 xmod){
  Rcpp::RNGScope scope;
  Rcpp::S4 model(clone(xmod)) ;
  // get model
  // get parameter estimates
  Rcpp::NumericMatrix theta = model.slot("theta_chd");
  double nu_0 = model.slot("nu.0chd");
  double sigma2_0 = model.slot("sigma2.0_chd");
  
  // get data and size attributes
  Rcpp::NumericVector x = model.slot("data");
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  Rcpp::NumericVector xc = x[child_ind];
  Rcpp::IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zo = z[child_ind];
  int n = xc.size();
  int K = theta.ncol();
  int B = theta.nrow();
  NumericVector u = model.slot("u") ;
  double df = getDf(model.slot("hyperparams")) ;
  
  //IntegerVector nn = model.slot("zfreq");
  // get batch info
  Rcpp::NumericMatrix tabz = tableBatchZchd(model);
  Rcpp::IntegerVector batch = model.slot("batch");
  Rcpp::IntegerVector ub = unique_batch(batch);
  Rcpp::NumericMatrix ss(B, K);
  
  for (int i = 0; i < n; ++i) {
    for (int b = 0; b < B; ++b) {
      if (batch[i] != ub[b]) {
        continue;
      }
      
      for (int k = 0; k < K; ++k){
        if(zo[i] == k+1 & batch[i] == b+1){
          ss(b, k) += u[i] * pow(x[i] - theta(b, k), 2);
          //ss(b, k) += pow(x[i] - theta(b, k), 2);
        }
      }
    }
  }
  
  //NumericMatrix sigma2_nh(B, K);
  double shape;
  double rate;
  double sigma2_nh;
  double nu_n;
  Rcpp::NumericMatrix sigma2_tilde(B, K);
  Rcpp::NumericMatrix sigma2_(B, K);
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) {
      nu_n = nu_0 + tabz(b, k);
      sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k)/df);
      // sigma2_nh = 1.0/nu_n*(nu_0*sigma2_0 + ss(b, k));
      shape = 0.5 * nu_n;
      rate = shape * sigma2_nh;
      sigma2_tilde(b, k) = Rcpp::as<double>(rgamma(1, shape, 1.0/rate));
      sigma2_(b, k) = 1.0 / sigma2_tilde(b, k);
    }
  }
  
  return sigma2_;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_vars2(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  NumericVector x = model.slot("data") ;
  int n = x.size() ;
  IntegerVector z = model.slot("z") ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  // line below not referenced anywhere here
  // IntegerVector nn = model.slot("zfreq") ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  int B = ub.size() ;
  NumericMatrix vars(B, K) ;
  NumericMatrix tabz = tableBatchZpar(model) ;
  NumericMatrix mn = model.slot("data.mean") ;
  NumericVector this_batch(n) ;
  NumericVector is_z(n) ;
  NumericVector tau2 = model.slot("tau2") ;
  NumericVector total(1) ;
  for(int b = 0; b < B; ++b){
    this_batch = batch == ub[b] ;
    for(int k = 0; k < K; ++k){
      is_z = z == (k + 1) ;
      total[0] = sum(is_z * this_batch) ;
      if(total[0] <= 1){
        vars(b, k) = tau2[k] ;
      } else {
        vars(b, k) = sum(pow(x - mn(b,k), 2.0) * this_batch * is_z) / (tabz(b, k) - 1) ;
      }
    }
  }
  return vars ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compute_prec2(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  NumericMatrix vars = compute_vars2(model) ;
  int B = vars.nrow() ;
  int K = vars.ncol() ;
  NumericMatrix prec(B, K) ;
  for(int b = 0; b < B; ++b){
    for(int k = 0; k < K; ++k){
      prec(b, k) = 1.0/vars(b, k) ;
    }
  }
  return prec ;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix update_probzpar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = hypp.slot("k") ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zp = z[!child_ind];
  NumericMatrix theta = model.slot("theta") ;
  int N = zp.size() ;
  IntegerMatrix pZ = model.slot("probz_par") ;
  //
  // update probz such that the z value corresponding to the lowest
  // mean is 1, the second lowest mean is 2, etc.
  //
  // Assume that all batches have the same ordering and so here we
  // just use the first batch
  //
  NumericVector means(K) ;
  means = theta(0, _) ;
  NumericVector cn(K) ;
  cn = order_(means) ;
  for(int k = 0; k < K; ++k) cn[k] = cn[k] - 1 ;
  for(int i = 0; i < N; ++i){
    for(int k = 0; k < K; ++k){
      if(zp[i] == (k + 1)){
        pZ(i, cn[k]) += 1;
      }
    }
  }
  return pZ ;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix update_probzchd(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = hypp.slot("k") ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zo = z[child_ind];
  NumericMatrix theta = model.slot("theta_chd") ;
  int N = zo.size() ;
  IntegerMatrix pZ = model.slot("probz_chd") ;
  //
  // update probz such that the z value corresponding to the lowest
  // mean is 1, the second lowest mean is 2, etc.
  //
  // Assume that all batches have the same ordering and so here we
  // just use the first batch
  //
  NumericVector means(K) ;
  means = theta(0, _) ;
  NumericVector cn(K) ;
  cn = order_(means) ;
  for(int k = 0; k < K; ++k) cn[k] = cn[k] - 1 ;
  for(int i = 0; i < N; ++i){
    for(int k = 0; k < K; ++k){
      if(zo[i] == (k + 1)){
        pZ(i, cn[k]) += 1;
      }
    }
  }
  return pZ ;
}

Rcpp::IntegerVector update_z2(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  NumericVector x = model.slot("data") ;
  NumericMatrix theta = model.slot("theta") ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector z = model.slot("z");
  int B = theta.nrow() ;
  int n = x.size() ;
  NumericMatrix p(n, K);
  p = update_multinomialPr(xmod) ;
  //NumericMatrix cumP(n, K) ;
  //  Make more efficient
  //return cumP ;
  NumericVector u = runif(n) ;
  IntegerVector zz_(n) ;
  IntegerVector zz = clone(zz_) ;
  IntegerMatrix freq(B, K) ;
  int b ;
  for(int i=0; i < n; i++){
    //initialize accumulator ;
    double acc = 0 ;
    for(int k = 0; k < K; k++){
      acc += p(i, k) ;
      if( u[i] < acc ) {
        zz[i] = k + 1 ;
        b = batch[i] - 1 ;
        freq(b, k) += 1 ;
        break ;
      }
    }
  }
  if(is_true(all(freq > 1))){
    return zz ;
  } else {
    return z ;
  }

}


// [[Rcpp::export]]
Rcpp::S4 trios_burnin(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int S = params.slot("burnin") ;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;
  double df = getDf(model.slot("hyperparams")) ;
  if( S < 1 ){
    return model ;
  }
  for(int s = 0; s < S; ++s){
    if(up[7] > 0){
      //model.slot("z") = update_z2(model) ;
      model.slot("z") = update_zparents(model) ;
      model.slot("z") = update_zchild(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      model.slot("zfreq_parents") = tableZpar(model) ;
      model.slot("zfreq_chd") = tableZchd(model) ;
    }
    if(up[0] > 0){
        model.slot("theta") = update_thetapar(model) ;
      model.slot("theta_chd") = update_thetachd(model) ;
      }
      if(up[1] > 0)
        model.slot("sigma2") = update_sigma2par(model) ;
      model.slot("sigma2_chd") = update_sigma2chd(model) ;
      if(up[3] > 0)
        model.slot("mu") = update_mupar(model) ;
      model.slot("mu_chd") = update_muchd(model) ;
      if(up[4] > 0)
        model.slot("tau2") = update_tau2(model) ;
      model.slot("tau2_chd") = update_tau2chd(model) ;
      if(up[6] > 0)
        model.slot("sigma2.0") = update_sigma20(model) ;
      model.slot("sigma2.0_chd") = update_sigma20chd(model) ;
      if(up[5] > 0)
        model.slot("nu.0") = update_nu0(model) ;
        model.slot("nu.0chd") = update_nu0chd(model) ;
      if(up[2] > 0)
        model.slot("pi") = update_pp(model) ;
      model.slot("pi_chd") = update_pc(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
  }
  NumericVector lls2(1);
  NumericVector ll(1);
  lls2 = stageTwoLogLikBatch(model);
  ll = compute_loglik(model);
  ll = ll + lls2;
  model.slot("loglik") = ll;
  model.slot("logprior") = compute_logprior(model) ;
  return model ;
  // return vars ;
}

// [[Rcpp::export]]
Rcpp::S4 trios_mcmc(Rcpp::S4 object, Rcpp::S4 mcmcp) {
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  Rcpp::S4 chain(model.slot("mcmc.chains")) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  Rcpp::S4 params(mcmcp) ;
  IntegerVector up = params.slot("param_updates") ;
  int K = getK(hypp) ;
  int T = params.slot("thin") ;
  int S = params.slot("iter") ;
  if( S < 1 ) return model ;
  NumericVector x = model.slot("data") ;
  int N = x.size() ;
  double df = getDf(model.slot("hyperparams")) ;
  NumericMatrix thetac = chain.slot("theta") ;
  NumericMatrix sigma2c = chain.slot("sigma2") ;
  NumericMatrix thetachdc = chain.slot("theta_chd") ;
  NumericMatrix sigma2chdc = chain.slot("sigma2_chd") ;
  NumericMatrix th = model.slot("theta");
  NumericMatrix thchd = model.slot("theta_chd");
  int B = th.nrow();
  NumericMatrix s2(B, K);
  NumericMatrix pmix = chain.slot("pi") ;
  NumericMatrix pmix_chd = chain.slot("pi_chd") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  NumericMatrix zfreq_parents = chain.slot("zfreq_parents") ;
  NumericMatrix zfreq_chd = chain.slot("zfreq_chd") ;
  NumericMatrix mu = chain.slot("mu") ;
  NumericMatrix mu_chd = chain.slot("mu_chd") ;
  NumericMatrix tau2 = chain.slot("tau2") ;
  NumericMatrix tau2_chd = chain.slot("tau2_chd") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector nu0_chd = chain.slot("nu.0chd") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector sigma2_0chd = chain.slot("sigma2.0_chd") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  //NumericVector th(K) ;
  //NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector pc(K) ;
  NumericVector m(K) ; //mu
  NumericVector mc(K) ; //mu
  NumericVector t2(K) ;//tau2
  NumericVector t2c(K) ;//tau2
  NumericVector n0(1) ;//nu0
  NumericVector n0c(1) ;//nu0
  IntegerVector z(N) ;
  NumericVector u(N) ;
  NumericVector s20(1) ; //sigma2_0
  NumericVector s20c(1) ; //sigma2_0
  NumericVector lls2(1) ;  // stage 2 log lik
  NumericVector ll(1) ;
  NumericVector lp(1) ;
  IntegerVector tmp(K) ;
  IntegerVector zf(K) ;
  IntegerVector zp(K) ;
  IntegerVector zc(K) ;
  // Initial values
  p = model.slot("pi") ;
  pc = model.slot("pi_chd") ;
  m = model.slot("mu") ;
  mc = model.slot("mu_chd") ;
  t2 = model.slot("tau2") ;
  t2c = model.slot("tau2_chd") ;
  n0 = model.slot("nu.0") ;
  n0c = model.slot("nu.0chd") ;
  s20 = model.slot("sigma2.0") ;
  s20c = model.slot("sigma2.0_chd") ;
  zf = model.slot("zfreq") ;
  zp = model.slot("zfreq_parents") ;
  zc = model.slot("zfreq_chd") ;
  z = model.slot("z") ;
  u = model.slot("u") ;
  ll = model.slot("loglik") ;
  lp = model.slot("logprior") ;
  // Record initial values in chains
  mu(0, _) = m ;
  mu_chd(0, _) = mc ;
  nu0[0] = n0[0] ;
  nu0_chd[0] = n0c[0] ;
  sigma2_0[0] = s20[0] ;
  sigma2_0chd[0] = s20c[0] ;
  loglik_[0] = ll[0] ;
  logprior_[0] = lp[0] ;
  thetac(0, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
  thetachdc(0, _) = as<Rcpp::NumericVector>(model.slot("theta_chd")) ;
  tau2(0, _) = t2 ;
  tau2_chd(0, _) = t2c ;
  sigma2c(0, _) = as<Rcpp::NumericVector>(model.slot("sigma2")) ;
  sigma2chdc(0, _) = as<Rcpp::NumericVector>(model.slot("sigma2_chd")) ;
  pmix(0, _) = p ;
  pmix_chd (0, _) = pc ;
  zfreq(0, _) = zf ;
  zfreq_parents(0, _) = zp ;
  zfreq_chd(0, _) = zc ;
  for(int s = 1; s < S; ++s){
    if(up[7] > 0){
      z = update_zparents(model) ;
      model.slot("z") = z ;
      tmp = tableZpar(model) ;
      model.slot("zfreq_parents") = tmp ;
      //model.slot("zfreq") = tmp ;
      //model.slot("probz") = update_probz(model) ;
    } else {
      z = model.slot("z") ;
      tmp = model.slot("zfreq") ;
      //tmp = model.slot("zfreq") ;
    }
    //zfreq_parents(s, _) = tmp ;
    zfreq_parents(s, _) = tmp ;
    
    if(up[8] > 0){
      z = update_zchild(model) ;
      model.slot("z") = z ;
      tmp = tableZchd(model) ;
      model.slot("zfreq_chd") = tmp ;
      model.slot("probz") = update_probz(model) ;
      model.slot("probz_par") = update_probzpar(model) ;
      model.slot("probz_chd") = update_probzchd(model) ;
    } else {
      z = model.slot("z") ;
      tmp = model.slot("zfreq_chd") ;
    }
    //zfreq(s, _) = tmp ;
    zfreq_chd(s, _) = tmp ;
    
    tmp = tableZ(K, model.slot("z")) ; ;
    model.slot("zfreq") = tmp ;
    zfreq(s, _) = tmp ;
    
    if(up[9] > 0){
      p = update_pp(model) ;
      model.slot("pi") = p ;
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[2] > 0){
      pc = update_pc(model) ;
      model.slot("pi_chd") = pc ;
    } else {
      pc = model.slot("pi_chd") ;
    }
    pmix_chd(s, _) = pc ;
    if(up[0] > 0) {
      model.slot("theta") = update_thetapar(model) ;
      model.slot("theta_chd") = update_thetachd(model) ;
    }
    thetac(s, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
    thetachdc(s, _) = as<Rcpp::NumericVector>(model.slot("theta_chd")) ;
    if(up[1] > 0){
      model.slot("sigma2") = update_sigma2par(model) ;
      model.slot("sigma2_chd") = update_sigma2chd(model) ;
    }
    sigma2c(s, _) = as<Rcpp::NumericVector>(model.slot("sigma2"));
    sigma2chdc(s, _) = as<Rcpp::NumericVector>(model.slot("sigma2_chd"));
    if(up[3] > 0){
      m = update_mupar(model) ;
      model.slot("mu") = m ;
      mc = update_muchd(model) ;
      model.slot("mu_chd") = mc ;
      
    } else {
      m = model.slot("mu") ;
      mc = model.slot("mu_chd") ;
    }
    mu(s, _) = m ;
    mu_chd(s, _) = mc ;
    if(up[4] > 0){
      t2 = update_tau2(model) ;
      model.slot("tau2") = t2 ;
      t2c = update_tau2chd(model) ;
      model.slot("tau2_chd") = t2c ;
    } else {
      t2 = model.slot("tau2") ;
      t2c = model.slot("tau2_chd") ;
    }
    tau2(s, _) = t2 ;
    tau2_chd(s, _) = t2c ;
    if(up[5] > 0){
      n0 = update_nu0(model) ;
      model.slot("nu.0") = n0 ;
      n0c = update_nu0chd(model) ;
      model.slot("nu.0chd") = n0c ;
    } else {
      n0 = model.slot("nu.0") ;
      n0c = model.slot("nu.0chd") ;
    }
    nu0[s] = n0[0] ;
    nu0_chd[s] = n0c[0] ;
    if(up[6] > 0){
      s20 = update_sigma20(model) ;
      model.slot("sigma2.0") = s20 ;
      s20c = update_sigma20chd(model) ;
      model.slot("sigma2.0_chd") = s20c ;
    } else {
      s20 = model.slot("sigma2.0") ;
      s20c = model.slot("sigma2.0_chd") ;
    }
    sigma2_0[s] = s20[0] ;
    sigma2_0chd[s] = s20c[0] ;
    ll = compute_loglik(model) ;
    lls2 = stageTwoLogLikBatch(model) ;
    ll = ll + lls2 ;
    loglik_[s] = ll[0] ;
    model.slot("loglik") = ll ;
    lp = compute_logprior(model) ;
    logprior_[s] = lp[0] ;
    model.slot("logprior") = lp ;
    u = Rcpp::rchisq(N, df) ;
    model.slot("u") = u;
    // Thinning
    for(int t = 0; t < T; ++t){
      if(up[7] > 0){
        model.slot("z") = update_zparents(model) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
        model.slot("zfreq_parents") = tableZpar(model) ;
        //tmp = model.slot("zfreq_parents") ;
        //model.slot("zfreq") = tmp ;
      }
      if(up[8] > 0){
        model.slot("z") = update_zchild(model) ;
        model.slot("zfreq_chd") = tableZchd(model) ;
      }
      if(up[9] > 0)
        model.slot("pi") = update_pp(model) ;
      if(up[2] > 0)
        model.slot("pi_chd") = update_pc(model) ;
      if(up[0] > 0)
        model.slot("theta") = update_thetapar(model) ;
      model.slot("theta_chd") = update_thetachd(model) ;
      if(up[1] > 0)
        model.slot("sigma2") = update_sigma2par(model) ;
      model.slot("sigma2_chd") = update_sigma2chd(model) ;
      if(up[3] > 0)
        model.slot("mu") = update_mupar(model) ;
      model.slot("mu_chd") = update_muchd(model) ;
      if(up[4] > 0)
        model.slot("tau2") = update_tau2(model) ;
      model.slot("tau2_chd") = update_tau2chd(model) ;
      if(up[5] > 0)
        model.slot("nu.0") = update_nu0(model) ;
      model.slot("nu.0chd") = update_nu0chd(model) ;
      if(up[6] > 0)
        model.slot("sigma2.0") = update_sigma20(model) ;
      model.slot("sigma2.0_chd") = update_sigma20chd(model) ;
      model.slot("u") = Rcpp::rchisq(N, df) ;
    }
  }
  //
  // assign chains back to object
  //
  chain.slot("theta") = thetac ;
  chain.slot("theta_chd") = thetachdc ;
  chain.slot("sigma2") = sigma2c ;
  chain.slot("sigma2_chd") = sigma2chdc ;
  chain.slot("pi") = pmix ;
  chain.slot("pi_chd") = pmix_chd ;
  chain.slot("mu") = mu ;
  chain.slot("mu_chd") = mu_chd ;
  chain.slot("tau2") = tau2 ;
  chain.slot("tau2_chd") = tau2_chd ;
  chain.slot("nu.0") = nu0 ;
  chain.slot("nu.0chd") = nu0_chd;
  chain.slot("sigma2.0") = sigma2_0 ;
  chain.slot("sigma2.0_chd") = sigma2_0chd ;
  chain.slot("zfreq") = zfreq ;
  chain.slot("zfreq_parents") = zfreq_parents ;
  chain.slot("zfreq_chd") = zfreq_chd ;
  chain.slot("loglik") = loglik_ ;
  chain.slot("logprior") = logprior_ ;
  model.slot("mcmc.chains") = chain ;
  return model ;
}

// [[Rcpp::export]]
Rcpp::S4  z2cn(Rcpp::S4 xmod, Rcpp::IntegerVector map){
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::IntegerVector ztrio = model.slot("z");
  //Rcpp::IntegerVector map = model.slot("maplabel");
  //Rcpp::IntegerVector map2 = map.sort();
  // map2.erase(std::unique(map2.begin(), map2.end()), map2.end());
  //int map2_max = max(map);
  for (int i = 0; i < ztrio.size(); i++){
    int zind = ztrio[i] - 1;
    int maplab = map[zind];
    ztrio[i] = maplab;
  }
  model.slot("z") = ztrio ;
  return model;
}
