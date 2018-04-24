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
  Rcpp::LogicalVector index(ind.size());
  index = ind==TRUE & ind2==TRUE;
  int nr=f.size();
  int j;
  for(int i = 0; i < nr; i ++){
    if(index[i]){
      j = i ;
      break;
    }
  }
  Rcpp::NumericVector result(mprob.ncol());
  result=mprob(j, _);
  return result;
}

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
Rcpp::IntegerVector update_offspring(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  Rcpp::DataFrame triodat(model.slot("triodata"));
  NumericMatrix theta = model.slot("theta") ;
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
  p = update_trioPr(xmod) ;  // number trios x K
  
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
  }
  
  // avoid the control may reach of non-void function
  // will prob need to fix later to return subset of slot z
 // int counter = model.slot(".internal.counter");
// counter++;
//  model.slot(".internal.counter") = counter;
//  return model.slot("z") ;
}  

// [[Rcpp::export]]
Rcpp::IntegerVector z2cn(Rcpp::IntegerVector ztrio, Rcpp::IntegerVector map){
  //Rcpp::S4 model(clone(xmod)) ;
  //Rcpp::IntegerVector map = model.slot("maplabel");
  //Rcpp::IntegerVector map2 = map.sort();
  // map2.erase(std::unique(map2.begin(), map2.end()), map2.end());
  //int map2_max = max(map);
  
  //Rcpp::IntegerVector ztemp = model.slot("z");

  for (int i = 0; i < ztrio.size(); i++){
    int zind = ztrio[i] - 1;
    int maplab = map[zind];
      ztrio[i] = maplab;
    }
  
  return ztrio;
}
  
// need to map components to copy number
// [[Rcpp::export]]
Rcpp::IntegerVector update_ztrio(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::IntegerVector map = model.slot("maplabel");
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  Rcpp::IntegerVector sts = getSt(hypp);
  Rcpp::DataFrame triodat(model.slot("triodata"));
  int n = triodat.size() ;
  NumericVector u = runif(n) ;
  IntegerVector ztrio_(n) ;
  IntegerVector ztrio = clone(ztrio_) ;
  ztrio = update_z(model) ;
  
  //Rcpp::CharacterVector family_member=triodat["family_member"];
  // Rcpp::LogicalVector is_offspring;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector zz_offspring;
  zz_offspring = update_offspring(model); // length 300
  
  n = zz_offspring.size();
  int j = 0;
  for(int i = 0; i < n; i++){
  if(child_ind[i] == TRUE){
  ztrio[i] = zz_offspring[j];
  j++;
    }
  }
  
   return ztrio;
  
  Rcpp::IntegerVector ztrio2;
  ztrio2 = z2cn(ztrio, map);
  //return ztrio2;
  
  //
  // Don't update ztrio if there are states with zero frequency.
  //
  int counter = model.slot(".internal.counter");
  counter++;
  model.slot(".internal.counter") = counter;
  return model.slot("z") ;
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
      model.slot("z") = update_ztrio(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    }
    if(up[0] > 0)
      model.slot("theta") = update_theta(model) ;
    if(up[1] > 0)
      model.slot("sigma2") = update_sigma2(model) ;
    if(up[3] > 0)
      model.slot("mu") = update_mu(model) ;
    if(up[4] > 0)
      model.slot("tau2") = update_tau2(model) ;
    if(up[6] > 0)
      model.slot("sigma2.0") = update_sigma20(model) ;
    if(up[5] > 0)
      model.slot("nu.0") = update_nu0(model) ;
    if(up[2] > 0)
      model.slot("pi") = update_p(model) ;
    model.slot("u") = Rcpp::rchisq(N, df) ;
  }
  // compute log prior probability from last iteration of burnin
  // compute log likelihood from last iteration of burnin
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
  NumericMatrix th = model.slot("theta");
  int B = th.nrow();
  NumericMatrix s2(B, K);
  NumericMatrix pmix = chain.slot("pi") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  NumericMatrix mu = chain.slot("mu") ;
  NumericMatrix tau2 = chain.slot("tau2") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  //NumericVector th(K) ;
  //NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector m(K) ; //mu
  NumericVector t2(K) ;//tau2
  NumericVector n0(1) ;//nu0
  IntegerVector z(N) ;
  NumericVector u(N) ;
  NumericVector s20(1) ; //sigma2_0
  NumericVector lls2(1) ;  // stage 2 log lik
  NumericVector ll(1) ;
  NumericVector lp(1) ;
  IntegerVector tmp(K) ;
  IntegerVector zf(K) ;
  // Initial values
  p = model.slot("pi") ;
  m = model.slot("mu") ;
  t2 = model.slot("tau2") ;
  n0 = model.slot("nu.0") ;
  s20 = model.slot("sigma2.0") ;
  zf = model.slot("zfreq") ;
  z = model.slot("z") ;
  u = model.slot("u") ;
  ll = model.slot("loglik") ;
  lp = model.slot("logprior") ;
  // Record initial values in chains
  mu(0, _) = m ;
  nu0[0] = n0[0] ;
  sigma2_0[0] = s20[0] ;
  loglik_[0] = ll[0] ;
  logprior_[0] = lp[0] ;
  thetac(0, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
  tau2(0, _) = t2 ;
  sigma2c(0, _) = as<Rcpp::NumericVector>(model.slot("sigma2")) ;
  pmix(0, _) = p ;
  zfreq(0, _) = zf ;

  // Is accessing a slot in an object expensive?
  // Currently, there is no alternative as the current values are
  // stored in the object.  Hence, the entire object has to be passed
  // to the updating functions.
  // start at 1 instead of zero. Initial values are as above
  //up[7] = 0;
  //up[6] = 0;
  //up[5] = 0;
  //up[4] = 0;
  //up[3] = 0;
  //up[2] = 0;
  //up[1] = 0;
  //up[0] = 0;
  for(int s = 1; s < S; ++s){
    if(up[7] > 0){
      z = update_ztrio(model) ;
      model.slot("z") = z ;
      tmp = tableZ(K, z) ;
      model.slot("probz") = update_probz(model) ;
      model.slot("zfreq") = tmp ;
    } else {
      tmp = model.slot("zfreq") ;
    }
    zfreq(s, _) = tmp ;
    if(up[0] > 0) {
      model.slot("theta") = update_theta(model) ;
    }
    thetac(s, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
    if(up[1] > 0){
      model.slot("sigma2") = update_sigma2(model) ;
    }
    sigma2c(s, _) = as<Rcpp::NumericVector>(model.slot("sigma2"));
    if(up[2] > 0){
      p = update_p(model) ;
      model.slot("pi") = p ;
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[3] > 0){
      m = update_mu(model) ;
      model.slot("mu") = m ;
    } else {
      m = model.slot("mu") ;
    }
    mu(s, _) = m ;
    if(up[4] > 0){
      t2 = update_tau2(model) ;
      model.slot("tau2") = t2 ;
    } else {
      t2 = model.slot("tau2") ;
    }
    tau2(s, _) = t2 ;
    if(up[5] > 0){
      n0 = update_nu0(model) ;
      model.slot("nu.0") = n0 ;
    } else {
      n0 = model.slot("nu.0") ;
    }
    nu0[s] = n0[0] ;
    if(up[6] > 0){
      s20 = update_sigma20(model) ;
      model.slot("sigma2.0") = s20 ;
    } else {
      s20 = model.slot("sigma2.0") ;
    }
    sigma2_0[s] = s20[0] ;
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
        model.slot("z") = update_ztrio(model) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      }
      if(up[0] > 0)
        model.slot("theta") = update_theta(model) ;
      if(up[1] > 0)
        model.slot("sigma2") = update_sigma2(model) ;
      if(up[2] > 0)
        model.slot("pi") = update_p(model) ;
      if(up[3] > 0)
        model.slot("mu") = update_mu(model) ;
      if(up[4] > 0)
        model.slot("tau2") = update_tau2(model) ;
      if(up[5] > 0)
        model.slot("nu.0") = update_nu0(model) ;
     if(up[6] > 0)
       model.slot("sigma2.0") = update_sigma20(model) ;
     model.slot("u") = Rcpp::rchisq(N, df) ;
    }
  }
  //
  // assign chains back to object
  //
  chain.slot("theta") = thetac ;
  chain.slot("sigma2") = sigma2c ;
  chain.slot("pi") = pmix ;
  chain.slot("mu") = mu ;
  chain.slot("tau2") = tau2 ;
  chain.slot("nu.0") = nu0 ;
  chain.slot("sigma2.0") = sigma2_0 ;
  chain.slot("zfreq") = zfreq ;
  chain.slot("loglik") = loglik_ ;
  chain.slot("logprior") = logprior_ ;
  model.slot("mcmc.chains") = chain ;
  return model ;
}

