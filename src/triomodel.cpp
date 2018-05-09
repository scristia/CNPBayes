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
Rcpp::NumericMatrix update_multinomialPrPar(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  NumericVector p = model.slot("pi") ;
  NumericVector pp = model.slot("pi_parents") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;
  IntegerVector nb = model.slot("batchElements") ;
  int N = x.size() ;
  NumericMatrix lik(N, K) ;
  NumericVector this_batch(N) ;
  NumericVector tmp(N) ;
  NumericVector rowtotal(N) ;
  double df = getDf(hypp) ;
  
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  
  Rcpp::NumericVector xp = x[!child_ind];
  int M = xp.size() ;
  
  // added an additional p[k] weight for experimentation/
  // more feedback
  for(int k = 0; k < K; ++k){
    NumericVector dens(N) ;
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
  return zpar ;
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

// this is update for pi_parents only

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
Rcpp::NumericMatrix update_multinomialPrChild(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  NumericMatrix ptrio = update_trioPr(xmod) ;
  NumericVector p = model.slot("pi") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;
  IntegerVector nb = model.slot("batchElements") ;
  int N = ptrio.size() ;
  NumericVector this_batch(N) ;
  NumericVector tmp ;
  NumericVector rowtotal(N) ;
  double df = getDf(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  
  Rcpp::NumericVector xo = x[child_ind];
  int M = xo.size() ;
  NumericMatrix lik(M, K) ;
  
  for(int k = 0; k < K; ++k){
    NumericVector dens(N) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      double sigma = sqrt(sigma2(b, k));
      tmp = ptrio(_,k) * dlocScale_t(xo, df, theta(b, k), sigma) * this_batch ;
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

// works but too slow - use different function/ approach
// deprecate once update_multinomialPrChild works
// [[Rcpp::export]]
Rcpp::NumericMatrix update_multinomialPrOff(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  NumericMatrix ptrio = update_trioPr(xmod) ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;
  IntegerVector nb = model.slot("batchElements") ;
  int N = ptrio.size() ;
  NumericVector this_batch(N) ;
  NumericVector tmp ;
  NumericVector rowtotal(N) ;
  double df = getDf(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  
  Rcpp::NumericVector xo = x[child_ind];
  int M = xo.size() ;
  NumericMatrix poff(M,K) ;
  NumericMatrix lik(M, K) ;
  
  for(int i = 0; i < M; ++i) {
    NumericVector ptrio2 = ptrio(i, _ ) ;
    for(int k = 0; k < K; ++k){
      NumericVector dens(N) ;
      for(int b = 0; b < B; ++b){
        this_batch = batch == ub[b] ;
        double sigma = sqrt(sigma2(b, k));
        tmp = ptrio2[k] * dlocScale_t(xo, df, theta(b, k), sigma) * this_batch ;
        dens += tmp ;
      }
      lik(_, k) = tmp;
      //rowtotal += dens ;
    }
    // NumericMatrix PP(M, K) ;
    //for(int k=0; k<K; ++k){
    //  PP(_, k) = lik(_, k)/rowtotal ;
    //} 
    NumericVector PPI = lik(i,_) ;
    poff.row(i) = PPI ;
  }
  NumericVector sumrow = rowSums(poff);
  for (int j = 0; j < M; ++j) {
    poff(j,_) = poff(j,_) / sumrow[j];
  }
  return poff;
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
  return zc ;
}  

// this selectively updates z for offspring

// [[Rcpp::export]]
Rcpp::IntegerVector update_zchild(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  // Rcpp::IntegerVector map = model.slot("maplabel");
  // Rcpp::S4 hypp(model.slot("hyperparams")) ;
  // int K = getK(hypp) ;
  // Rcpp::IntegerVector sts = getSt(hypp);
  // Rcpp::DataFrame triodat(model.slot("triodata"));
  //int n = triodat.size() ;
  //  NumericVector u = runif(n) ;
  //  IntegerVector ztrio_(n) ;
  //  IntegerVector ztrio = clone(ztrio_) ;
  // updates z for everyone
  //  ztrio = update_z(model) ;
  
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
Rcpp::NumericVector update_mu2(Rcpp::S4 xmod){
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
Rcpp::NumericMatrix update_sigma22(Rcpp::S4 xmod){
  Rcpp::RNGScope scope;
  Rcpp::S4 model(clone(xmod)) ;
  // get model
  // get parameter estimates
  Rcpp::NumericMatrix theta = model.slot("theta");
  Rcpp::IntegerVector z = model.slot("z");
  double nu_0 = model.slot("nu.0");
  double sigma2_0 = model.slot("sigma2.0");
  
  // get data and size attributes
  Rcpp::NumericVector x = model.slot("data");
  int n = x.size();
  int K = theta.ncol();
  int B = theta.nrow();
  NumericVector u = model.slot("u") ;
  double df = getDf(model.slot("hyperparams")) ;
  
  //IntegerVector nn = model.slot("zfreq");
  // get batch info
  Rcpp::NumericMatrix tabz = tableBatchZpar(model);
  Rcpp::IntegerVector batch = model.slot("batch");
  Rcpp::IntegerVector ub = unique_batch(batch);
  Rcpp::NumericMatrix ss(B, K);
  
  for (int i = 0; i < n; ++i) {
    for (int b = 0; b < B; ++b) {
      if (batch[i] != ub[b]) {
        continue;
      }
      for (int k = 0; k < K; ++k){
        if(z[i] == k+1 & batch[i] == b+1){
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

// From stackoverflow http://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder


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
    //model.slot("z") = update_z(model) ;
    //model.slot("zfreq") = tableZ(K, model.slot("z"));
    if(up[7] > 0){
      model.slot("z") = update_zparents(model) ;
      model.slot("zfreq_parents") = tableZpar(model);
      //model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      model.slot("zfreq") = tableZpar(model);
    }
    if(up[8] > 0){
      model.slot("z") = update_zchild(model) ;
    }
    if(up[9] > 0)
      model.slot("pi_parents") = update_pp(model) ;
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
  NumericMatrix pmix_parents = chain.slot("pi_parents") ;
  NumericMatrix zfreq = chain.slot("zfreq") ;
  NumericMatrix zfreq_parents = chain.slot("zfreq_parents") ;
  NumericMatrix mu = chain.slot("mu") ;
  NumericMatrix tau2 = chain.slot("tau2") ;
  NumericVector nu0 = chain.slot("nu.0") ;
  NumericVector sigma2_0 = chain.slot("sigma2.0") ;
  NumericVector loglik_ = chain.slot("loglik") ;
  NumericVector logprior_ = chain.slot("logprior") ;
  //NumericVector th(K) ;
  //NumericVector s2(K) ;
  NumericVector p(K) ;
  NumericVector pp(K) ;
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
  IntegerVector zp(K) ;
  // Initial values
  p = model.slot("pi") ;
  pp = model.slot("pi_parents") ;
  m = model.slot("mu") ;
  t2 = model.slot("tau2") ;
  n0 = model.slot("nu.0") ;
  s20 = model.slot("sigma2.0") ;
  zf = model.slot("zfreq") ;
  zp = model.slot("zfreq_parents") ;
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
  pmix_parents (0, _) = pp ;
  zfreq(0, _) = zf ;
  zfreq_parents(0, _) = zp ;
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
      tmp = model.slot("zfreq_parents") ;
      //tmp = model.slot("zfreq") ;
    }
    zfreq_parents(s, _) = tmp ;
    //zfreq(s, _) = tmp ;
    if(up[8] > 0){
      z = update_zchild(model) ;
      model.slot("z") = z ;
      tmp = tableZ(K, model.slot("z")) ;
      model.slot("zfreq") = tmp ;
      model.slot("probz") = update_probz(model) ;
    } else {
      z = model.slot("z") ;
      tmp = model.slot("zfreq") ;
    }
    zfreq(s, _) = tmp ;
    if(up[9] > 0){
      pp = update_pp(model) ;
      model.slot("pi_parents") = pp ;
    } else {
      pp = model.slot("pi_parents") ;
    }
    pmix_parents(s, _) = pp ;
    if(up[2] > 0){
      p = update_p(model) ;
      model.slot("pi") = p ;
    } else {
      p = model.slot("pi") ;
    }
    pmix(s, _) = p ;
    if(up[0] > 0) {
      model.slot("theta") = update_theta(model) ;
    }
    thetac(s, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
    if(up[1] > 0){
      model.slot("sigma2") = update_sigma2(model) ;
    }
    sigma2c(s, _) = as<Rcpp::NumericVector>(model.slot("sigma2"));
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
        model.slot("z") = update_zparents(model) ;
        model.slot("zfreq_parents") = tableZpar(model) ;
        //tmp = model.slot("zfreq_parents") ;
        //model.slot("zfreq") = tmp ;
      }
      if(up[8] > 0){
        model.slot("z") = update_zchild(model) ;
        model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      }
      if(up[9] > 0)
        model.slot("pi_parents") = update_pp(model) ;
      if(up[2] > 0)
        model.slot("pi") = update_p(model) ;
      if(up[0] > 0)
        model.slot("theta") = update_theta(model) ;
      if(up[1] > 0)
        model.slot("sigma2") = update_sigma2(model) ;
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
  chain.slot("pi_parents") = pmix_parents ;
  chain.slot("mu") = mu ;
  chain.slot("tau2") = tau2 ;
  chain.slot("nu.0") = nu0 ;
  chain.slot("sigma2.0") = sigma2_0 ;
  chain.slot("zfreq") = zfreq ;
  chain.slot("zfreq_parents") = zfreq_parents ;
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
