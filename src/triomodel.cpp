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
  return family_member;
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
  is_parental_cn = (ind==TRUE) & (ind2==TRUE) ;
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
  // TODO: repeated with each MCMC
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
Rcpp::LogicalVector is_father(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  CharacterVector fam = family_member(model);
  int N = fam.size() ;
  LogicalVector is_f(N) ;
  for (int i = 0; i < N; i++){
    is_f[i] = (fam[i] == "f");
  }
  return is_f ;
}

// [[Rcpp::export]]
Rcpp::LogicalVector is_mother(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  CharacterVector fam = family_member(model);
  int N = fam.size() ;
  LogicalVector is_f(N) ;
  for (int i = 0; i < N; i++){
    is_f[i] = (fam[i] == "m");
  }
  return is_f ;
}

// [[Rcpp::export]]
Rcpp::LogicalVector is_child(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  CharacterVector fam = family_member(model);
  int N = fam.size() ;
  LogicalVector is_f(N) ;
  for (int i = 0; i < N; i++){
    is_f[i] = (fam[i] == "o");
  }
  return is_f ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix update_trioPr2(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::DataFrame triodat(model.slot("triodata"));
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  IntegerVector is_mendelian=model.slot("is_mendelian") ;
  int K = getK(hypp) ;
  IntegerVector z = model.slot("z");
  NumericVector p = model.slot("pi");
  CharacterVector fam = family_member(model);
  LogicalVector is_f=is_father(model) ;
  LogicalVector is_m=is_mother(model) ;
  //LogicalVector is_o=is_child(model) ;
  Rcpp::IntegerVector zf = z[is_f];
  Rcpp::IntegerVector zm = z[is_m];
  int T = zf.size();
  Rcpp::NumericMatrix zo_prob(T, K);
  Rcpp::NumericVector probs(K);
  LogicalVector is_mendel(T);
  is_mendel = is_mendelian == 1 ;
  for (int i = 0; i < T; i++){
    if(is_mendel[i]){
      probs = lookup_mprobs(model, zf[i], zm[i]) ;
    } else {
      probs = p ;
    }
    zo_prob(i, _) = probs ;
  }
  return zo_prob;
}

// [[Rcpp::export]]
Rcpp::IntegerVector update_mendelian(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  IntegerVector z = model.slot("z");
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  Rcpp::LogicalVector fat_ind(fam.size());
  Rcpp::LogicalVector mat_ind(fam.size());
  // TODO: repeated with each MCMC
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  Rcpp::IntegerVector zo = z[child_ind];
  // TODO: prior probability for mendelian indicator
  double m_prior = 0.9 ;
  double numer;
  double denom;
  NumericMatrix ptrio = update_trioPr(model) ;
  NumericVector p = model.slot("pi_parents") ;
  int T=zo.size() ;
  IntegerVector mendel(T);
  double prob_mendel ;
  NumericVector u(1) ;
  int cn;
  // Pr(mendel = 1 |...) = p(z_0 | z_m, z_f, M=1) x
  //                       p(M=1)/(p(z_0 | z_m, z_f, M=1)P(M=1) +
  //                       p(z_0 | M=0)P(M=0))
  for(int i=0; i < T; ++i){
    cn = zo[i] - 1;
    numer=ptrio(i, _)[ cn ] * m_prior;
    denom=numer + p[ cn ] * (1-m_prior) ;
    prob_mendel = numer/denom ;
    u = runif(1) ;
    if(u[0] <= prob_mendel){
      mendel[i] = 1 ;
    } else {
      mendel[i] = 0 ;
    }
  }
  return mendel;
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
  IntegerVector nb = model.slot("batchElements") ;;
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

  for(int k = 0; k < K; ++k){
    NumericVector dens(M) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      double sigma = sqrt(sigma2(b, k));
      //tmp = p[k] * pp[k] * dlocScale_t(xp, df, theta(b, k), sigma) * this_batch ;
      //tmp = ((p[k]+pp[k])/2) * dlocScale_t(xp, df, theta(b, k), sigma) * this_batch ;
      tmp = p[k] * dlocScale_t(xp, df, theta(b, k), sigma) * this_batch ;
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



//NumericVector mendelprobs(Rcpp::S4 xmod, int father, int mother, int offspring, NumericVector trioprob){
//  RNGScope scope ;
//  Rcpp::S4 model(clone(xmod)) ;
//  NumericVector p = model.slot("pi") ;
//  NumericVector pp = model.slot("pi_parents") ;
//  NumericMatrix ptrio = update_trioPr(xmod) ;
//  double fatp = p[father-1] ;
//  double motp = p[mother-1] ;
//  double offp = p[offspring-1] ;
//  double fatpp = pp[father-1] ;
//  double motpp = pp[mother-1] ;
//  NumericVector beta = Rcpp::rbeta(1, 10, 1) ;
//  double beta1 = beta[0] ;
//  double denom = (fatpp * motpp * offp * beta1) + (fatp * motp * offp * (1-beta1)) ;
//  NumericVector mendel = (trioprob * fatpp * motpp * fatp * motp * beta1) / denom ;
//  return mendel;
//}


//NumericMatrix update_mendelPr(Rcpp::S4 xmod){
//  RNGScope scope ;
//  Rcpp::S4 model(clone(xmod)) ;
//  Rcpp::S4 hypp(model.slot("hyperparams")) ;
//  int K = getK(hypp) ;
//  IntegerVector z = model.slot("z");
//  NumericMatrix ptrio = update_trioPr(xmod) ;
//  CharacterVector fam = family_member(xmod);
//  Rcpp::LogicalVector fat_ind(fam.size());
//  Rcpp::LogicalVector mat_ind(fam.size());
//  Rcpp::LogicalVector off_ind(fam.size());
//  for (int i = 0; i < fam.size(); i++){
//    fat_ind[i] = (fam[i] == "f");
//    mat_ind[i] = (fam[i] == "m");
//    off_ind[i] = (fam[i] == "o");
//  }
//  Rcpp::IntegerVector zf = z[fat_ind];
//  Rcpp::IntegerVector zm = z[mat_ind];
//  Rcpp::IntegerVector zo = z[off_ind];
//  int trio_size = zf.size();
//  Rcpp::NumericVector mendel_probs;
//  Rcpp::NumericMatrix mendelmat(trio_size, K);
//  NumericVector ptrios;
//  for (int i = 0; i < trio_size; i++){
//   ptrios = ptrio(i,_) ;
//   mendel_probs = mendelprobs(xmod, zf[i], zm[i], zo[i], ptrios);
//   mendelmat(i,_) = mendel_probs;
//  }
//  return mendelmat;
//}

// double mendelprobs2(Rcpp::S4 xmod, int father, int mother, int offspring){
//   RNGScope scope ;
//   Rcpp::S4 model(clone(xmod)) ;
//   NumericVector p = model.slot("pi") ;
//   NumericVector pp = model.slot("pi_parents") ;
//   NumericMatrix ptrio = update_trioPr(xmod) ;
//   double fatp = p[father-1] ;
//   double motp = p[mother-1] ;
//   double offp = p[offspring-1] ;
//   double fatpp = pp[father-1] ;
//   double motpp = pp[mother-1] ;
//   double offpp = pp[offspring-1] ;
//   NumericVector beta = Rcpp::rbeta(1, 10, 1) ;
//   double beta1 = beta[0] ;
//   double denom = (fatpp * motpp * offpp * beta1) + (fatp * motp * offp * (1-beta1)) ;
//   double mendel = (offpp * fatpp * motpp * beta1) / denom ;
//   return mendel;
// }
// 
// NumericMatrix update_mendelPr2(Rcpp::S4 xmod){
//   RNGScope scope ;
//   Rcpp::S4 model(clone(xmod)) ;
//   Rcpp::S4 hypp(model.slot("hyperparams")) ;
//   //int K = getK(hypp) ;
//   IntegerVector z = model.slot("z");
//   //NumericMatrix ptrio = update_trioPr(xmod) ;
//   CharacterVector fam = family_member(xmod);
//   Rcpp::LogicalVector fat_ind(fam.size());
//   Rcpp::LogicalVector mat_ind(fam.size());
//   Rcpp::LogicalVector off_ind(fam.size());
//   for (int i = 0; i < fam.size(); i++){
//     fat_ind[i] = (fam[i] == "f");
//     mat_ind[i] = (fam[i] == "m");
//     off_ind[i] = (fam[i] == "o");
//   }
//   Rcpp::IntegerVector zf = z[fat_ind];
//   Rcpp::IntegerVector zm = z[mat_ind];
//   Rcpp::IntegerVector zo = z[off_ind];
//   int trio_size = zf.size();
//   Rcpp::NumericVector mendel_probs;
//   Rcpp::NumericMatrix mendelmat(trio_size,1);
//   for (int i = 0; i < trio_size; i++){
//     mendel_probs = mendelprobs2(xmod, zf[i], zm[i], zo[i]);
//     mendelmat(i,_) = mendel_probs;
//   }
//   return mendelmat;
// }

// [[Rcpp::export]]
Rcpp::NumericMatrix update_multinomialPrChild(Rcpp::S4 xmod) {
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = getK(hypp) ;
  IntegerVector batch = model.slot("batch") ;
  IntegerVector ub = unique_batch(batch) ;
  // Mendelian transmission probability matrix
  NumericMatrix ptrio = update_trioPr2(xmod) ;
  NumericVector p = model.slot("pi") ;
  // commented by RS: pp is not unused
  //NumericVector pp = model.slot("pi_parents") ;
  NumericMatrix sigma2 = model.slot("sigma2") ;
  NumericMatrix theta = model.slot("theta") ;
  int B = sigma2.nrow() ;
  NumericVector x = model.slot("data") ;
  IntegerVector nb = model.slot("batchElements") ;
  double df = getDf(hypp) ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  // TODO:
  // double p_mendel=0.9;
  // TODO: we do this with each MCMC iteration
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  Rcpp::NumericVector xo = x[child_ind];
  // Mendelian observations
  int M = xo.size() ;
  NumericMatrix lik(M, K) ;
  NumericVector this_batch(M) ;
  NumericVector tmp(M) ;
  NumericVector tmp2(M) ;
  NumericVector rowtotal(M) ;
  // MC:  why do we multiply ptrio by p[k]?
  for(int k = 0; k < K; ++k){
    NumericVector dens(M) ;
    for(int b = 0; b < B; ++b){
      this_batch = batch == ub[b] ;
      //tmp = p[k] * ptrio(_,k) *
      //      dlocScale_t(xo, df, theta(b, k), sigma) * this_batch ;
      double sigma = sqrt(sigma2(b, k));
      NumericVector phi=dlocScale_t(xo, df, theta(b, k), sigma);
      tmp = ptrio(_, k) * phi * this_batch ;
      //tmp = ptrio(_, k) * phi * this_batch * (1 - p_mendel) ;
      //tmp2 = p[k] * phi * this_batch * p_mendel ;
      //tmp = tmp + tmp2 ;
      dens += tmp ;
    }
    lik(_, k) = dens;
    rowtotal += dens ;
  }
  NumericMatrix PC(M, K) ;
  for(int k=0; k < K; ++k){
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
  Rcpp::DataFrame triodat(model.slot("triodata"));
  NumericMatrix theta = model.slot("theta") ;
  IntegerVector batch = model.slot("batch") ;
  int B = theta.nrow() ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  // TODO: This is unchanging, yet we do it with each MCMC iteration
  for (int i = 0; i < fam.size(); i++){
    child_ind[i] = (fam[i] == "o");
  }
  IntegerVector z = model.slot("z");
  Rcpp::IntegerVector zo = z[child_ind];
  int child_size = zo.size();
  NumericMatrix p(child_size, K);
  p = update_multinomialPrChild(xmod) ;
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


// [[Rcpp::export]]
Rcpp::IntegerVector update_zchild(Rcpp::S4 xmod) {
  RNGScope scope ;
  S4 model(clone(xmod)) ;
  IntegerVector ztrio = model.slot("z");
  IntegerVector zz_offspring ;
  CharacterVector fam = family_member(xmod);
  LogicalVector child_ind(fam.size());
  // TODO:  This is unchanging, yet we do it every MCMC iteration
  for (int i = 0; i < fam.size(); i++){
     child_ind[i] = (fam[i] == "o");
  }
  //
  zz_offspring = update_offspring(model); // length 300
  int n = ztrio.size() ;
  int j = 0;
  for(int i = 0; i < n; i++){
    if(child_ind[i] == TRUE){
      ztrio[i] = zz_offspring[j];
      ++j;
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
Rcpp::IntegerMatrix update_probzpar(Rcpp::S4 xmod){
  RNGScope scope ;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::S4 hypp(model.slot("hyperparams")) ;
  int K = hypp.slot("k") ;
  CharacterVector fam = family_member(xmod);
  Rcpp::LogicalVector child_ind(fam.size());
  // TODO:  this is re-computed at every MCMC iteration
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
        if((z[i] == k+1) & (batch[i] == b+1)){
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

//[[Rcpp::export]]
Rcpp::IntegerVector sample_trio_components(Rcpp::IntegerVector x, int size, Rcpp::NumericVector prob){
  int n = x.size() ;
  Rcpp::IntegerVector z = clone(x);
  for(int i=0; i < n; i++){
    //initialize accumulator ;
    double accept = 0 ;
    Rcpp::NumericVector u=runif(1);
    for(int j = 0; j < n; j++){
      accept += prob[j] ;
      if( u[0] < accept ) {
        z[i] = x[j] ;
        break ;
      }
    }
  }
  return(z);
}

//[[Rcpp::export]]
Rcpp::S4 predictive_trios(Rcpp::S4 xmod){
  Rcpp::RNGScope scope;
  Rcpp::S4 model(clone(xmod)) ;
  Rcpp::NumericMatrix theta = model.slot("theta");
  Rcpp::NumericMatrix sigma2 = model.slot("sigma2");
  Rcpp::NumericVector prob = model.slot("pi");
  Rcpp::IntegerVector batch = model.slot("batch") ;
  int K = theta.ncol();
  int B = theta.nrow();
  Rcpp::IntegerVector components=seq_len(K);
  double df = getDf(model.slot("hyperparams")) ;
  Rcpp::IntegerVector z(K);
  Rcpp::NumericMatrix ystar(B, K) ;
  Rcpp::IntegerMatrix zstar(B, K) ;
  Rcpp::NumericVector u=Rcpp::rchisq(K*B, df) ;;
  // sample components according to mixture probabilities
  // mixture probabilities are assumed to be the same for each batch
  z=sample_trio_components(components, K, prob);
  z=z-1;  // Need to subtract 1 for this to index the right column
  Rcpp::NumericVector yy(K*B);
  Rcpp::IntegerVector zz(K*B);
  double sigma;
  int index;
  int j=0;
  for(int k=0; k < K; ++k){
    for(int b = 0; b < B; b++){
      index = z[k];
      sigma = sqrt(sigma2(b, index)) ;
      zstar(b, k) = index ;
      ystar(b, k) = (rlocScale_t(1, theta(b, index), sigma, df, u[j]))[0];
      yy[j] = ystar(b, k);
      zz[j] = zstar(b, k);
      j++;
    }
  }
  model.slot("predictive") = yy;
  model.slot("zstar") = zz ;
  return model ;
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
  for(int s = 1; s < S; ++s){
    model.slot("z") = update_z(model) ;
    model.slot("zfreq_parents") = tableZpar(model) ;
    model.slot("sigma2") = update_sigma2(model) ;
    model.slot("nu.0") = update_nu0(model) ;
    model.slot("sigma2.0") = update_sigma20(model) ;
    model.slot("z") = update_zchild(model) ;
    model.slot("is_mendelian") = update_mendelian(model) ;
    model.slot("zfreq") = tableZ(K, model.slot("z")) ;
    model.slot("theta") = update_theta(model) ;
    model.slot("mu") = update_mu(model) ;
    model.slot("tau2") = update_tau2(model) ;
    model.slot("pi_parents") = update_pp(model) ;
    model.slot("pi") = update_p(model) ;
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
}

// [[Rcpp::export]]
Rcpp::IntegerVector test_trio(Rcpp::S4 object) {
  RNGScope scope ;
  Rcpp::S4 model(clone(object)) ;
  IntegerVector is_mendel=model.slot("is_mendelian");
  Rcpp::S4 chains=model.slot("mcmc.chains") ;
  IntegerVector mendelian=chains.slot("is_mendelian");
  mendelian = mendelian + is_mendel  + is_mendel;
  return mendelian;
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
  S = S - 1;
  T = T - 1;
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
  NumericMatrix predictive_ = chain.slot("predictive") ;
  IntegerMatrix zstar_ = chain.slot("zstar") ;
  IntegerVector mendelian_ = chain.slot("is_mendelian") ;
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
  int ntrio = mendelian_.size();
  IntegerVector temp( ntrio ) ;
  IntegerVector zf(K) ;
  IntegerVector zp(K) ;
  NumericVector ystar = NumericVector(B*K);
  IntegerVector zstar = IntegerVector(B*K);
  for(int s = 0; s < (S + 1); ++s){
    // update z for every one (independence)
    z = update_z(model) ;
    model.slot("z") = z ;
    // z frequency of parents
    tmp = tableZpar(model) ;
    model.slot("zfreq_parents") = tmp ;
    zfreq_parents(s, _) = tmp ;
    // updates integer matrix of slot probz for only the parents
    model.slot("probz_par") = update_probzpar(model) ;
    // updates z slot only for the offspring
    z = update_zchild(model) ;
    model.slot("z") = z ;
    model.slot("probz") = update_probz(model) ;
    tmp = tableZ(K, model.slot("z")) ;
    // updates integer matrix of slot probz for all individuals
    model.slot("zfreq") = tmp ;
    zfreq(s, _) = tmp ;
    temp = update_mendelian(model) ;
    model.slot("is_mendelian") = temp ;
    mendelian_ = mendelian_ + temp ;

    model.slot("sigma2") = update_sigma2(model) ;
    sigma2c(s, _) = as<Rcpp::NumericVector>(model.slot("sigma2"));
    n0 = update_nu0(model) ;
    model.slot("nu.0") = n0 ;
    nu0[s] = n0[0] ;
    s20 = update_sigma20(model) ;
    model.slot("sigma2.0") = s20 ;
    sigma2_0[s] = s20[0] ;
    model.slot("theta") = update_theta(model) ;
    thetac(s, _) = as<Rcpp::NumericVector>(model.slot("theta")) ;
    t2 = update_tau2(model) ;
    model.slot("tau2") = t2 ;
    tau2(s, _) = t2 ;
    m = update_mu(model) ;
    model.slot("mu") = m ;
    mu(s, _) = m ;
    pp = update_pp(model) ;
    model.slot("pi_parents") = pp ;
    //pmix_parents(s, _) = pp ;
    p = update_p(model) ;
    model.slot("pi") = p ;
    pmix(s, _) = p ;
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
    model = predictive_trios(model);
    ystar = model.slot("predictive");
    zstar = model.slot("zstar");
    predictive_(s, _) = ystar ;
    zstar_(s, _) = zstar ;
    // Thinning
    for(int t = 0; t < T; ++t){
      model.slot("z") = update_z(model) ;
      model.slot("zfreq_parents") = tableZpar(model) ;
      model.slot("sigma2") = update_sigma2(model) ;
      model.slot("nu.0") = update_nu0(model) ;
      model.slot("sigma2.0") = update_sigma20(model) ;
      model.slot("z") = update_zchild(model) ;
      model.slot("is_mendelian") = update_mendelian(model) ;
      model.slot("zfreq") = tableZ(K, model.slot("z")) ;
      model.slot("theta") = update_theta(model) ;
      model.slot("tau2") = update_tau2(model) ;
      model.slot("mu") = update_mu(model) ;
      model.slot("pi_parents") = update_pp(model) ;
      model.slot("pi") = update_p(model) ;
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
  chain.slot("is_mendelian") = mendelian_ ;
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
