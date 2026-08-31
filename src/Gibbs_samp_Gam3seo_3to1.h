#ifndef Gibbs_samp_hpp
#define Gibbs_samp_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>

 
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]



struct ObjGibbs_Gam3seo_3to1{
  vec gamma1;
  vec beta2;
  vec gamma3;
  vec Beta1res;
  vec Beta4res;
  vec Sg12Res;
  vec Sg22Res;
  vec Sg32Res;
};

ObjGibbs_Gam3seo_3to1 MRGibbs_Gam3seo_3to1(
    arma::vec &gammah1,
    arma::vec &gammah3,
    arma::vec &Gammah1,
    arma::vec &Gammah3,
    arma::mat &invS1,
    arma::mat &invS2,
    arma::mat &invS3,
    arma::mat &invS4,
    arma::mat &R,
    double &rho_1,
    double &rho_2,
    double &mup,
    int maxIter,
    int burnin,
    int thin
);




#endif /* GibbsAlpGamEta_ptr_hpp */