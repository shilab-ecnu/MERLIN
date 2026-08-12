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



struct ObjGibbs_Gamseo_fixb1{
  vec gamma1;
  vec beta2;
  vec gamma3;
  vec Beta1res;
  vec Beta4res;
  vec Sg12Res;
  vec Sg22Res;
  vec Sg32Res;
  vec s3rs3_temp;
};

ObjGibbs_Gamseo_fixb1 MRGibbs_Gamseo_fixb1(arma::vec &gammah1, arma::vec &gammah3, arma::vec&Gammah1, 
                          arma::mat &S1, arma::mat &S2, arma::mat &S3,
                          arma::mat &L2, arma::mat &U2, arma::mat &R, double &rho, double &b1, int maxIter, int burnin, int thin);










#endif /* GibbsAlpGamEta_ptr_hpp */
