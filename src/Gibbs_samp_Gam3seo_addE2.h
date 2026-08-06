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



struct ObjGibbs_Gam3seo_addE2{
  vec gamma1;
  vec beta2;
  vec Beta1res;
  vec Beta4res;
  vec Sg12Res;
  vec Sg22Res;
  vec Sg32Res;
  vec s3rs3_temp;
};

ObjGibbs_Gam3seo_addE2 MRGibbs_Gam3seo_addE2(arma::vec &gammah1,arma::vec &gammah3, arma::vec &Gammah1, arma::vec &Gammah3, 
                                             arma::vec &se1, arma::vec &se2, arma::vec &se3, arma::vec &se4, arma::mat &R, 
                                             double &rho_1, double &rho_2, int maxIter, int burnin, int thin);









#endif /* GibbsAlpGamEta_ptr_hpp */
