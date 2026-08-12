#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "Gibbs_samp_Gamseo_fixb1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

inline double Updatesig2(arma::vec &gamk, double ak,double bk, int &p ){
    double invga,invgb,sig2;
    invga = ak + p / 2;
    invgb = as_scalar(gamk.t()*gamk) / 2 + bk;

    //Rprintf("Have a look of a = %4f \n",invga);
    //Rprintf("Have a look of b = %4f \n",invgb);
    sig2 =  1 / randg<double>(distr_param(invga, 1/invgb));
    //Rprintf("Have a look of siga = %4f \n",invga);

    return sig2;

}

// [[Rcpp::depends(RcppArmadillo)]]

ObjGibbs_Gamseo_fixb1 MRGibbs_Gamseo_fixb1(arma::vec &gammah1, arma::vec &gammah3, arma::vec&Gammah1, 
                          arma::mat &S1, arma::mat &S2, arma::mat &S3,
                          arma::mat &L2, arma::mat &U2, arma::mat &R, double &rho, double &b1,
                          int maxIter, int burnin, int thin)
{

//initial valiues---------------------------
int p = gammah1.n_elem;

double sigma12 = 1; double sigma22 = 1;double sigma32 = 1;
double beta1 = 0.1;double beta4 = 0.1;
double ak1 = 0.001 ; double bk1 =0.001;
double ak2 = 0.001 ; double bk2 =0.001;
double ak3 = 0.001 ; double bk3 =0.001;

vec gamma1 = 0.01*ones(p, 1);
vec gamma3 = 0.01*ones(p, 1);
vec beta2 = 0.01*ones(p, 1);

mat s1rs1 = inv(S1) * R * inv(S1);
mat s2rs2 = inv(S2) * R * inv(S2);
mat s3rs3 = inv(S3) * R * inv(S3);
mat s1rs3 = inv(S1) * R * inv(S3);
vec s1rs1_temp = diagvec(s1rs1);
vec s3rs3_temp = diagvec(s3rs3);
vec s1rs3_temp = diagvec(s1rs3);
vec s2rs2_temp = diagvec(s2rs2);

vec g3_til = inv(L2) * gammah3;
vec s1s1g1_hat = inv(S1) * inv(S1) * gammah1;
vec s1s3g1_hat = inv(S1) * inv(S3) * gammah1;
vec s3s3G1_hat = inv(S3) * inv(S3) * Gammah1;
vec s1s3G1_hat = inv(S1) * inv(S3) * Gammah1;
vec s2s2g3_hat = inv(S2) * inv(S2) * gammah3;

int numsave = maxIter / thin;

// for calculate prepared-----------------------

vec v12 = zeros(p, 1);
vec v32 = zeros(p, 1);
vec vb22 = zeros(p, 1);
double vb12 = 0.0;
double vb42 = 0.0;

vec mu11 = zeros(p, 1);
vec mu33 = zeros(p, 1);
vec mub22 = zeros(p, 1);
double mub11 = 0.0;
double mub44 = 0.0;
double sig12 = 1.0;
double sig22 = 1.0;
double sig32 = 1.0;

int l = 0;

//prepare for res para--------------------------
vec Beta1res = ones(numsave, 1);
vec Beta4res = ones(numsave, 1);
vec Sg12Res = ones(numsave, 1);
vec Sg22Res = ones(numsave, 1);
vec Sg32Res = ones(numsave, 1);


double rho1 = 1 / (1 - rho * rho);
for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin) ; iter ++){
//for(int iter = 0; iter < 500 ; iter ++){
//-----------------------------------Begin to give some middle value-----------

  double invsg12 = 1. / sigma12;
  double invsg22 = 1. / sigma22;
  double invsg32 = 1. / sigma32;

  double b12 = beta1*beta1;
  double b42 = beta4*beta4;

 for(int j = 0; j < p; j++){
    //update gamma1-----------------------------------------------
    v12[j] = 1. / (rho1 * s1rs1_temp[j] - 2 * rho * rho1 * beta1 * s1rs3_temp[j] + rho1 * b12 * s3rs3_temp[j] + invsg12) ;
//    Rprintf("-------------------------------------------------------------------------------Have a look of v12 = %4f \nsgx", v12[j] );

    double mu1;
   //Rprintf("Have a look of G1tu3 = %4f \n",G1tu3);
    mu1 = rho1 * (- (dot(s1rs1.row(j), gamma1) - gamma1[j] * s1rs1_temp[j] )
                  + 2 * rho * beta1 * (dot(s1rs3.row(j), gamma1) - gamma1[j] * s1rs3_temp[j])
                  - b12 * (dot(s3rs3.row(j), gamma1) - gamma1[j] * s3rs3_temp[j] )
                  - rho * beta1 * s1s3g1_hat[j] + beta1 * s3s3G1_hat[j]
                  + s1s1g1_hat[j] - rho * s1s3G1_hat[j]
                  + rho * dot(s1rs3.row(j), beta2) + rho * beta4 * dot(s1rs3.row(j), gamma3)
                  - beta1 * dot(s3rs3.row(j), beta2) - beta1 * beta4 * dot(s3rs3.row(j), gamma3)
                 ) * v12[j];
    mu11[j] = mu1 + randn()*sqrt(v12[j]);
    gamma1[j] = mu11[j];
    //gamma1[j] = 0 * mu11[j];
    //Rprintf("Have a look of mu1_p3 = %4f \n", gamma1[j] * s3rs3_temp[j] );
}

 for(int j = 0; j < p; j++){
      //update  gamma3-----------------------------------------------------
    v32[j] = 1. / (b42 * rho1 * s3rs3_temp[j] + s2rs2_temp[j] + invsg32) ;     
    double mu3;
    mu3 = rho1 * (- b42 * (dot(s3rs3.row(j), gamma3) - gamma3[j] * s3rs3_temp[j] )
                  - rho * beta4 * s1s3g1_hat[j] + beta4 * s3s3G1_hat[j]
                  + rho * beta4 * dot(s1rs3.row(j), gamma1) - beta1 * beta4 * dot(s3rs3.row(j), gamma1) - beta4 * dot(s3rs3.row(j), beta2)
                 ) * v32[j]
                  + s2s2g3_hat[j] * v32[j]
                  - (dot(s2rs2.row(j), gamma3) - gamma3[j] * s2rs2_temp[j] ) * v32[j];
    mu33[j] = mu3 + randn()*sqrt(v32[j]);
    //Rprintf("Have a look of gamma3 = %4f \n", mu3 );
    gamma3[j] = mu33[j];
    //gamma3[j] = 0 * gamma3[j];
    
}
 for(int j = 0; j < p; j++){
      //update beta2--------------------------------------------------
    vb22[j] = 1. / (rho1 * s3rs3_temp[j] + invsg22) ;     
    double mub2;
    mub2 = rho1 * (- (dot(s3rs3.row(j), beta2) - beta2[j] * s3rs3_temp[j] )
                   - rho * s1s3g1_hat[j] + s3s3G1_hat[j]
                   + rho * dot(s1rs3.row(j), gamma1) - beta1 * dot(s3rs3.row(j), gamma1) - beta4 * dot(s3rs3.row(j), gamma3)
                  ) * vb22[j];
    mub22[j] = mub2 + randn()*sqrt(vb22[j]);
    beta2[j] = mub22[j];
    //Rprintf("Have a look of beta2 = %4f \n", - beta1 * dot(s3rs3.row(j), gamma1) );
    //beta2[j] = 0 * beta2[j];
}
    //update beta1-----------------------------------------------------------
    vb12 = 1. / as_scalar( rho1 * gamma1.t() * s3rs3 * gamma1);
    double mub1;
    mub1 = rho1 * as_scalar(-rho * gamma1.t() * s1s3g1_hat + gamma1.t() * s3s3G1_hat
                            + rho * gamma1.t() * s1rs3 * gamma1 - gamma1.t() * s3rs3 * beta2 - beta4 * gamma1.t() * s3rs3 * gamma3
                           )* vb12 ;
    mub11 = mub1 + randn()*sqrt(vb12);
    beta1 = b1;  
    //beta1 = 0.0;
    //Rprintf("Have a look of b1 = %4f \n", beta1 );

    //update beta4-----------------------------------------------------------
    vb42 = 1. / as_scalar( rho1 * gamma3.t() * s3rs3 * gamma3);
    //Rprintf("Have a look of vb12 = %4f \n",vb12);
    double mub4;
    mub4 = rho1 * as_scalar(-rho* gamma3.t() * s1s3g1_hat + gamma3.t() * s3s3G1_hat
                            + rho * gamma3.t() * s1rs3 * gamma1 - beta1 * gamma1.t() * s3rs3 * gamma3 - beta2.t() * s3rs3 * gamma3
                           ) * vb42 ;
    mub44 = mub4 + randn()*sqrt(vb42);
    beta4 = mub44;
    //beta4 = 0.0;  
    //update sigma123------------------------------------------------------------
    sig12 = Updatesig2(gamma1, ak1, bk1, p );
    sigma12 = sig12;
    //sigma12 = sg1;
    sig22 = Updatesig2(beta2, ak2, bk2, p );
    sigma22 = sig22;
    //sigma22 = sg2;
    sig32 = Updatesig2(gamma3, ak3, bk3, p );
    sigma32 = sig32;
    //sigma32 = sg3;
    //prepare for next iter------------------------------------------------------------------------------------------

 if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
     Sg12Res[l] = sigma12;
        Sg22Res[l] = sigma22;
        Sg32Res[l] = sigma32;
        Beta1res[l] = beta1;
        Beta4res[l] = beta4;
        l += 1;
      }
    }

}
ObjGibbs_Gamseo_fixb1 obj;
  obj.gamma1 = gamma1;
  obj.beta2 = beta2;
  obj.gamma3 = gamma3;
  obj.Beta1res = Beta1res;
  obj.Beta4res = Beta4res;
  obj.Sg12Res = Sg12Res;
  obj.Sg22Res = Sg22Res;
  obj.Sg32Res = Sg32Res;
  obj.s3rs3_temp = s3rs3_temp;
  return obj;
}



//[[Rcpp::export]]
List MRGEI_Gamseo_fixb1(arma::vec &gammah1,arma::vec &gammah3, arma::vec &Gammah1, 
                        arma::vec &se1, arma::vec &se2, arma::vec &se3, 
                        arma::mat &R, double &rho, double &b1,
                        int maxIter, int burnin, int thin)
{
    //mat R_sh = cal_blockcor(R) ;
    
    mat S1 = diagmat(se1);
    mat S2 = diagmat(se2);
    mat S3 = diagmat(se3);

    Rprintf("Have a look1");
    mat L2 = chol( S2 * R * S2, "lower");
    mat U2 = inv(L2) * S2 * R * inv(S2);
    Rprintf("Have a look2");
    ObjGibbs_Gamseo_fixb1 obj = MRGibbs_Gamseo_fixb1( gammah1, gammah3, Gammah1, S1, S2, S3, L2, U2, R, rho, b1, maxIter, burnin, thin);

    double bhat1 = b1;
    double b1se1 = 0;
    double pvalue1 = 0;

    double bhat4 = mean(obj.Beta4res);
    double b4se4 = stddev(obj.Beta4res);
    double pvalue4 = 2 * (R::pnorm(abs(bhat4 / b4se4),0,1,0,0));


 // ObjGibbs3 obj = MRcorrObj(gammah, Gammah, se1, se2, lp_opt);
  List output = List::create(
    Rcpp::Named("Beta1.hat") =bhat1,
    Rcpp::Named("Beta1.se") =b1se1,
    Rcpp::Named("Beta1.pval") = pvalue1,
    Rcpp::Named("Beta4.hat") =bhat4,
    Rcpp::Named("Beta4.se") =b4se4,
    Rcpp::Named("Beta4.pval") = pvalue4,
    Rcpp::Named("gamma1") = Rcpp::wrap(obj.gamma1),
    Rcpp::Named("beta2") = Rcpp::wrap(obj.beta2),
    Rcpp::Named("gamma3") = Rcpp::wrap(obj.gamma3),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Beta4res") = Rcpp::wrap(obj.Beta4res),
    Rcpp::Named("Sg12Res") = Rcpp::wrap(obj.Sg12Res),
    Rcpp::Named("Sg22Res") = Rcpp::wrap(obj.Sg22Res),
    Rcpp::Named("Sg32Res") = Rcpp::wrap(obj.Sg32Res),
    Rcpp::Named("S3RS3") = Rcpp::wrap(obj.s3rs3_temp),
    Rcpp::Named("se2") = se2,
    Rcpp::Named("L2") = L2,
    Rcpp::Named("U2") = U2

  );
  return output;
}

