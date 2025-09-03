#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "Gibbs_samp_Gam3seo.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

double Updatesig2(arma::vec &gamk, double ak,double bk, int &p ){
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

ObjGibbs_Gam3seo MRGibbs_Gam3seo(arma::vec &gammah1, arma::vec &gammah3, arma::vec&Gammah1, arma::vec&Gammah3,
                          arma::mat &invS1, arma::mat &invS2, arma::mat &invS3, arma::mat &invS4,
                          arma::mat &R, double &rho_1, double &rho_2)
{
    int maxIter = 12000;
    int thin = 10;
    int burnin = 5000;
    int p = gammah1.n_elem;

    double sigma12 = 1; double sigma22 = 1;double sigma32 = 1;
    double beta1 = 0.1;double beta4 = 0.1;
    double ak1 = 0.001 ; double bk1 =0.001;
    double ak2 = 0.001 ; double bk2 =0.001;
    double ak3 = 0.001 ; double bk3 =0.001;

    vec gamma1 = 0.01*ones(p, 1);
    vec gamma3 = 0.01*ones(p, 1);
    vec beta2 = 0.01*ones(p, 1);

    mat s1rs1 = invS1 * R * invS1;
    mat s3rs3 = invS3 * R * invS3;
    mat s1rs3 = invS1 * R * invS3;

    mat s2rs2 = invS2 * R * invS2;
    mat s4rs4 = invS4 * R * invS4;
    mat s2rs4 = invS2 * R * invS4;
    vec s1rs1_temp = diagvec(s1rs1);
    vec s3rs3_temp = diagvec(s3rs3);
    vec s1rs3_temp = diagvec(s1rs3);
    vec s2rs2_temp = diagvec(s2rs2);
    vec s4rs4_temp = diagvec(s4rs4);
    vec s2rs4_temp = diagvec(s2rs4);


    vec s1s1g1_hat = invS1 * invS1 * gammah1;
    vec s1s3g1_hat = invS1 * invS3 * gammah1;
    vec s3s3G1_hat = invS3 * invS3 * Gammah1;
    vec s1s3G1_hat = invS1 * invS3 * Gammah1;

    vec s4s4G3_hat = invS4 * invS4 * Gammah3;
    vec s2s4g3_hat = invS2 * invS4 * gammah3;
    vec s2s2g3_hat = invS2 * invS2 * gammah3;
    vec s2s4G3_hat = invS2 * invS4 * Gammah3;

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

    double rho1 = 1 / (1 - rho_1 * rho_1);
    double rho2 = 1 / (1 - rho_2 * rho_2);
    for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin) ; iter ++){
    //------------Begin to give some middle value-----------
        double invsg12 = 1. / sigma12;
        double invsg22 = 1. / sigma22;
        double invsg32 = 1. / sigma32;
        double b12 = beta1*beta1;
        double b42 = beta4*beta4;
        //update gamma1
        for(int j = 0; j < p; j++){
            v12[j] = 1. / (rho1 * s1rs1_temp[j] - 2 * rho_1 * rho1 * beta1 * s1rs3_temp[j] + rho1 * b12 * s3rs3_temp[j] 
                   + rho2 * b42 * s4rs4_temp[j] + invsg12) ;
            double mu1;
            mu1 =  rho1 * (- (dot(s1rs1.row(j), gamma1) - gamma1[j] * s1rs1_temp[j] )
                   + 2 * rho_1 * beta1 * (dot(s1rs3.row(j), gamma1) - gamma1[j] * s1rs3_temp[j])
                   - b12 * (dot(s3rs3.row(j), gamma1) - gamma1[j] * s3rs3_temp[j] )
                   - rho_1 * beta1 * s1s3g1_hat[j] + beta1 * s3s3G1_hat[j]
                   + s1s1g1_hat[j] - rho_1 * s1s3G1_hat[j]
                   + rho_1 * dot(s1rs3.row(j), beta2) + rho_1 * beta4 * dot(s1rs3.row(j), gamma3)
                   - beta1 * dot(s3rs3.row(j), beta2) - beta1 * beta4 * dot(s3rs3.row(j), gamma3)
                  ) * v12[j]
                    +rho2 * (- b42 * (dot(s4rs4.row(j), gamma1) - gamma1[j] * s4rs4_temp[j] )
                   - rho_2 * beta4 * s2s4g3_hat[j] + beta4 * s4s4G3_hat[j]
                   - beta1 * beta4 * dot(s4rs4.row(j), gamma3) + rho_2 * beta4 * dot(s2rs4.row(j), gamma3)
                  ) * v12[j] ;
            mu11[j] = mu1 + randn()*sqrt(v12[j]);
            gamma1[j] = mu11[j];
        }
        //update  gamma3
        for(int j = 0; j < p; j++){
            v32[j] = 1. / (b42 * rho1 * s3rs3_temp[j] + invsg32
                  + rho2 * s2rs2_temp[j] - 2 * rho_2 * rho2 * beta1 * s2rs4_temp[j] + rho2 * b12 * s4rs4_temp[j] ) ;     
            double mu3;
            mu3 = rho1 * (- b42 * (dot(s3rs3.row(j), gamma3) - gamma3[j] * s3rs3_temp[j] )
                  - rho_1 * beta4 * s1s3g1_hat[j] + beta4 * s3s3G1_hat[j]
                  + rho_1 * beta4 * dot(s1rs3.row(j), gamma1) - beta1 * beta4 * dot(s3rs3.row(j), gamma1) - beta4 * dot(s3rs3.row(j), beta2)
                 ) * v32[j]
            +rho2 * (- (dot(s2rs2.row(j), gamma3) - gamma3[j] * s2rs2_temp[j] )
                   + 2 * rho_2 * beta1 * (dot(s2rs4.row(j), gamma3) - gamma3[j] * s2rs4_temp[j])
                   - b12 * (dot(s4rs4.row(j), gamma3) - gamma3[j] * s4rs4_temp[j] )
                   - rho_2 * beta1 * s2s4g3_hat[j] + beta1 * s4s4G3_hat[j]
                   + s2s2g3_hat[j] - rho_2 * s2s4G3_hat[j]
                   + rho_2 * beta4 * dot(s2rs4.row(j), gamma1) - beta1 * beta4 * dot(s4rs4.row(j), gamma1)
                 ) * v32[j];
            mu33[j] = mu3 + randn()*sqrt(v32[j]);
            gamma3[j] = mu33[j];
        }
        //update beta2
        for(int j = 0; j < p; j++){
            vb22[j] = 1. / (rho1 * s3rs3_temp[j] + invsg22) ;     
            double mub2;
            mub2 = rho1 * (- (dot(s3rs3.row(j), beta2) - beta2[j] * s3rs3_temp[j] )
                   - rho_1 * s1s3g1_hat[j] + s3s3G1_hat[j]
                   + rho_1 * dot(s1rs3.row(j), gamma1) - beta1 * dot(s3rs3.row(j), gamma1) - beta4 * dot(s3rs3.row(j), gamma3)
                  ) * vb22[j];
            mub22[j] = mub2 + randn()*sqrt(vb22[j]);
            beta2[j] = mub22[j];
        }
        //update beta1
        vb12 = 1. / as_scalar( rho1 * gamma1.t() * s3rs3 * gamma1 + rho2 * gamma3.t() * s4rs4 * gamma3);
        double mub1;
        mub1 = rho1 * as_scalar(-rho_1 * gamma1.t() * s1s3g1_hat + gamma1.t() * s3s3G1_hat
                            + rho_1 * gamma1.t() * s1rs3 * gamma1 - gamma1.t() * s3rs3 * beta2 - beta4 * gamma1.t() * s3rs3 * gamma3
                           )* vb12
               +rho2 * as_scalar(-rho_2 * gamma3.t() * s2s4g3_hat + gamma3.t() * s4s4G3_hat
                            + rho_2 * gamma3.t() * s2rs4 * gamma3 - beta4 * gamma1.t() * s4rs4 * gamma3
                           )* vb12;
        mub11 = mub1 + randn()*sqrt(vb12);
        beta1 = mub11;  
        //update beta4
        vb42 = 1. / as_scalar( rho1 * gamma3.t() * s3rs3 * gamma3 +  rho2 * gamma1.t() * s4rs4 * gamma1);
        double mub4;
        mub4 = rho1 * as_scalar(-rho_1* gamma3.t() * s1s3g1_hat + gamma3.t() * s3s3G1_hat
                            + rho_1 * gamma3.t() * s1rs3 * gamma1 - beta1 * gamma1.t() * s3rs3 * gamma3 - beta2.t() * s3rs3 * gamma3
                           ) * vb42
               +rho2 * as_scalar(-rho_2* gamma1.t() * s2s4g3_hat + gamma1.t() * s4s4G3_hat
                            + rho_2 * gamma1.t() * s2rs4 * gamma3 - beta1 * gamma1.t() * s4rs4 * gamma3 
                           ) * vb42 ;
        mub44 = mub4 + randn()*sqrt(vb42);
        beta4 = mub44;
        //update sigma123
        sig12 = Updatesig2(gamma1, ak1, bk1, p );
        sigma12 = sig12;
        sig22 = Updatesig2(beta2, ak2, bk2, p );
        sigma22 = sig22;
        sig32 = Updatesig2(gamma3, ak3, bk3, p );
        sigma32 = sig32;
        //prepare for next iter
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
    ObjGibbs_Gam3seo obj;
    obj.gamma1 = gamma1;
    obj.beta2 = beta2;
    obj.gamma3 = gamma3;
    obj.Beta1res = Beta1res;
    obj.Beta4res = Beta4res;
    obj.Sg12Res = Sg12Res;
    obj.Sg22Res = Sg22Res;
    obj.Sg32Res = Sg32Res;
    return obj;
}



//[[Rcpp::export]]
List MRGEI_Gam3seo(arma::vec &gammah1,arma::vec &gammah3, arma::vec &Gammah1, arma::vec &Gammah3,
                   arma::vec &se1, arma::vec &se2, arma::vec &se3, arma::vec &se4, arma::mat &R, double &rho_1, double &rho_2)
{
    //mat R_sh = cal_blockcor(R) ;
    
    mat invS1 = inv(diagmat(se1));
    mat invS2 = inv(diagmat(se2));
    mat invS3 = inv(diagmat(se3));
    mat invS4 = inv(diagmat(se4));
    
    ObjGibbs_Gam3seo obj = MRGibbs_Gam3seo( gammah1, gammah3, Gammah1, Gammah3, invS1, invS2, invS3, invS4, R, rho_1, rho_2);

    double bhat1 = mean(obj.Beta1res);
    double b1se1 = stddev(obj.Beta1res);
    double pvalue1 = 2 * (R::pnorm(abs(bhat1 / b1se1),0,1,0,0));

    double bhat4 = mean(obj.Beta4res);
    double b4se4 = stddev(obj.Beta4res);
    double pvalue4 = 2 * (R::pnorm(abs(bhat4 / b4se4),0,1,0,0));

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
        Rcpp::Named("Sg32Res") = Rcpp::wrap(obj.Sg32Res)
        );
  return output;
}

