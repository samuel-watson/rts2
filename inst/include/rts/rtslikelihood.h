#ifndef RTSLIKELIHOOD_H
#define RTSLIKELIHOOD_H

#include <rbobyqa.h>
#include <RcppEigen.h>
#include "rtsmodel.h"

using namespace rminqa;


// [[Rcpp::depends(RcppEigen)]]


namespace rts{

template<typename model>
class rho_likelihood : public Functor<std::vector<double> > {
  model* M_;
public:
  rho_likelihood(model* M) :  
  M_(M){}
  double operator()(const std::vector<double> &par) {
    double parrho = par[0];
    M_->update_rho(parrho);
    double ll = M_->log_likelihood();
    return -1*ll;
  }
};

}


#endif