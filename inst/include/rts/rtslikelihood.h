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
  model& M_;
  double parrho;
  double ll;
public:
  rho_likelihood(model& M) :  
  M_(M), parrho(0.0), ll(0.0) {}
  double operator()(const std::vector<double> &par) {
    parrho = par[0];
    M_.update_rho(parrho);
    ll = M_.log_likelihood();
    return -1*ll;
  }
};

}


#endif