#ifndef RTSOPTIM_H
#define RTSOPTIM_H

#define _USE_MATH_DEFINES

#include <cmath> 
#include <glmmr.h>
#include <glmmrMCML.h>
#include <RcppEigen.h>
#include "rtslikelihood.h"

// [[Rcpp::depends(RcppEigen)]]

namespace rts {

template<typename model>
class rtsoptim {
  public:
    model* M_;
    int trace_;
    
    rtsoptim(model* M,
             int trace): M_(M), trace_(trace) {}
    
    double rho_optim(double rhostart){
      rts::rho_likelihood<model> ldl(M_);
      Rbobyqa<rts::rho_likelihood<model>,std::vector<double> > opt;
      opt.control.iprint = trace_;
      std::vector<double> start_b;
      start_b.push_back(rhostart);
      opt.minimize(ldl, start_b);
      std::vector<double> rho = opt.par();
      return rho[0];
    }
};

}

#endif