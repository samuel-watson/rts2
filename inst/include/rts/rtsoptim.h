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
    model& M_;
    int trace_;
    std::vector<double> start_b;
    std::vector<double> rho;
    
    rtsoptim(model& M,
             int trace): M_(M), trace_(trace), start_b(1,0.0), rho(1,0.0) {}
    
    double rho_optim(double rhostart){
      rts::rho_likelihood<model> ldl(M_);
      Rbobyqa<rts::rho_likelihood<model>,std::vector<double> > opt;
      opt.control.iprint = trace_;
      start_b[0] = rhostart;
      opt.minimize(ldl, start_b);
      rho = opt.par();
      return rho[0];
    }
};

}

#endif