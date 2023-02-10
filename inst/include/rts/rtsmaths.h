#ifndef RTSMATHS_H
#define RTSMATHS_H

#define _USE_MATH_DEFINES

#include <cmath> 
#include <glmmr.h>
#include <RcppEigen.h>
#include "eigenext.h"
#ifdef _OPENMP
#include <omp.h>     
#else
// for machines with compilers void of openmp support
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#define omp_set_nested(a)   // empty statement to remove the call
#define omp_get_wtime()        0
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

namespace rts{

inline Eigen::MatrixXd inv_ldlt_AD(const Eigen::MatrixXd &A, 
                         const Eigen::VectorXd &D,
                         const Eigen::ArrayXXi &NN){
  int n = A.cols();
  int m = A.rows();
  Eigen::MatrixXd y = Eigen::MatrixXd::Zero(n,n);
#pragma omp parallel for  
  for(int k=0; k<n; k++){
    int idxlim;
    for (int i = 0; i < n; i++) {
      idxlim = i<=m ? i-1 : m;
      double lsum = 0;
      for (int j = 0; j < idxlim; j++) {
        lsum += -1.0 * A(j,i) * y(NN(j,i),k);
      }
      y(i,k) = i==k ? (1-lsum)  : (-1.0*lsum);
    }
  }
  
  return y*D.cwiseSqrt().asDiagonal();
}

}

#endif