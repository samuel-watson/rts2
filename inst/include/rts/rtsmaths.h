#ifndef RTSMATHS_H
#define RTSMATHS_H

#define _USE_MATH_DEFINES

#include <cmath> 
#include <glmmr/openmpheader.h>
#include <RcppEigen.h>
#include <queue>
#include <glmmr/algo.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

namespace rts{

using namespace Eigen;

inline MatrixXd inv_ldlt_AD(const MatrixXd &A, 
                         const VectorXd &D,
                         const ArrayXXi &NN){
  int n = A.cols();
  int m = A.rows();
  MatrixXd y = MatrixXd::Zero(n,n);
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

inline ArrayXi top_i_pq(ArrayXd v, int n) {
  typedef std::pair<double, int> Elt;
  std::priority_queue< Elt, std::vector<Elt>, std::greater<Elt> > pq;
  std::vector<int> result;
  
  for (int i = 0; i != v.size(); ++i) {
    if (pq.size() < n)
      pq.push(Elt(v[i], i));
    else {
      Elt elt = Elt(v[i], i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }
  
  //result.reserve(pq.size());
  ArrayXi res(pq.size());
  int iter = 0;
  while (!pq.empty()) {
    res(iter) = pq.top().second + 1;
    //result.push_back(pq.top().second + 1);
    pq.pop();
    iter++;
  }
  
  return res;
}

}

#endif