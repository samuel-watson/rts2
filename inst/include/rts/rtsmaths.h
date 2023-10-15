#pragma once

#include <glmmr.h> 

namespace rts{

using namespace Eigen;

//kroenecker product
inline MatrixXd kronecker(const MatrixXd& A, const MatrixXd& B){
  MatrixXd result = MatrixXd::Zero(A.rows()*B.rows(), A.cols()*B.cols());
#pragma omp parallel for collapse(2)
  for(int i = 0; i < A.rows(); i ++){
    for(int j = 0; j < A.cols(); j++){
      if(A(i,j)!=0) result.block(i*B.rows(),j*B.cols(),B.rows(),B.cols()) = A(i,j)*B;
    }
  }
  return result;
}

inline sparse ar_factor_inv_to_sparse(const MatrixXd& a, int n){
  int t = a.cols();
  sparse x(n*t,n*t);
  if(t == 1){
    x = identity(n);
  } else {
    x.Ap = intvec(n*t+1);
    for(int j = 0; j < t; j++){
      for(int i = 0; i < n; i++){
        x.Ap[i+j*n] = (i+j*n)*t;
        for(int k = 0; k < t; k++){
          x.Ai.push_back(i+n*k);
          x.Ax.push_back(a(j,k));
        }
      }
    }
    x.Ap[n*t] = x.Ai.size();
  }
  return x;
}

inline void cholesky(MatrixXd& B, const MatrixXd& A){
  int n = A.rows();
  std::vector<double> L(n * n, 0.0);
  
  for (int j = 0; j < n; j++) {
    double s = glmmr::algo::inner_sum(&L[j * n], &L[j * n], j);
    L[j * n + j] = sqrt(A(j,j) - s);
    for (int i = j + 1; i < n; i++) {
      double s = glmmr::algo::inner_sum(&L[j * n], &L[i * n], j);
      L[i * n + j] = (1.0 / L[j * n + j] * (A(j,i) - s));
    }
  }
  B = Map<MatrixXd>(L.data(), n, n);
  B = B.transpose();
}

}