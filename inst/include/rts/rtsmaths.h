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

//template<typename Derived>
inline MatrixXd sparse_matrix_mult(const sparse& A, const Eigen::Ref<const Eigen::MatrixXd>& X, bool transpose = false){
  if(!transpose){
    MatrixXd AX = MatrixXd::Zero(A.n,X.cols());
    double val;
    for(int i = 0; i < A.n; i++){
      for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
        val = A.Ax[j];
        for(int k = 0; k<X.cols(); k++){
          AX(i,k) += val*X(A.Ai[j],k);
        }
      }
    }
    return AX;
  } else {
    MatrixXd AX = MatrixXd::Zero(A.n,X.rows());
    double val;
    for(int i = 0; i < A.n; i++){
      for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
        val = A.Ax[j];
        for(int k = 0; k<X.rows(); k++){
          AX(i,k) += val*X(k,A.Ai[j]);
        }
      }
    }
    return AX.transpose();
  }
  
}

inline MatrixXd sparse_matrix_t_mult(const sparse& A, const MatrixXd& X, int T){
  MatrixXd AX = MatrixXd::Zero(A.n*T,X.cols());
  int N = X.rows()/T;
  for(int t = 0; t < T; t++){
    double val;
    for(int i = 0; i < A.n; i++){
      for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
        val = A.Ax[j];
        for(int k = 0; k<X.cols(); k++){
          AX(i+t*A.n,k) += val*X(A.Ai[j]+t*N,k);
        }
      }
    }
  }
  return AX;
}

inline VectorXd sparse_vector_t_mult(const sparse& A, const VectorXd& X, int T){
  VectorXd AX = VectorXd::Zero(A.n*T);
  int N = X.size()/T;
  for(int t = 0; t < T; t++){
    for(int i = 0; i < A.n; i++){
      for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
        AX(i+t*A.n) += A.Ax[j]*X(A.Ai[j]+t*N);
      }
    }
  }
  return AX;
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

inline MatrixXd cholesky(const MatrixXd& A){
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
  MatrixXd B = Map<MatrixXd>(L.data(), n, n);
  B = B.transpose();
  return B;
}

}