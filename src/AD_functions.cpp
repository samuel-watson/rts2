#include "../inst/include/rtsheader.h"
using namespace Rcpp;

//' Generate the matrices A and D for the nearest neighbour approximation
//' 
//' @details
//' The neatest neighbour Gaussian process approximation to the inverse covariance matrix
//' is \eqn{(I-A)D^-1(I-A)^T}. This function generates matrix D and a "collapsed" version of matrix A.
//' If there are `m` nearest neighbours then A will have `m` rows and columns equal to the number of observations, `n`.
//' The "full" A matrix is n x n. See \link[rts2]{full_A} for converting A to full representation.
//' @param cov An integer matrix with columns of block identifier, dimension of block, function definition, number of variables
//' in the argument to the funciton, and index of the parameters, respectively. Rows are specific functions of each block.
//' @param data Vector of data. Created by flattening the matrices in column-major order of the data used in each block.
//' @param NN A matrix of integers with rows equal to the number of nearest neighbours and columns equal to the number of observations.
//' Each element is the set of nearest neighbours. See \link[rts2]{genNN}.
//' @param gamma Vector of parameters used to generate the matrix D. 
//' @return A named list containing matrices A and D.
// [[Rcpp::export]]
Rcpp::List get_AD(const Eigen::ArrayXXi &cov,
                  const Eigen::ArrayXd &data,
                  const Eigen::ArrayXXi &NN,
                  Eigen::ArrayXd theta){
  
  Eigen::ArrayXd eff_range = Eigen::ArrayXd::Zero(10);
  rts::NNGPDmatrix dmat(cov,data,eff_range, NN, theta,1);
  
  Rcpp::List out = Rcpp::List::create(_["A"] = dmat.A_, _["D"] = dmat.D_);
  return out;
}

//' Coverts a collapsed matrix A to full matrix A
//' 
//' The output of the function \link{rts2}[get_AD] is a m x n matrix A, where 
//' m is the number of nearest neighbours. This function coverts A into a full 
//' n x n matrix.
//' @param A A "collapsed" matrix A output from \link[rts2]{get_AD}
//' @param NN A matrix of integers with rows equal to the number of nearest neighbours and columns equal to the number of observations.
//' Each element is the set of nearest neighbours. See \link[rts2]{genNN}.
//' @return A matrix
 // [[Rcpp::export]]
Eigen::MatrixXd full_A(const Eigen::MatrixXd &A,
                       const Eigen::ArrayXXi &NN){
  int m = A.rows();
  int n = A.cols();
  int idxlim;
  Eigen::MatrixXd Af = Eigen::MatrixXd::Zero(n,n);
  for(int i = 1; i<n; i++){
    idxlim = i<=m ? i-1 : m;
    for(int j=0; j<idxlim; j++){
      Af(i,NN(j,i)) = A(j,i);
    }
  }
  return Af;
}

//' Generates an approximate Cholesky decomposition of the covariance matrix
//' 
//' Takes the output from the function \link[rts2]{get_AD} and generates the approximate
//' Cholesky decomposition of the covariance matrix of the random effects.
//' 
//' @details
//' For the nearest neighbour Gaussian process approximation, the covariance matrix of the random effects/
//' latent fieldis approximated by \eqn{\Sigma = (I-A)^{-T}D(I-A)^{-1}} where matrices 
//' A and D are determined by the set of nearest neighbours in `NN` (see \link[rts2]{genNN}). 
//' The matrices A and D are returned by the function \link[rts2]{get_AD}.
//' Many of the calculations for the maximum likelihood model fitting methods 
//' require the Cholesky decomposition of the covariance matrix \eqn{\Sigma = LL^T}, which is approximated here 
//' by \eqn{(I-A)^{-T}D^{1/2}}. This function generates this lower triangular matrix
//' using a forward substitution algorithm as A is strictly lower triangular and D is diagonal.
//' @param A A "collapsed" matrix A output from \link[rts2]{get_AD}
//' @param D A diagonal matrix for the approximate LDL decomposition from \link[rts2]{get_AD}. 
//' @param NN A matrix of integers with rows equal to the number of nearest neighbours and columns equal to the number of observations.
//' Each element is the set of nearest neighbours. See \link[rts2]{genNN}.
//' @return A lower triangular matrix  
// [[Rcpp::export]]
Eigen::MatrixXd inv_ldlt(const Eigen::MatrixXd &A, 
                         const Eigen::VectorXd &D,
                         const Eigen::ArrayXXi &NN){
  return rts::inv_ldlt_AD(A,D,NN);
}

//' Generates the matrix ZL
//' 
//' For the MCMCML algorithms the matrix ZL multiplies the random effects u~N(0,I) to produce
//' the latent field.
//' 
//' @details
//' The LGCP model specified in this package uses an auto-regressive specification for the random effect 
//' terms. If S(t) are the random effects for time period t, then these are specified as S(t) = rho*S(t-1) + v(t)
//' where v(t) is a spatial innovation term specified as a spatial Gaussian process. The Cholesky decomposition 
//' of the covariance of v(t) is L and is the same in each time period, such that Lu = v where u~N(0,I). The matrix
//' Z generates the weighted sums of v(t) such that S = ZLu. This function produces ZL.
//' @param L Lower triangular matrix representing the Cholesky decomposition of the spatial terms in one time period.
//' @param nT Integer specifying the number of time periods
//' @param rho Scalar, the auto-regressive parameter value.
//' @return A matrix ZL
// [[Rcpp::export]]
Eigen::MatrixXd get_ZL(const Eigen::MatrixXd &L,
                       int nT,
                       double rho){
  int n = L.rows();
  Eigen::MatrixXd ZL(n*nT,n*nT);
  if(nT==1){
    ZL = L;
  } else {
    for(int t=0; t<nT;t++){
      for(int s=t; s<nT;s++){
        if(t==0){
          ZL.block(s*n,t*n,n,n) = (pow(rho,s)/(1-rho*rho))*L;
        } else {
          ZL.block(s*n,t*n,n,n) = pow(rho,s-t)*L;
        }
      }
      
    }
  }
  return ZL;
}

//' Generates the random effect terms from independent samples
//' 
//' For the MCMCML algorithms the spatial innovation terms are v. The Cholesky decomposition of the covariance matrix
//' of v is L. Random effects are sampled as standard normals u~N(0,I). This function generates v from samples u.
//' 
//' @details
//' The LGCP model specified in this package uses an auto-regressive specification for the random effect 
//' terms. If S(t) are the random effects for time period t, then these are specified as S(t) = rho*S(t-1) + v(t)
//' where v(t) is a spatial innovation term specified as a spatial Gaussian process. The Cholesky decomposition 
//' of the covariance of v(t) is L and is the same in each time period, such that Lu = v where u~N(0,I). This function 
//' produces v from L and u. The MCMCML algorithms also generate multiple samples of u so multiple samples of v are returned.
//' @param L Lower triangular matrix representing the Cholesky decomposition of the spatial terms in one time period.
//' @param u A matrix of samples of u
//' @return A matrix Lu.
// [[Rcpp::export]]
Eigen::MatrixXd get_Lu(const Eigen::MatrixXd &L,
                       const Eigen::MatrixXd &u){
  int nT = u.rows()/L.rows();
  int n = L.rows();
  if(nT==1){
    return L*u;
  } else {
    Eigen::MatrixXd Lu(u.rows(),u.cols());
    for(int t=0; t<nT; t++){
      Lu.block(t*n,0,n,u.cols()).noalias() = L*u.block(t*n,0,n,u.cols());
    }
    return Lu;
  }
}

//' Generates samples of the latent field from samples of u
//' 
//' For the MCMCML algorithms the matrice Z and L multiply the random effects u~N(0,I) to produce
//' the latent field.
//' 
//' @details
//' The LGCP model specified in this package uses an auto-regressive specification for the random effect 
//' terms. If S(t) are the random effects for time period t, then these are specified as S(t) = rho*S(t-1) + v(t)
//' where v(t) is a spatial innovation term specified as a spatial Gaussian process. The Cholesky decomposition 
//' of the covariance of v(t) is L and is the same in each time period, such that Lu = v where u~N(0,I). The matrix
//' Z generates the weighted sums of v(t) such that S = ZLu. This function produces S given Z, L, and u.
//' @param L Lower triangular matrix representing the Cholesky decomposition of the spatial terms in one time period.
//' @param u A matrix of samples of u
//' @param nT Integer specifying the number of time periods
//' @param rho Scalar, the auto-regressive parameter value.
//' @return A matrix ZL
// [[Rcpp::export]]
Eigen::MatrixXd get_ZLu(const Eigen::MatrixXd &L,
                        const Eigen::MatrixXd &u,
                        int nT,
                        double rho){
  
  int niter = u.cols();
  int n = L.rows();
  Eigen::MatrixXd ZLu(n*nT,niter);
  
  if(nT == 1){
    ZLu = L * u;
  } else {
    for(int t=0; t<nT;t++){
      for(int s=t; s<nT;s++){
        if(t==0){
          ZLu.block(s*n,0,n,niter) = (pow(rho,s)/(1-rho*rho))*L*u.block(s*n,0,n,niter);
        } else {
          ZLu.block(s*n,0,n,niter) += pow(rho,s-t)*L*u.block(s*n,0,n,niter);
        }
      }
    }
  }
  
  return ZLu;
}

//' Generates predictions of the intensity of each region for the region model
//' 
//' The "region model" specifies the intensity for each region as a combination of the
//' offset, linear predictor, and a weighted combination of the intensity of overlapping 
//' computational grid cells. This function returns the region intensities given above 
//' components.
//' @param xb A vector of values of the linear predictor X*beta
//' @param y Vector of outcome counts 
//' @param zu Matrix with samples of the random field, see \link[rts2]{get_ZLu}.
//' @param offset Vector of offset values
//' @param nT Integer. Number of time periods.
//' @param n_cell Vector of integers specifying the number of grid cells intersecting each region, 
//' see the `region_data()` function in the \link[rts2]{grid} class.
//' @param cell_id Vector of indices of grid cells overlapping each region, in order,
//' see the `region_data()` function in the \link[rts2]{grid} class.
//' @param q_weights Vector specifying the proportion of the area of each region covered by overlapping 
//' grid cells, see the `region_data()` function in the \link[rts2]{grid} class.
//' @return A matrix with samples of the intensity in each region at each time point.
// [[Rcpp::export]]
Eigen::ArrayXXd region_intensity(const Eigen::ArrayXd &xb,
                                 const Eigen::ArrayXd &y,
                                 const Eigen::ArrayXXd &zu,
                                 const Eigen::ArrayXd &offset,
                                 int nT,
                                 const Eigen::ArrayXi &n_cell,
                                 const Eigen::ArrayXi &cell_id,
                                 const Eigen::ArrayXd &q_weights){
  Eigen::ArrayXXd intens = Eigen::ArrayXXd::Zero(y.size(),zu.cols());
  for(int i = 0; i< zu.cols();i++)intens.col(i) = xb.exp();
  int nRegion_ = y.size()/nT;
  int nCell_ = zu.rows()/nT;
  
#pragma omp parallel for
  for(int r=0; r<nRegion_;r++){
    for(int t=0; t<nT; t++){
      int nInter = n_cell(r+1)-n_cell(r);
      for(int j=0; j<zu.cols(); j++){
        double accum = 0;
        for(int l=0; l<nInter; l++){
          accum += q_weights(n_cell(r)-1+l)*exp(zu(cell_id(n_cell(r)-1+l) + t*nCell_,j));
        }
        intens(r + t*nRegion_,j) *= accum;
      }
    }
  }
  return intens;
}

// // [[Rcpp::export]]
// Eigen::MatrixXd sample_u(const Eigen::ArrayXXi &cov,
//                           const Eigen::ArrayXd &data,
//                           const Eigen::ArrayXd &eff_range,
//                           const Eigen::VectorXd& gamma,
//                           int nT,
//                           double rho,
//                           int nsample){
//   glmmr::DData dat(cov,data,eff_range);
//   glmmr::DMatrix dmat(&dat, gamma);
//   Eigen::MatrixXd L = dmat.genD(0,true,false);
//   int n = L.rows();
//   Eigen::MatrixXd samps(n*nT,nsample);
//   for(int t=0; t<nT; t++){
//     for(int i =0; i<nsample; i++){
//       samps.col(i).segment(t*n,n) = dmat.sim_re();
//     }
//   }
//   Eigen::MatrixXd ZLu = get_ZLu(L,samps,nT,n,rho);
//   return ZLu;
// }