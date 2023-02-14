#include "../inst/include/rtsheader.h"
using namespace Rcpp;

//' Maximum likelihood optimisation step for MCMCML algorithm using nearest neighbour approximation
//' 
//' Maximum likelihood optimisation step for linear predictor parameters and covariance 
//' parameters for MCMCML algorithm. The model is a point data LGCP using the nearest neighbour 
//' Gaussian process approximation. This function is  used internally, 
//' use `lgcp_fit_ml()` function in the \link[rts2]{grid} class.
//' @param cov An integer matrix with columns of block identifier, dimension of block, function definition, number of variables
//' in the argument to the funciton, and index of the parameters, respectively. Rows are specific functions of each block.
//' @param data Vector of data. Created by flattening the matrices in column-major order of the data used in each block.
//' @param X Matrix of covariates for the linear predictor
//' @param y Vector of outcome counts for the computational grid cells
//' @param NN A matrix of integers with rows equal to the number of nearest neighbours and columns equal to the number of observations.
//' Each element is the set of nearest neighbours. See \link[rts2]{genNN}.
//' @param u Matrix of samples of the random effects terms with rows equal to the number of grid cells times by number of time 
//' periods, and columns equal to the number of samples.
//' @param nT Integer specifying the number of time periods
//' @param start Vector of starting values for the model parameters in order beta, theta (covariance parameters), rho
//' @param offset Vector of offset values
//' @param trace Either 0, 1, or 2 indicating the level of detail to print to the console
//' @param mcnr Logical indicating whether to use newton-raphson (`TRUE`) or expectation maximisation (`FALSE`) for the estimation 
//' of the parameters of the linear predictor
//' @param known_cov Logical indicating whether to treat the covariance parameters as known.
//' @return A named list with beta, theta, and rho parameters
// [[Rcpp::export]]
Rcpp::List nngp_optim(const Eigen::ArrayXXi &cov,
                      const Eigen::ArrayXd &data,
                      const Eigen::MatrixXd &X,
                      const Eigen::VectorXd &y, 
                      const Eigen::ArrayXXi &NN,
                      Eigen::MatrixXd u,
                      int nT,
                      Eigen::ArrayXd start,
                      const Eigen::VectorXd &offset,
                      int trace,
                      bool mcnr = true,
                      bool known_cov = false){
  
  //check start is the right length
  int theta_expect = glmmr::algo::get_n_cov_pars(cov);
  int npars_expect = theta_expect + X.cols();
  if(nT>1) npars_expect++;
  if(start.size() != npars_expect){
    Rcpp::stop("Wrong number of starting parameter values.");
  }
  
  Eigen::ArrayXd eff_range = Eigen::ArrayXd::Zero(10);
  double rho = nT>1 ? start(start.size()-1) : 0.0;
  Eigen::ArrayXd thetapars = start.segment(X.cols(),theta_expect);
  Eigen::VectorXd beta = start.segment(0,X.cols());
  
  rts::NNGPDmatrix dmat(cov,data,eff_range, NN, thetapars, nT);
  rts::rtsModel model(u,X,y,beta,offset,nT);
  int startlen = start.size();
  if(nT>1) startlen--;
  glmmr::mcmloptim<rts::NNGPDmatrix, rts::rtsModel> mc(dmat,model, start.segment(0,startlen),trace);
  rts::rtsoptim<rts::rtsModel> rmc(model, trace);
  
  if(!mcnr){
    mc.l_optim();
  } else {
    mc.mcnr();
  }
  if(nT>1)rho = rmc.rho_optim(rho);
  if(!known_cov)mc.d_optim();
  
  beta = mc.get_beta();
  thetapars = mc.get_theta();
  
  Rcpp::List out = Rcpp::List::create(_["beta"] = beta, 
                                      _["theta"] = thetapars,
                                      _["rho"] = rho);
  return out;
}


//' Maximum likelihood model fitting using Laplace approximation to likelihood and nearest neighbour Gaussian process approximation
//' 
//' Model fitting of the full Log Gaussian Cox Process model using a Laplace approximation to the likelihood.
//' This function is called internally, see the `fit_lgcp_la` in the \link[rts2]{grid} class.
//' @param cov An integer matrix with columns of block identifier, dimension of block, function definition, number of variables
//' in the argument to the funciton, and index of the parameters, respectively. Rows are specific functions of each block.
//' @param data Vector of data. Created by flattening the matrices in column-major order of the data used in each block.
//' @param X Matrix of covariates for the linear predictor
//' @param y Vector of outcome counts for the computational grid cells
//' @param NN A matrix of integers with rows equal to the number of nearest neighbours and columns equal to the number of observations.
//' Each element is the set of nearest neighbours. See \link[rts2]{genNN}.
//' @param start Vector of starting values for the model parameters in order beta, theta (covariance parameters), rho
//' @param offset Vector of offset values
//' @param nT Integer specifying the number of time periods
//' @param known_cov Logical indicating whether to treat the covariance parameters as known.
//' @param usehess Logical indicating whether to use the Hessian matrix to estimate standard errors
//' @param tol Maximum difference between parameter values between iterations at which the algorithm is considered to have converged
//' @param verbose Logical indicating whether to provide output to the console
//' @param trace Either 0, 1, or 2 indicating the level of detail to print to the console
//' @param maxiter Integer. Maximum number of iterations of the algorithm.
//' @return A named list with beta, theta, and rho parameters with standard errors and random effect estimates
// [[Rcpp::export]]
Rcpp::List nngp_optim_la(const Eigen::ArrayXXi &cov,
                         const Eigen::ArrayXd &data,
                         const Eigen::MatrixXd &X,
                         const Eigen::VectorXd &y, 
                         const Eigen::ArrayXXi &NN,
                         Eigen::ArrayXd start,
                         const Eigen::VectorXd &offset,
                         int nT,
                         bool known_cov = false,
                         bool usehess = false,
                         double tol = 1e-3,
                         bool nr = false,
                         bool verbose = true,
                         int trace = 0,
                         int maxiter = 10){
  
  //check start is the right length
  int theta_expect = glmmr::algo::get_n_cov_pars(cov);
  int npars_expect = theta_expect + X.cols();
  if(nT>1) npars_expect++;
  if(start.size() != npars_expect){
    Rcpp::stop("Wrong number of starting parameter values.");
  }
  
  double rho = nT>1 ? start(start.size()-1) : 0.0;
  Eigen::ArrayXd eff_range = Eigen::ArrayXd::Zero(10);
  Eigen::VectorXd theta = start.segment(X.cols(),theta_expect);
  Eigen::VectorXd beta = start.segment(0,X.cols());
  
  rts::NNGPDmatrix dmat(cov,data,eff_range, NN, theta, nT);
  Eigen::MatrixXd L = dmat.chol();
  rts::rtsModel model(L,1,X,y,beta,offset,nT);
  int startlen = start.size();
  if(nT>1) startlen--;
  glmmr::mcmloptim<rts::NNGPDmatrix,rts::rtsModel> mc(dmat,model, start.segment(0,startlen),trace);
  rts::rtsoptim<rts::rtsModel> rmc(model, trace);
  
  Eigen::ArrayXd diff = Eigen::ArrayXd::Zero(start.size());
  double maxdiff = 1;
  int iter = 1;
  double newrho = rho;
  Eigen::VectorXd newbeta = Eigen::VectorXd::Zero(beta.size());
  Eigen::VectorXd newtheta = Eigen::VectorXd::Zero(theta.size());
  if(known_cov)newtheta = theta;
  bool converged = false;
  if(trace > 0 ) Rcpp::Rcout << "\n STARTING LA \n " ;
  
  while(maxdiff > tol && iter <= maxiter){
    if(verbose)Rcpp::Rcout << "\n\nIter " << iter << "\n" << std::string(40, '-');
    if(nr){
      mc.mcnr_b();
    } else {
      mc.la_optim();
    }
    newbeta = mc.get_beta();
    model.update_beta(newbeta);
    if(nT>1)newrho = rmc.rho_optim(rho);
    if(!known_cov){
      mc.la_optim_cov();
      newtheta = mc.get_theta();
    }
    
    diff.segment(0,beta.size()) = (beta - newbeta).cwiseAbs();
    diff.segment(beta.size(),theta.size()) = (theta - newtheta).cwiseAbs();
    if(nT>1)diff(diff.size()-1) = abs(rho- newrho);
    maxdiff = diff.maxCoeff();
    
    if(maxdiff < tol) converged = true;
    
    //update all the parameters
    beta = newbeta;
    if(!known_cov)theta = newtheta;
    rho = newrho;
    if(!converged){
      if(!known_cov){
        dmat.update_parameters(newtheta);
        dmat.genAD();
        model.update_L(dmat.chol());
      }
      model.update_rho(rho);
      model.update_beta(beta);
    }
    
    iter++;
    
    if(verbose){
      
      Rcpp::Rcout << "\nbeta: " << beta.transpose() << "\ntheta: " << theta.transpose();
      if(nT>1)Rcpp::Rcout << "\nrho: " << rho;
      Rcpp::Rcout << "\n Max. diff: " << maxdiff;
      if(converged)Rcpp::Rcout << " CONVERGED!";
      Rcpp::Rcout << "\n" << std::string(40, '-');
    }
  }
  
  if((trace > 0) & !converged & (iter > maxiter))Rcpp::Rcout << "\nWARNING: Reached maximum iterations and not converged.";
  
  if(!known_cov){
    mc.la_optim_bcov();
    beta = mc.get_beta();
    theta = mc.get_theta();
    
    if(verbose){
      Rcpp::Rcout << "\n" << std::string(40, '-') << "\nCompleted: \nbeta: " << beta.transpose() << "\ntheta: " << theta.transpose();
      Rcpp::Rcout << "\n" << std::string(40, '-');
    }
  }
  
  //estimate standard errors
  Eigen::VectorXd se = Eigen::VectorXd::Zero(start.size());
  
  if(usehess){
    Eigen::MatrixXd hess = mc.hess_la();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(hess.rows(),hess.cols());
    
    hess = hess.llt().solve(I);
    for(int i = 0; i < hess.cols();i++)se(i) = sqrt(hess(i,i));
  } 
  
  
  Rcpp::List res = Rcpp::List::create(_["beta"] = beta, _["theta"] = theta,
                                      _["rho"] = rho, _["se"] = se,
                                      _["u"] = model.u_);
  return res;
}