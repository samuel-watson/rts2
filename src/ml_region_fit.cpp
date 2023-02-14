#include "../inst/include/rtsheader.h"
using namespace Rcpp;


//' Maximum likelihood optimisation step for MCMCML algorithm for the region LGCP model
//' 
//' Maximum likelihood optimisation step for linear predictor parameters and covariance 
//' parameters for MCMCML algorithm. The model is region data LGCP model. This function is 
//' used internally, use `lgcp_fit_ml()` function in the \link[rts2]{grid} class.
//' @param cov An integer matrix with columns of block identifier, dimension of block, function definition, number of variables
//' in the argument to the funciton, and index of the parameters, respectively. Rows are specific functions of each block.
//' @param data Vector of data. Created by flattening the matrices in column-major order of the data used in each block.
//' @param X Matrix of covariates for the linear predictor
//' @param y Vector of outcome counts for the computational grid cells
//' @param u Matrix of samples of the random effects terms with rows equal to the number of grid cells times by number of time 
//' periods, and columns equal to the number of samples.
//' @param nT Integer specifying the number of time periods
//' @param start Vector of starting values for the model parameters in order beta, theta (covariance parameters), rho
//' @param offset Vector of offset values
//' @param n_cell Vector of integers specifying the number of grid cells intersecting each region, 
//' see the `region_data()` function in the \link[rts2]{grid} class.
//' @param cell_id Vector of indices of grid cells overlapping each region, in order,
//' see the `region_data()` function in the \link[rts2]{grid} class.
//' @param q_weights Vector specifying the proportion of the area of each region covered by overlapping 
//' grid cells, see the `region_data()` function in the \link[rts2]{grid} class.
//' @param trace Either 0, 1, or 2 indicating the level of detail to print to the console
//' @param mcnr Logical indicating whether to use newton-raphson (`TRUE`) or expectation maximisation (`FALSE`) for the estimation 
//' of the parameters of the linear predictor
//' @param known_cov Logical indicating whether to treat the covariance parameters as known.
//' @return A named list with beta, theta, and rho parameters
// [[Rcpp::export]]
Rcpp::List lgcp_region_optim(const Eigen::ArrayXXi &cov,
                      const Eigen::ArrayXd &data,
                      const Eigen::MatrixXd &X,
                      const Eigen::VectorXd &y, 
                      Eigen::MatrixXd u,
                      int nT,
                      Eigen::ArrayXd start,
                      const Eigen::VectorXd &offset,
                      const Eigen::ArrayXi &n_cell,
                      const Eigen::ArrayXi &cell_id,
                      const Eigen::ArrayXd &q_weights,
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
  
  rts::rtsDMatrix dmat(cov,data,eff_range, thetapars, nT);
  Eigen::VectorXd beta = start.segment(0,X.cols());
  rts::rtsRegionModel model(u,X,y,beta,offset,nT,n_cell,cell_id,q_weights);
  int startlen = start.size();
  if(nT>1) startlen--;
  glmmr::mcmloptim<rts::rtsDMatrix, rts::rtsRegionModel> mc(dmat,model, start.segment(0,startlen),trace);
  rts::rtsoptim<rts::rtsRegionModel> rmc(model, trace);
  
  if(!mcnr){
    mc.l_optim();
  } else {
    mc.mcnr();
  }
  if(nT>1)rho = rmc.rho_optim(rho);
  if(!known_cov)mc.d_optim();
  
  beta = mc.get_beta();
  if(!known_cov)thetapars = mc.get_theta();
  
  Rcpp::List out = Rcpp::List::create(_["beta"] = beta, _["theta"] = thetapars,
                                      _["rho"] = rho);
  return out;
}