#include "rtsheader.h"

using namespace Rcpp;

// [[Rcpp::export]]
SEXP nngpCovariance__new(SEXP formula_, SEXP data_, SEXP colnames_,
                         SEXP T_,SEXP m_, SEXP gptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  XPtr<rts::griddata> gptr(gptr_);
  XPtr<rts::nngpCovariance> ptr(new rts::nngpCovariance(formula,data,colnames,T,m,*gptr));
  return ptr;
}

// [[Rcpp::export]]
SEXP GridData__new(SEXP x_, SEXP t_){
  Eigen::ArrayXXd x = as<Eigen::ArrayXXd>(x_);
  int t = as<int>(t_);
  XPtr<rts::griddata> ptr(new rts::griddata(x,t), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP GridData_nn__new(SEXP x_, SEXP t_, SEXP m_){
  Eigen::ArrayXXd x = as<Eigen::ArrayXXd>(x_);
  int t = as<int>(t_);
  int m = as<int>(m_);
  XPtr<rts::griddata> ptr(new rts::griddata(x,t,m), true);
  return ptr;
}

// [[Rcpp::export]]
SEXP RegionData__new(SEXP n_cell_, SEXP cell_id_, SEXP q_weights_, SEXP N_, SEXP T_){
  Eigen::ArrayXi n_cell = as<Eigen::ArrayXi>(n_cell_);
  Eigen::ArrayXi cell_id = as<Eigen::ArrayXi>(cell_id_);
  Eigen::ArrayXd q_weights = as<Eigen::ArrayXd>(q_weights_);
  int N = as<int>(N_);
  int T = as<int>(T_);
  XPtr<rts::RegionData> rptr(new rts::RegionData(n_cell,cell_id,q_weights,N,T));
  return rptr;
}

// [[Rcpp::export]]
SEXP Model_ar_lp__new(SEXP formula_, 
                      SEXP data_, 
                      SEXP grid_data_, 
                      SEXP colnames_,
                      SEXP beta_,
                      SEXP theta_, 
                      int T){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  XPtr<ModelAR>ptr(new ModelAR(formula,data,grid_data,colnames,T),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_nngp_lp__new(SEXP formula_, 
                        SEXP data_,
                        SEXP grid_data_, 
                        SEXP colnames_,
                        SEXP beta_,
                        SEXP theta_, 
                        int T, 
                        int m, 
                        SEXP gptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  XPtr<rts::griddata> gptr(gptr_);
  XPtr<ModelNNGP> ptr(new ModelNNGP(formula,data,grid_data,colnames,T,m,*gptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_hsgp_lp__new(SEXP formula_, 
                        SEXP data_,
                        SEXP grid_data_, 
                        SEXP colnames_,
                        SEXP beta_, 
                        SEXP theta_, 
                        int T, 
                        int m, 
                        SEXP L_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  Eigen::ArrayXd L = as<Eigen::ArrayXd>(L_);
  XPtr<ModelHSGP> ptr(new ModelHSGP(formula,data,grid_data,colnames,T,m,L),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_ar_region_grid__new(SEXP formula_region_, 
                               SEXP formula_grid_,
                               SEXP data_region_, 
                               SEXP data_grid_, 
                               SEXP colnames_region_, 
                               SEXP colnames_grid_,
                               SEXP beta_,
                               SEXP theta_, 
                               SEXP rptr_,
                               int T){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelARRegionG> ptr(new ModelARRegionG(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,T,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_nngp_region_grid__new(SEXP formula_region_, 
                                 SEXP formula_grid_,
                                 SEXP data_region_, 
                                 SEXP data_grid_, 
                                 SEXP colnames_region_, 
                                 SEXP colnames_grid_,
                                 SEXP beta_,
                                 SEXP theta_, SEXP rptr_,
                                 SEXP gptr_,
                                 int T,
                                 int m){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<rts::griddata> gptr(gptr_);
  XPtr<ModelNNGPRegionG> ptr(new ModelNNGPRegionG(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,*rptr,*gptr,T,m),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_hsgp_region_grid__new(SEXP formula_region_, 
                                 SEXP formula_grid_,
                                 SEXP data_region_, 
                                 SEXP data_grid_, 
                                 SEXP colnames_region_, 
                                 SEXP colnames_grid_,
                                 SEXP beta_,
                                 SEXP theta_, SEXP rptr_,
                                 int T,
                                 int m,
                                 SEXP L_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  Eigen::ArrayXd L = as<Eigen::ArrayXd>(L_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelHSGPRegionG> ptr(new ModelHSGPRegionG(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,T,m,L,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_ar_region__new(SEXP formula_, 
                          SEXP data_, 
                          SEXP grid_data_, 
                          SEXP colnames_,
                          SEXP beta_,
                          SEXP theta_, 
                          int T, 
                          SEXP rptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelARRegion> ptr(new ModelARRegion(formula,data,grid_data,colnames,T,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_nngp_region__new(SEXP formula_, 
                            SEXP data_,
                            SEXP grid_data_, 
                            SEXP colnames_,
                            SEXP beta_,
                            SEXP theta_, 
                            int T,
                            int m,
                            SEXP rptr_, 
                            SEXP gptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<rts::griddata> gptr(gptr_);
  XPtr<ModelNNGPRegion> ptr(new ModelNNGPRegion(formula,data,grid_data,colnames,T,m,*rptr,*gptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_hsgp_region__new(SEXP formula_, 
                            SEXP data_,
                            SEXP grid_data_, 
                            SEXP colnames_,
                            SEXP beta_,
                            SEXP theta_, 
                            int T,
                            int m,
                            SEXP rptr_, 
                            SEXP L_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  Eigen::ArrayXd L = as<Eigen::ArrayXd>(L_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelHSGPRegion> ptr(new ModelHSGPRegion(formula,data,grid_data,colnames,T,m,L,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}


