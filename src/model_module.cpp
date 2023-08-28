#include "rtsheader.h"

using namespace Rcpp;

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
SEXP RegionData__new(SEXP n_cell_, SEXP cell_id_, SEXP q_weights_, SEXP g_ptr_){
  Eigen::ArrayXi n_cell = as<Eigen::ArrayXi>(n_cell_);
  Eigen::ArrayXi cell_id = as<Eigen::ArrayXi>(cell_id_);
  Eigen::ArrayXd q_weights = as<Eigen::ArrayXd>(q_weights_);
  XPtr<rts::griddata> ptr(g_ptr_);
  XPtr<rts::RegionData> rptr(new rts::RegionData(n_cell,cell_id,q_weights,*ptr));
  return rptr;
}

// [[Rcpp::export]]
SEXP GridData__NN(SEXP ptr_){
  XPtr<rts::griddata> ptr(ptr_);
  Eigen::ArrayXXi nn = ptr->NN;
  return wrap(nn);
}

// [[Rcpp::export]]
void GridData__gen_NN(SEXP ptr_, SEXP m_){
  XPtr<rts::griddata> ptr(ptr_);
  int m = as<int>(m_);
  ptr->genNN(m);
}

// [[Rcpp::export]]
SEXP ModelBits_ar_lp__new(SEXP formula_, SEXP data_, SEXP colnames_,
                        SEXP family_, SEXP link_, SEXP beta_,
                        SEXP theta_, SEXP T_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::string family = as<std::string>(family_);
  std::string link = as<std::string>(link_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  XPtr<BitsAR>ptr(new BitsAR(formula,data,colnames,family,link,T),true);
  ptr->linear_predictor.update_parameters(beta);
  ptr->covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP ModelBits_nngp_lp__new(SEXP formula_, SEXP data_, SEXP colnames_,
                          SEXP family_, SEXP link_, SEXP beta_,
                          SEXP theta_, SEXP T_, SEXP m_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::string family = as<std::string>(family_);
  std::string link = as<std::string>(link_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  XPtr<BitsNNGP>ptr(new BitsNNGP(formula,data,colnames,family,link,T,m),true);
  ptr->linear_predictor.update_parameters(beta);
  ptr->covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP ModelBits_ar_region__new(SEXP formula_region_, SEXP formula_grid_,
                              SEXP data_region_, SEXP data_grid_, 
                              SEXP colnames_region_, SEXP colnames_grid_,
                              SEXP family_, SEXP link_, 
                              SEXP beta_,
                              SEXP theta_, SEXP rptr_,
                              SEXP T_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::string family = as<std::string>(family_);
  std::string link = as<std::string>(link_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<BitsARRegion> ptr(new BitsARRegion(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,family,link,T,*rptr),true);
  ptr->linear_predictor.update_parameters(beta);
  ptr->covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP ModelBits_nngp_region__new(SEXP formula_region_, SEXP formula_grid_,
                              SEXP data_region_, SEXP data_grid_, 
                              SEXP colnames_region_, SEXP colnames_grid_,
                              SEXP family_, SEXP link_, 
                              SEXP beta_,
                              SEXP theta_, SEXP rptr_,
                              SEXP T_, SEXP m_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::string family = as<std::string>(family_);
  std::string link = as<std::string>(link_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<BitsNNGPRegion>ptr(new BitsNNGPRegion(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,family,link,*rptr,T,m),true);
  ptr->linear_predictor.update_parameters(beta);
  ptr->covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
void rtsModelBits__update_beta(SEXP xp, SEXP beta_, SEXP covtype_, SEXP lptype_){
  // covtype 1 = ar, 2 = nn, lptype 1 = linearpredictor, 2 = regionpredictor
  std::vector<double> beta = as<std::vector<double> >(beta_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<BitsAR> ptr(xp);
    ptr->linear_predictor.update_parameters(beta);
  } else if(covtype == 2 && lptype == 1){
    XPtr<BitsNNGP> ptr(xp);
    ptr->linear_predictor.update_parameters(beta);
  } else if(covtype == 1 && lptype == 2){
    XPtr<BitsARRegion> ptr(xp);
    ptr->linear_predictor.update_parameters(beta);
  } else if(covtype == 2 && lptype == 2){
    XPtr<BitsNNGPRegion> ptr(xp);
    ptr->linear_predictor.update_parameters(beta);
  }
  
}

// [[Rcpp::export]]
void rtsModelBits__update_theta(SEXP xp, SEXP theta_, SEXP covtype_, SEXP lptype_){
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<BitsAR> ptr(xp);
    ptr->covariance.update_parameters(theta);
  } else if(covtype == 2 && lptype == 1){
    XPtr<BitsNNGP> ptr(xp);
    ptr->covariance.update_parameters(theta);
  } else if(covtype == 1 && lptype == 2){
    XPtr<BitsARRegion> ptr(xp);
    ptr->covariance.update_parameters(theta);
  } else if(covtype == 2 && lptype == 2){
    XPtr<BitsNNGPRegion> ptr(xp);
    ptr->covariance.update_parameters(theta);
  }
  
}

// [[Rcpp::export]]
void rtsModelBits__update_rho(SEXP xp, SEXP rho_, SEXP covtype_, SEXP lptype_){
  double rho = as<double>(rho_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<BitsAR> ptr(xp);
    ptr->covariance.update_rho(rho);
  } else if(covtype == 2 && lptype == 1){
    XPtr<BitsNNGP> ptr(xp);
    ptr->covariance.update_rho(rho);
  } else if(covtype == 1 && lptype == 2){
    XPtr<BitsARRegion> ptr(xp);
    ptr->covariance.update_rho(rho);
  } else if(covtype == 2 && lptype == 2){
    XPtr<BitsNNGPRegion> ptr(xp);
    ptr->covariance.update_rho(rho);
  }
  
}

// [[Rcpp::export]]
SEXP rtsModel__new_from_bits(SEXP bptr_, SEXP covtype_){
  int covtype = as<int>(covtype_);
  if(covtype == 1){
    XPtr<BitsAR> bptr(bptr_);
    XPtr<ModelAR> ptr(new ModelAR(*bptr),true);
    return ptr;
  } else {
    XPtr<BitsNNGP> bptr(bptr_);
    XPtr<ModelNNGP> ptr(new ModelNNGP(*bptr),true);
    return ptr;
  }
}

// [[Rcpp::export]]
SEXP rtsRegionModel__new_from_bits(SEXP bptr_, SEXP covtype_){
  int covtype = as<int>(covtype_);
  if(covtype == 1){
    XPtr<BitsARRegion> bptr(bptr_);
    XPtr<ModelARRegion> ptr(new ModelARRegion(*bptr),true);
    return ptr;
  } else {
    XPtr<BitsNNGPRegion> bptr(bptr_);
    XPtr<ModelNNGPRegion> ptr(new ModelNNGPRegion(*bptr),true);
    return ptr;
  }
}

// [[Rcpp::export]]
void rtsModel__set_y(SEXP xp, SEXP y_,SEXP covtype_, SEXP lptype_){
  Eigen::VectorXd y = as<Eigen::VectorXd>(y_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->set_y(y);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->set_y(y);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->set_y(y);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->set_y(y);
  }
}

// [[Rcpp::export]]
void rtsModel__set_offset(SEXP xp, SEXP offset_,SEXP covtype_, SEXP lptype_){
  Eigen::VectorXd offset = as<Eigen::VectorXd>(offset_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->set_offset(offset);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->set_offset(offset);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->set_offset(offset);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->set_offset(offset);
  }
}

// [[Rcpp::export]]
void rtsModel__set_weights(SEXP xp, SEXP weights_,SEXP covtype_, SEXP lptype_){
  Eigen::ArrayXd weights = as<Eigen::ArrayXd>(weights_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->set_weights(weights);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->set_weights(weights);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->set_weights(weights);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->set_weights(weights);
  }
}

// [[Rcpp::export]]
void rtsModel__update_beta(SEXP xp, SEXP beta_,SEXP covtype_, SEXP lptype_){
  std::vector<double> beta = as<std::vector<double> >(beta_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->update_beta(beta);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->update_beta(beta);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->update_beta(beta);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->update_beta(beta);
  }
}

// [[Rcpp::export]]
void rtsModel__update_theta(SEXP xp, SEXP theta_,SEXP covtype_, SEXP lptype_){
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->update_theta(theta);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->update_theta(theta);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->update_theta(theta);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->update_theta(theta);
  }
}

// [[Rcpp::export]]
void rtsModel__update_u(SEXP xp, SEXP u_,SEXP covtype_, SEXP lptype_){
  Eigen::MatrixXd u = as<Eigen::MatrixXd>(u_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->update_u(u);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->update_u(u);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->update_u(u);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->update_u(u);
  }
}

// [[Rcpp::export]]
void rtsModel__use_attenuation(SEXP xp, SEXP use_,SEXP covtype_, SEXP lptype_){
  bool use = as<bool>(use_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->matrix.W.attenuated = use;
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->matrix.W.attenuated = use;
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->matrix.W.attenuated = use;
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->matrix.W.attenuated = use;
  }
}

// [[Rcpp::export]]
SEXP rtsModel__get_W(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    VectorXd W = ptr->matrix.W.W();
    return wrap(W);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    VectorXd W = ptr->matrix.W.W();
    return wrap(W);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    VectorXd W = ptr->matrix.W.W();
    return wrap(W);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    VectorXd W = ptr->matrix.W.W();
    return wrap(W);
  }
}

// [[Rcpp::export]]
void rtsModel__ml_theta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.ml_theta();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.ml_theta();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.ml_theta();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.ml_theta();
  }
}

// [[Rcpp::export]]
void rtsModel__ml_beta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.ml_beta();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.ml_beta();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.ml_beta();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.ml_beta();
  }
}

// [[Rcpp::export]]
void rtsModel__ml_all(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.ml_all();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.ml_all();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.ml_all();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.ml_all();
  }
}

// [[Rcpp::export]]
void rtsModel__laplace_ml_beta_u(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.laplace_ml_beta_u();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.laplace_ml_beta_u();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.laplace_ml_beta_u();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.laplace_ml_beta_u();
  }
}

// [[Rcpp::export]]
void rtsModel__laplace_ml_theta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.laplace_ml_theta();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.laplace_ml_theta();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.laplace_ml_theta();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.laplace_ml_theta();
  }
}

// [[Rcpp::export]]
void rtsModel__laplace_ml_beta_theta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.laplace_ml_beta_theta();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.laplace_ml_beta_theta();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.laplace_ml_beta_theta();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.laplace_ml_beta_theta();
  }
}

// [[Rcpp::export]]
void rtsModel__nr_beta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.nr_beta();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.nr_beta();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.nr_beta();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.nr_beta();
  }
}

// [[Rcpp::export]]
void rtsModel__laplace_nr_beta_u(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->optim.laplace_nr_beta_u();
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->optim.laplace_nr_beta_u();
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->optim.laplace_nr_beta_u();
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->optim.laplace_nr_beta_u();
  }
}

// [[Rcpp::export]]
SEXP rtsModel__Sigma(SEXP xp, bool inverse, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd S = ptr->matrix.Sigma(inverse);
    return wrap(S);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd S = ptr->matrix.Sigma(inverse);
    return wrap(S);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd S = ptr->matrix.Sigma(inverse);
    return wrap(S);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd S = ptr->matrix.Sigma(inverse);
    return wrap(S);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__information_matrix(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix();
    return wrap(M);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix();
    return wrap(M);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix();
    return wrap(M);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix();
    return wrap(M);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__u(SEXP xp, bool scaled_, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd u = ptr->re.u(scaled_);
    return wrap(u);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd u = ptr->re.u(scaled_);
    return wrap(u);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd u = ptr->re.u(scaled_);
    return wrap(u);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd u = ptr->re.u(scaled_);
    return wrap(u);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__X(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd X = ptr->model.linear_predictor.X();
    return wrap(X);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd X = ptr->model.linear_predictor.X();
    return wrap(X);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd X = ptr->model.linear_predictor.X();
    return wrap(X);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd X = ptr->model.linear_predictor.X();
    return wrap(X);
  }
}

// [[Rcpp::export]]
void rtsModel__set_trace(SEXP xp, SEXP trace_, SEXP covtype_, SEXP lptype_){
  int trace = as<int>(trace_);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    ptr->set_trace(trace);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    ptr->set_trace(trace);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    ptr->set_trace(trace);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    ptr->set_trace(trace);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__get_beta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::VectorXd beta = ptr->model.linear_predictor.parameter_vector();
    return wrap(beta);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::VectorXd beta = ptr->model.linear_predictor.parameter_vector();
    return wrap(beta);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::VectorXd beta = ptr->model.linear_predictor.parameter_vector();
    return wrap(beta);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::VectorXd beta = ptr->model.linear_predictor.parameter_vector();
    return wrap(beta);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__get_theta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    std::vector<double> theta = ptr->model.covariance.parameters_;
    return wrap(theta);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    std::vector<double> theta = ptr->model.covariance.parameters_;
    return wrap(theta);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    std::vector<double> theta = ptr->model.covariance.parameters_;
    return wrap(theta);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    std::vector<double> theta = ptr->model.covariance.parameters_;
    return wrap(theta);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__L(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd L = ptr->model.covariance.D(true,false);
    return wrap(L);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd L = ptr->model.covariance.D(true,false);
    return wrap(L);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd L = ptr->model.covariance.D(true,false);
    return wrap(L);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd L = ptr->model.covariance.D(true,false);
    return wrap(L);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__ZL(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd ZL = ptr->model.covariance.ZL();
    return wrap(ZL);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd ZL = ptr->model.covariance.ZL();
    return wrap(ZL);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd ZL = ptr->model.covariance.ZL();
    return wrap(ZL);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd ZL = ptr->model.covariance.ZL();
    return wrap(ZL);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__xb(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_); 
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::VectorXd xb = ptr->model.xb();
    return wrap(xb);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::VectorXd xb = ptr->model.xb();
    return wrap(xb);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::VectorXd xb = ptr->model.xb();
    return wrap(xb);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::VectorXd xb = ptr->model.xb();
    return wrap(xb);
  }
}
