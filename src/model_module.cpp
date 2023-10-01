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
SEXP Model_ar_lp__new(SEXP formula_, SEXP data_, SEXP grid_data_, SEXP colnames_,
                        SEXP beta_,
                        SEXP theta_, SEXP T_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  XPtr<ModelAR>ptr(new ModelAR(formula,data,grid_data,colnames,T),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_nngp_lp__new(SEXP formula_, SEXP data_,SEXP grid_data_, SEXP colnames_,
                          SEXP beta_,
                          SEXP theta_, SEXP T_, SEXP m_, SEXP gptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  XPtr<rts::griddata> gptr(gptr_);
  XPtr<ModelNNGP> ptr(new ModelNNGP(formula,data,grid_data,colnames,T,m,*gptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_hsgp_lp__new(SEXP formula_, SEXP data_,SEXP grid_data_, SEXP colnames_,
                       SEXP beta_, SEXP theta_, SEXP T_, SEXP m_, SEXP L_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
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
                               SEXP theta_, SEXP rptr_,
                               SEXP T_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
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
                                 SEXP T_, SEXP m_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
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
                                 SEXP T_, SEXP m_,
                                 SEXP L_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  Eigen::ArrayXd L = as<Eigen::ArrayXd>(L_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelHSGPRegionG> ptr(new ModelHSGPRegionG(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,T,m,L,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_ar_region__new(SEXP formula_, SEXP data_, SEXP grid_data_, SEXP colnames_,
                          SEXP beta_,SEXP theta_, SEXP T_, SEXP rptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelARRegion> ptr(new ModelARRegion(formula,data,grid_data,colnames,T,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_nngp_region__new(SEXP formula_, SEXP data_,SEXP grid_data_, SEXP colnames_,
                        SEXP beta_,SEXP theta_, SEXP T_, SEXP m_, SEXP rptr_, SEXP gptr_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<rts::griddata> gptr(gptr_);
  XPtr<ModelNNGPRegion> ptr(new ModelNNGPRegion(formula,data,grid_data,colnames,T,m,*rptr,*gptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_hsgp_region__new(SEXP formula_, SEXP data_,SEXP grid_data_, SEXP colnames_,
                             SEXP beta_,SEXP theta_, SEXP T_, SEXP m_, SEXP rptr_, SEXP L_){
  std::string formula = as<std::string>(formula_);
  Eigen::ArrayXXd data = as<Eigen::ArrayXXd>(data_);
  Eigen::ArrayXXd grid_data = as<Eigen::ArrayXXd>(grid_data_);
  std::vector<std::string> colnames = as<std::vector<std::string> >(colnames_);
  std::vector<double> beta = as<std::vector<double> >(beta_);
  std::vector<double> theta = as<std::vector<double> >(theta_);
  int T = as<int>(T_);
  int m = as<int>(m_);
  Eigen::ArrayXd L = as<Eigen::ArrayXd>(L_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<ModelHSGPRegion> ptr(new ModelHSGPRegion(formula,data,grid_data,colnames,T,m,L,*rptr),true);
  ptr->model.linear_predictor.update_parameters(beta);
  ptr->model.covariance.update_parameters(theta);
  return ptr;
}

// [[Rcpp::export]]
void rtsModel__set_y(SEXP xp, SEXP y_,int covtype_, int lptype_){
  Eigen::VectorXd y = as<Eigen::VectorXd>(y_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&y](auto mptr){mptr->set_y(y);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__set_offset(SEXP xp, SEXP offset_,int covtype_, int lptype_){
  Eigen::VectorXd offset = as<Eigen::VectorXd>(offset_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&offset](auto mptr){mptr->set_offset(offset);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__set_weights(SEXP xp, SEXP weights_,int covtype_, int lptype_){
  Eigen::ArrayXd weights = as<Eigen::ArrayXd>(weights_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&weights](auto mptr){mptr->set_weights(weights);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__update_beta(SEXP xp, SEXP beta_,int covtype_, int lptype_){
  std::vector<double> beta = as<std::vector<double> >(beta_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&beta](auto mptr){mptr->update_beta(beta);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__update_rho(SEXP xp, double rho_,int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&rho_](auto mptr){mptr->update_rho(rho_);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__update_theta(SEXP xp, SEXP theta_,int covtype_, int lptype_){
  std::vector<double> theta = as<std::vector<double> >(theta_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&theta](auto mptr){mptr->update_theta(theta);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__update_u(SEXP xp, SEXP u_,int covtype_, int lptype_){
  Eigen::MatrixXd u = as<Eigen::MatrixXd>(u_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&u](auto mptr){mptr->update_u(u);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__use_attenuation(SEXP xp, SEXP use_,int covtype_, int lptype_){
  bool use = as<bool>(use_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&use](auto mptr){mptr->matrix.W.attenuated = use;}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
SEXP rtsModel__get_W(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      VectorXd tmp = VectorXd::Constant(1);
      return tmp;
    }, 
    [](auto mptr){return mptr->matrix.W.W();}
  };
  VectorXd W = std::visit(functor,model.ptr);
  return wrap(W);
}

// [[Rcpp::export]]
void rtsModel__ml_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.ml_theta();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__ml_beta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.ml_beta();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__ml_rho(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.ml_rho();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__laplace_ml_beta_u(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.laplace_ml_beta_u();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__laplace_ml_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.laplace_ml_theta();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__laplace_ml_beta_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.laplace_ml_beta_theta();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__nr_beta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.nr_beta();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__laplace_nr_beta_u(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [](auto mptr){mptr->optim.laplace_nr_beta_u();}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
SEXP rtsModel__Sigma(SEXP xp, bool inverse, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [inverse](auto mptr){return mptr->matrix.Sigma(inverse);}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
SEXP rtsModel__information_matrix(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [](auto mptr){return mptr->matrix.information_matrix();}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
SEXP rtsModel__u(SEXP xp, bool scaled_, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [scaled_](auto mptr){return mptr->re.u(scaled_);}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
SEXP rtsModel__X(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [](auto mptr){return mptr->model.linear_predictor.X();}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
void rtsModel__set_trace(SEXP xp, SEXP trace_, int covtype_, int lptype_){
  int trace = as<int>(trace_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [trace](auto mptr){mptr->set_trace(trace);}
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
SEXP rtsModel__get_beta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::VectorXd S = Eigen::VectorXd::Zero(1);
      return S;
    }, 
    [](auto mptr){return mptr->model.linear_predictor.parameter_vector();}
  };
  Eigen::VectorXd beta = std::visit(functor,model.ptr);
  return wrap(beta);
}

// [[Rcpp::export]]
SEXP rtsModel__get_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      std::vector<double> theta = {0};
      return theta;
    }, 
    [](auto mptr){return mptr->model.covariance.parameters_;}
  };
  std::vector<double> theta = std::visit(functor,model.ptr);
  return wrap(theta);
}

// [[Rcpp::export]]
SEXP rtsModel__ZL(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [](auto mptr){return mptr->model.covariance.ZL();}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
SEXP rtsModel__D(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [](auto mptr){return mptr->model.covariance.D(false,false);}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
SEXP rtsModel__xb(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::VectorXd S = Eigen::VectorXd::Zero(1);
      return S;
    }, 
    [](auto mptr){return mptr->model.xb();}
  };
  Eigen::VectorXd xb = std::visit(functor,model.ptr);
  return wrap(xb);
}
