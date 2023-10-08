#include "rtsheader.h"

using namespace Rcpp;

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
    [](int) {return returns(0);  }, 
    [](auto mptr){return returns(mptr->matrix.W.W());}
  };
  auto W = std::visit(functor,model.ptr);
  return wrap(std::get<VectorXd>(W));
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
    [](int) { return returns(0);}, 
    [inverse](auto mptr){return returns(mptr->matrix.Sigma(inverse));}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::MatrixXd>(S));
}

// [[Rcpp::export]]
SEXP rtsModel__information_matrix(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->matrix.Sigma(true));}
  };
  auto functorx = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.linear_predictor.X());}
  };
  auto S = std::visit(functor,model.ptr);
  auto X = std::visit(functorx,model.ptr);
  Eigen::MatrixXd Sigma = std::get< Eigen::MatrixXd>(S);
  Eigen::MatrixXd Xm = std::get< Eigen::MatrixXd>(X);
  Eigen::MatrixXd M = Xm.transpose()*Sigma*Xm;
  M = M.llt().solve(Eigen::MatrixXd::Identity(M.rows(),M.cols()));
  return wrap(M);
}

// [[Rcpp::export]]
SEXP rtsModel__information_matrix_region(SEXP xp, int covtype, int lptype){
  if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelHSGPRegion> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__u(SEXP xp, bool scaled_, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {  return returns(0);}, 
    [scaled_](auto mptr){return returns(mptr->re.u(scaled_));}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::MatrixXd>(S));
}

// [[Rcpp::export]]
SEXP rtsModel__X(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.linear_predictor.X());}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::MatrixXd>(S));
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
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.linear_predictor.parameter_vector());}
  };
  auto beta = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::VectorXd>(beta));
}

// [[Rcpp::export]]
SEXP rtsModel__get_rho(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.rho);}
  };
  auto beta = std::visit(functor,model.ptr);
  return wrap(std::get<double>(beta));
}

// [[Rcpp::export]]
SEXP rtsModel__log_likelihood(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->optim.log_likelihood());}
  };
  auto ll = std::visit(functor,model.ptr);
  return wrap(std::get<double>(ll));
}

// [[Rcpp::export]]
SEXP rtsModel__get_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.parameters_);}
  };
  auto theta = std::visit(functor,model.ptr);
  return wrap(std::get<std::vector<double> >(theta));
}

// [[Rcpp::export]]
SEXP rtsModel__ZL(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.ZL());}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::MatrixXd>(S));
}

// [[Rcpp::export]]
SEXP rtsModel__D(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.D(false,false));}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::MatrixXd>(S));
}

// [[Rcpp::export]]
SEXP rtsModel__xb(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);}, 
    [](auto mptr){return returns(Eigen::VectorXd((mptr->model.xb()).matrix()));}
  };
  auto xb = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::VectorXd>(xb));
}
