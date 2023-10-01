#include "rtsheader.h"

namespace Rcpp {
template<>
SEXP wrap(const vector_matrix& x){
  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("vec") = Rcpp::wrap(x.vec),
      Rcpp::Named("mat") = Rcpp::wrap(x.mat)
  ));
}

template<>
SEXP wrap(const matrix_matrix& x){
  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("mat1") = Rcpp::wrap(x.mat1),
      Rcpp::Named("mat2") = Rcpp::wrap(x.mat2),
      Rcpp::Named("a") = Rcpp::wrap(x.a),
      Rcpp::Named("b") = Rcpp::wrap(x.b)
  ));
}

template<>
SEXP wrap(const kenward_data& x){
  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("vcov_beta") = Rcpp::wrap(x.vcov_beta),
      Rcpp::Named("vcov_theta") = Rcpp::wrap(x.vcov_theta),
      Rcpp::Named("dof") = Rcpp::wrap(x.dof)
  ));
}
}

using namespace Rcpp;

// [[Rcpp::export]]
SEXP rtsModel_nngp__A(SEXP ptr_, SEXP lptype_){
  int lptype = as<int>(lptype_);
  if(lptype == 1){
    XPtr<ModelNNGP> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.A;
    return wrap(A);
  } else if(lptype == 2){
    XPtr<ModelNNGPRegion> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.A;
    return wrap(A);
  } else if(lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.A;
    return wrap(A);
  }
}

// [[Rcpp::export]]
SEXP rtsModel_nngp__D(SEXP ptr_, SEXP lptype_){
  int lptype = as<int>(lptype_);
  if(lptype == 1){
    XPtr<ModelNNGP> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.Dvec;
    return wrap(A);
  } else if(lptype == 2){
    XPtr<ModelNNGPRegion> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.Dvec;
    return wrap(A);
  } else if(lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.Dvec;
    return wrap(A);
  }
}

// [[Rcpp::export]]
SEXP rtsModel_nngp__submatrix(SEXP ptr_, SEXP lptype_, SEXP i_){
  int lptype = as<int>(lptype_);
  int i = as<int>(i_);
  if(lptype == 1){
    XPtr<ModelNNGP> ptr(ptr_);
    vector_matrix A = ptr->model.covariance.submatrix(i);
    return wrap(A);
  } else if(lptype == 2){
    XPtr<ModelNNGPRegion> ptr(ptr_);
    vector_matrix A = ptr->model.covariance.submatrix(i);
    return wrap(A);
  } else if(lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(ptr_);
    vector_matrix A = ptr->model.covariance.submatrix(i);
    return wrap(A);
  }
}

// [[Rcpp::export]]
SEXP rtsModel_cov__log_likelihood(SEXP xp, int covtype_, int lptype_, SEXP u_){
  Eigen::VectorXd u = as<Eigen::VectorXd>(u_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      return 0.0;
    }, 
    [&u](auto mptr){return mptr->model.covariance.log_likelihood(u);}
  };
  double ll = std::visit(functor,model.ptr);
  return wrap(ll);
}

// [[Rcpp::export]]
void rtsModel_cov__set_sparse(SEXP ptr_, SEXP lptype_, SEXP sparse_){
  int lptype = as<int>(lptype_);
  bool sparse = as<bool>(sparse_);
  if(lptype == 1){
    XPtr<ModelAR> ptr(ptr_);
    ptr->model.covariance.set_sparse(sparse);
  } else if(lptype == 2){
    XPtr<ModelARRegion> ptr(ptr_);
    ptr->model.covariance.set_sparse(sparse);
  } else if(lptype == 3){
    XPtr<ModelARRegionG> ptr(ptr_);
    ptr->model.covariance.set_sparse(sparse);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__aic(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      return 0.0;
    }, 
    [](auto mptr){return mptr->optim.aic();}
  };
  double aic = std::visit(functor,model.ptr);
  return wrap(aic);
}

// [[Rcpp::export]]
SEXP rtsModel__beta_parameter_names(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      strvec val = {"a"};
      return val;
    }, 
    [](auto mptr){return mptr->model.linear_predictor.parameter_names();}
  };
  strvec parnames = std::visit(functor,model.ptr);
  return wrap(parnames);
}

// [[Rcpp::export]]
SEXP rtsModel__theta_parameter_names(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      strvec val = {"a"};
      return val;
    }, 
    [](auto mptr){return mptr->model.covariance.parameter_names();}
  };
  strvec parnames = std::visit(functor,model.ptr);
  return wrap(parnames);
}

// [[Rcpp::export]]
SEXP rtsModel__infomat_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      Eigen::MatrixXd S = Eigen::MatrixXd::Zero(1,1);
      return S;
    }, 
    [](auto mptr){return mptr->matrix.information_matrix_theta();}
  };
  Eigen::MatrixXd S = std::visit(functor,model.ptr);
  return wrap(S);
}

// [[Rcpp::export]]
SEXP rtsModel__hessian(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      vector_matrix hess(1);
      return hess;
    }, 
    [](auto mptr){return mptr->matrix.re_score();}
  };
  vector_matrix hess = std::visit(functor,model.ptr);
  return wrap(hess);
}

// [[Rcpp::export]]
SEXP rtsModel__region_intensity(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelHSGPRegion> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__region_grid_xb(SEXP xp, SEXP covtype_){
  int covtype = as<int>(covtype_);
  if(covtype == 1){
    XPtr<ModelARRegionG> ptr(xp);
    Eigen::VectorXd intens = ptr->model.linear_predictor.grid_predictor.xb();
    return wrap(intens);
  } else if(covtype == 2){
    XPtr<ModelNNGPRegionG> ptr(xp);
    Eigen::VectorXd intens = ptr->model.linear_predictor.grid_predictor.xb();
    return wrap(intens);
  } else if(covtype == 2){
    XPtr<ModelHSGPRegionG> ptr(xp);
    Eigen::VectorXd intens = ptr->model.linear_predictor.grid_predictor.xb();
    return wrap(intens);
  } 
}

// [[Rcpp::export]]
SEXP rtsModel__grid_to_region(SEXP xp, SEXP u_){
  Eigen::MatrixXd u = as<Eigen::MatrixXd>(u_);
  XPtr<rts::RegionData> rptr(xp);
  Eigen::MatrixXd out = rptr->grid_to_region(u);
  return wrap(out);
}

// [[Rcpp::export]]
SEXP rtsModel__predict(SEXP xp, SEXP newdata_,
                    SEXP newoffset_,
                    int m, int covtype_, int lptype_){
  Eigen::ArrayXXd newdata = Rcpp::as<Eigen::ArrayXXd>(newdata_);
  Eigen::ArrayXd newoffset = Rcpp::as<Eigen::ArrayXd>(newoffset_);
  Eigen::MatrixXd samps(newdata.rows(),m>0 ? m : 1);
  
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {
      vector_matrix res(1);
      return res;
    }, 
    [&](auto mptr){return mptr->re.predict_re(newdata,newoffset);}
  };
  auto functor2 = overloaded {
    [](int) {
      Eigen::VectorXd S = Eigen::VectorXd::Zero(1);
      return res;
    }, 
    [&](auto mptr){return mptr->model.linear_predictor.predict_xb(newdata,newoffset);}
  };
  
  vector_matrix res = std::visit(functor,model.ptr);
  Eigen::VectorXd xb = std::visit(functor2,model.ptr);
  if(m>0){
    samps = glmmr::maths::sample_MVN(res,m);
  } else {
    samps.setZero();
  }
  return Rcpp::List::create(
    Rcpp::Named("linear_predictor") = wrap(xb),
    Rcpp::Named("re_parameters") = wrap(res),
    Rcpp::Named("samples") = wrap(samps)
  );
  
}