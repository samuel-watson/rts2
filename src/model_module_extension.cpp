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
  }
}

// [[Rcpp::export]]
SEXP rtsModel__aic(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    double aic = ptr->optim.aic();
    return wrap(aic);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    double aic = ptr->optim.aic();
    return wrap(aic);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    double aic = ptr->optim.aic();
    return wrap(aic);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    double aic = ptr->optim.aic();
    return wrap(aic);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__beta_parameter_names(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    strvec parnames = ptr->model.linear_predictor.parameter_names();
    return wrap(parnames);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    strvec parnames = ptr->model.linear_predictor.parameter_names();
    return wrap(parnames);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    strvec parnames = ptr->model.linear_predictor.parameter_names();
    return wrap(parnames);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    strvec parnames = ptr->model.linear_predictor.parameter_names();
    return wrap(parnames);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__theta_parameter_names(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    strvec parnames = ptr->model.covariance.parameter_names();
    return wrap(parnames);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    strvec parnames = ptr->model.covariance.parameter_names();
    return wrap(parnames);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    strvec parnames = ptr->model.covariance.parameter_names();
    return wrap(parnames);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    strvec parnames = ptr->model.covariance.parameter_names();
    return wrap(parnames);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__sandwich(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd sandwich = ptr->matrix.sandwich_matrix();
    return wrap(sandwich);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd sandwich = ptr->matrix.sandwich_matrix();
    return wrap(sandwich);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd sandwich = ptr->matrix.sandwich_matrix();
    return wrap(sandwich);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd sandwich = ptr->matrix.sandwich_matrix();
    return wrap(sandwich);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__infomat_theta(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix_theta();
    return wrap(M);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix_theta();
    return wrap(M);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix_theta();
    return wrap(M);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix_theta();
    return wrap(M);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__kenward_roger(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    kenward_data M = ptr->matrix.kenward_roger();
    return wrap(M);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    kenward_data M = ptr->matrix.kenward_roger();
    return wrap(M);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    kenward_data M = ptr->matrix.kenward_roger();
    return wrap(M);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    kenward_data M = ptr->matrix.kenward_roger();
    return wrap(M);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__hessian(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    vector_matrix hess = ptr->matrix.re_score();
    return wrap(hess);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    vector_matrix hess = ptr->matrix.re_score();
    return wrap(hess);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    vector_matrix hess = ptr->matrix.re_score();
    return wrap(hess);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    vector_matrix hess = ptr->matrix.re_score();
    return wrap(hess);
  }
}

// [[Rcpp::export]]
SEXP rtsModel__predict(SEXP xp, SEXP newdata_,
                    SEXP newoffset_,
                    int m, SEXP covtype_, SEXP lptype_){
  Eigen::ArrayXXd newdata = Rcpp::as<Eigen::ArrayXXd>(newdata_);
  Eigen::ArrayXd newoffset = Rcpp::as<Eigen::ArrayXd>(newoffset_);
  Eigen::MatrixXd samps(newdata.rows(),m>0 ? m : 1);
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    vector_matrix res = ptr->re.predict_re(newdata,newoffset);
    if(m>0){
      samps = glmmr::maths::sample_MVN(res,m);
    } else {
      samps.setZero();
    }
    Eigen::VectorXd xb = ptr->model.linear_predictor.predict_xb(newdata,newoffset);
    return Rcpp::List::create(
      Rcpp::Named("linear_predictor") = wrap(xb),
      Rcpp::Named("re_parameters") = wrap(res),
      Rcpp::Named("samples") = wrap(samps)
    );
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    vector_matrix res = ptr->re.predict_re(newdata,newoffset);
    if(m>0){
      samps = glmmr::maths::sample_MVN(res,m);
    } else {
      samps.setZero();
    }
    Eigen::VectorXd xb = ptr->model.linear_predictor.predict_xb(newdata,newoffset);
    return Rcpp::List::create(
      Rcpp::Named("linear_predictor") = wrap(xb),
      Rcpp::Named("re_parameters") = wrap(res),
      Rcpp::Named("samples") = wrap(samps)
    );
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    vector_matrix res = ptr->re.predict_re(newdata,newoffset);
    if(m>0){
      samps = glmmr::maths::sample_MVN(res,m);
    } else {
      samps.setZero();
    }
    Eigen::VectorXd xb = ptr->model.linear_predictor.predict_xb(newdata,newoffset);
    return Rcpp::List::create(
      Rcpp::Named("linear_predictor") = wrap(xb),
      Rcpp::Named("re_parameters") = wrap(res),
      Rcpp::Named("samples") = wrap(samps)
    );
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    vector_matrix res = ptr->re.predict_re(newdata,newoffset);
    if(m>0){
      samps = glmmr::maths::sample_MVN(res,m);
    } else {
      samps.setZero();
    }
    Eigen::VectorXd xb = ptr->model.linear_predictor.predict_xb(newdata,newoffset);
    return Rcpp::List::create(
      Rcpp::Named("linear_predictor") = wrap(xb),
      Rcpp::Named("re_parameters") = wrap(res),
      Rcpp::Named("samples") = wrap(samps)
    );
  }
}