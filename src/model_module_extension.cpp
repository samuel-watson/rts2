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
SEXP rtsModel_cov__log_likelihood(SEXP xp, SEXP covtype_, SEXP lptype_, SEXP u_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  Eigen::VectorXd u = as<Eigen::VectorXd>(u_);
  if(covtype == 1 && lptype == 1){
    XPtr<ModelAR> ptr(xp);
    double aic = ptr->model.covariance.log_likelihood(u);
    return wrap(aic);
  } else if(covtype == 2 && lptype == 1){
    XPtr<ModelNNGP> ptr(xp);
    double aic = ptr->model.covariance.log_likelihood(u);
    return wrap(aic);
  } else if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    double aic = ptr->model.covariance.log_likelihood(u);
    return wrap(aic);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    double aic = ptr->model.covariance.log_likelihood(u);
    return wrap(aic);
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    double aic =ptr->model.covariance.log_likelihood(u);
    return wrap(aic);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
    double aic = ptr->model.covariance.log_likelihood(u);
    return wrap(aic);
  }
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
SEXP nngp_ldlt(SEXP A_, SEXP D_, SEXP NN_){
  Eigen::MatrixXd A = as<Eigen::MatrixXd>(A_);
  Eigen::VectorXd D = as<Eigen::VectorXd>(D_);
  Eigen::ArrayXXi NN = as<Eigen::ArrayXXi>(NN_);
  Eigen::MatrixXd L = rts::inv_ldlt_AD(A,D,NN);
  return wrap(L);
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
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    double aic = ptr->optim.aic();
    return wrap(aic);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
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
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    strvec parnames = ptr->model.linear_predictor.parameter_names();
    return wrap(parnames);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
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
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    strvec parnames = ptr->model.covariance.parameter_names();
    return wrap(parnames);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
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
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    Eigen::MatrixXd sandwich = ptr->matrix.sandwich_matrix();
    return wrap(sandwich);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
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
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix_theta();
    return wrap(M);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
    Eigen::MatrixXd M = ptr->matrix.information_matrix_theta();
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
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    vector_matrix hess = ptr->matrix.re_score();
    return wrap(hess);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
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
  }else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
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
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
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