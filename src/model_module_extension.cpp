#include "rtsheader.h"

namespace Rcpp {
template<>
SEXP wrap(const VectorMatrix& x){
  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("vec") = Rcpp::wrap(x.vec),
      Rcpp::Named("mat") = Rcpp::wrap(x.mat)
  ));
}

template<>
SEXP wrap(const MatrixMatrix& x){
  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("mat1") = Rcpp::wrap(x.mat1),
      Rcpp::Named("mat2") = Rcpp::wrap(x.mat2),
      Rcpp::Named("a") = Rcpp::wrap(x.a),
      Rcpp::Named("b") = Rcpp::wrap(x.b)
  ));
}

template<>
SEXP wrap(const sparse& x){
  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("Ap") = Rcpp::wrap(x.Ap),
      Rcpp::Named("Ai") = Rcpp::wrap(x.Ai),
      Rcpp::Named("Ax") = Rcpp::wrap(x.Ax)
  ));
}

template<typename T1, typename T2> SEXP wrap( const std::pair<T1,T2>& _v ) {
  return Rcpp::List::create(
    Rcpp::Named("first")  = Rcpp::wrap<T1>( _v.first ),
    Rcpp::Named("second") = Rcpp::wrap<T2>( _v.second )
  );
};

}

using namespace Rcpp;

// [[Rcpp::export]]
SEXP rtsModel__get_log_likelihood_values(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {  return returnType(0);}, 
    [](auto ptr){return returnType(ptr->optim.current_likelihood_values());}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<std::pair<double,double> >(S));
}

// [[Rcpp::export]]
SEXP rtsModel__u_diagnostic(SEXP xp, int covtype_, int lptype_){
 TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {  return returnType(0);}, 
    [](auto ptr){return returnType(ptr->optim.u_diagnostic());}
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<std::pair<double,double> >(S));
}

// [[Rcpp::export]]
void rtsModel__saem(SEXP xp, bool saem_, int block_size, double alpha, bool pr_average, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&](auto ptr){
      ptr->optim.control.saem = saem_;
      ptr->optim.control.alpha = alpha;
      ptr->re.mcmc_block_size = block_size;
      ptr->optim.control.pr_average = pr_average;
    }
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
SEXP rtsModel__ll_diff_variance(SEXP xp, bool beta, bool theta, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returnType(0);}, 
    [&](auto ptr){
      return returnType(ptr->optim.ll_diff_variance(beta,theta));
    }
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<double>(S));
}

// // [[Rcpp::export]]
// SEXP rtsModel__hess_and_grad(SEXP xp, int covtype_, int lptype_){
//   TypeSelector model(xp,covtype_,lptype_);
//   auto functor = overloaded {
//     [](int) {return returns(0);}, 
//     [](auto mptr){return returns(mptr->matrix.hess_and_grad());}
//   };
//   auto S = std::visit(functor,model.ptr);
//   return wrap(std::get<MatrixMatrix>(S));
// }

// [[Rcpp::export]]
void rtsModel__set_bobyqa_control(SEXP xp, int covtype_, int lptype_,
                                  SEXP npt_, SEXP rhobeg_, SEXP rhoend_){
  TypeSelector model(xp,covtype_,lptype_);
  int npt = as<int>(npt_);
  double rhobeg = as<double>(rhobeg_);
  double rhoend = as<double>(rhoend_);
  auto functor = overloaded {
    [](int) {}, 
    [&](auto mptr){
      mptr->optim.set_bobyqa_control(npt,rhobeg,rhoend);
    }
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__set_bound(SEXP xp, int covtype_, int lptype_, SEXP bound_, bool lower = true){
  TypeSelector model(xp,covtype_,lptype_);
  std::vector<double> bound = as<std::vector<double> >(bound_);
  auto functor = overloaded {
      [](int) {}, 
      [&](auto ptr){ptr->optim.set_bound(bound,lower);}
  };
  std::visit(functor,model.ptr);
}

// // [[Rcpp::export]]
// void rtsModel__set_cov_bobyqa_control(SEXP xp, int covtype_, int lptype_,
//                                       SEXP rhobeg_, SEXP rhoend_){
//   TypeSelector model(xp,covtype_,lptype_);
//   double rhobeg = as<double>(rhobeg_);
//   double rhoend = as<double>(rhoend_);
//   auto functor = overloaded {
//     [](int) {}, 
//     [&](auto mptr){
//       mptr->optim.set_cov_bobyqa_control(rhobeg,rhoend);
//     }
//   };
//   std::visit(functor,model.ptr);
// }

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
  } else {
    Rcpp::stop("Invalid lp type.");
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
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}

// [[Rcpp::export]]
SEXP rtsModel_nngp__submatrix(SEXP ptr_, SEXP lptype_, SEXP i_){
  int lptype = as<int>(lptype_);
  int i = as<int>(i_);
  if(lptype == 1){
    XPtr<ModelNNGP> ptr(ptr_);
    VectorMatrix A = ptr->model.covariance.submatrix(i);
    return wrap(A);
  } else if(lptype == 2){
    XPtr<ModelNNGPRegion> ptr(ptr_);
    VectorMatrix A = ptr->model.covariance.submatrix(i);
    return wrap(A);
  } else if(lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(ptr_);
    VectorMatrix A = ptr->model.covariance.submatrix(i);
    return wrap(A);
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}

// [[Rcpp::export]]
SEXP rtsModel_hsgp__Phi(SEXP ptr_, SEXP lptype_, bool lambda, bool inverse){
  int lptype = as<int>(lptype_);
  if(lptype == 1){
    XPtr<ModelHSGP> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.PhiSPD(lambda,inverse);
    return wrap(A);
  } else if(lptype == 2){
    XPtr<ModelHSGPRegion> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.PhiSPD(lambda,inverse);
    return wrap(A);
  } else if(lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(ptr_);
    Eigen::MatrixXd A = ptr->model.covariance.PhiSPD(lambda,inverse);
    return wrap(A);
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}

// [[Rcpp::export]]
SEXP rtsModel_hsgp__Lambda(SEXP ptr_, SEXP lptype_){
  int lptype = as<int>(lptype_);
  if(lptype == 1){
    XPtr<ModelHSGP> ptr(ptr_);
    Eigen::ArrayXd A = ptr->model.covariance.LambdaSPD();
    return wrap(A);
  } else if(lptype == 2){
    XPtr<ModelHSGPRegion> ptr(ptr_);
    Eigen::ArrayXd A = ptr->model.covariance.LambdaSPD();
    return wrap(A);
  } else if(lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(ptr_);
    Eigen::ArrayXd A = ptr->model.covariance.LambdaSPD();
    return wrap(A);
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}

// [[Rcpp::export]]
void rtsModel_hsgp__set_function(SEXP ptr_, SEXP lptype_, bool sqexp){
  int lptype = as<int>(lptype_);
  if(lptype == 1){
    XPtr<ModelHSGP> ptr(ptr_);
    ptr->model.covariance.set_function(sqexp);
  } else if(lptype == 2){
    XPtr<ModelHSGPRegion> ptr(ptr_);
    ptr->model.covariance.set_function(sqexp);
  } else if(lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(ptr_);
    ptr->model.covariance.set_function(sqexp);
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}

// [[Rcpp::export]]
SEXP rtsModel_cov__log_likelihood(SEXP xp, int covtype_, int lptype_, SEXP u_){
  Eigen::VectorXd u = as<Eigen::VectorXd>(u_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0); }, 
    [&u](auto mptr){return returns(mptr->model.covariance.log_likelihood(u));}
  };
  auto ll = std::visit(functor,model.ptr);
  return wrap(std::get<double>(ll));
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
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}

// [[Rcpp::export]]
SEXP rtsModel__aic(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);}, 
    [](auto mptr){return returns(mptr->optim.aic());}
  };
  auto aic = std::visit(functor,model.ptr);
  return wrap(std::get<double>(aic));
}

// [[Rcpp::export]]
SEXP rtsModel__beta_parameter_names(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);}, 
    [](auto mptr){return returns(mptr->model.linear_predictor.parameter_names());}
  };
  auto parnames = std::visit(functor,model.ptr);
  return wrap(std::get<std::vector<std::string> >(parnames));
}

// [[Rcpp::export]]
SEXP rtsModel__theta_parameter_names(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.parameter_names());}
  };
  auto parnames = std::visit(functor,model.ptr);
  return wrap(std::get<std::vector<std::string> >(parnames));
}

// [[Rcpp::export]]
SEXP rtsModel__infomat_theta(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){
#ifdef GLMMR11
      return returns(mptr->matrix.template information_matrix_theta<glmmr::IM::EIM>());
#else
      return returns(mptr->matrix.information_matrix_theta());
#endif
    }
  };
  auto S = std::visit(functor,model.ptr);
  return wrap(std::get<Eigen::MatrixXd>(S));
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
  } else if(covtype == 3 && lptype == 2){
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
  } else if(covtype == 3 && lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.region_intensity();
    return wrap(intens);
  } else {
    Rcpp::stop("Invalid lp type.");
  }
}


// [[Rcpp::export]]
SEXP rtsModel__y_pred(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.y_predicted(false);
    return wrap(intens);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.y_predicted(false);
    return wrap(intens);
  } else if(covtype == 3 && lptype == 2){
    XPtr<ModelHSGPRegion> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.y_predicted(false);
    return wrap(intens);
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.y_predicted(false);
    return wrap(intens);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.y_predicted(false);
    return wrap(intens);
  } else if(covtype == 3 && lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(xp);
    Eigen::ArrayXXd intens = ptr->optim.y_predicted(false);
    return wrap(intens);
  } else {
    Rcpp::stop("Invalid type.");
  }
}

// [[Rcpp::export]]
SEXP rtsModel__grid_to_region_multiplier_matrix(SEXP xp, SEXP covtype_, SEXP lptype_){
  int covtype = as<int>(covtype_);
  int lptype = as<int>(lptype_);
  if(covtype == 1 && lptype == 2){
    XPtr<ModelARRegion> ptr(xp);
    sparse P = ptr->grid_to_region_multiplier_matrix();
    return wrap(P);
  } else if(covtype == 2 && lptype == 2){
    XPtr<ModelNNGPRegion> ptr(xp);
    sparse P = ptr->grid_to_region_multiplier_matrix();
    return wrap(P);
  } else if(covtype == 3 && lptype == 2){
    XPtr<ModelHSGPRegion> ptr(xp);
    sparse P = ptr->grid_to_region_multiplier_matrix();
    return wrap(P);
  } else if(covtype == 1 && lptype == 3){
    XPtr<ModelARRegionG> ptr(xp);
    sparse P = ptr->grid_to_region_multiplier_matrix();
    return wrap(P);
  } else if(covtype == 2 && lptype == 3){
    XPtr<ModelNNGPRegionG> ptr(xp);
    sparse P = ptr->grid_to_region_multiplier_matrix();
    return wrap(P);
  } else if(covtype == 3 && lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(xp);
    sparse P = ptr->grid_to_region_multiplier_matrix();
    return wrap(P);
  } else {
    Rcpp::stop("Invalid type.");
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
  }  else {
    Rcpp::stop("Invalid cov type.");
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
  Rcpp::stop("Predicting at new locations is currently disabled for LGCPs in this package. It will \
               be updated for the new covariance functions in a future update.");
  // 
  // Eigen::ArrayXXd newdata = Rcpp::as<Eigen::ArrayXXd>(newdata_);
  // Eigen::ArrayXd newoffset = Rcpp::as<Eigen::ArrayXd>(newoffset_);
  // Eigen::MatrixXd samps(newdata.rows(),m>0 ? m : 1);
  // 
  // TypeSelector model(xp,covtype_,lptype_);
  // auto functor = overloaded {
  //   [](int) {return returns(0);}, 
  //   [&](auto mptr){return returns(mptr->re.predict_re(newdata,newoffset));}
  // };
  // auto functor2 = overloaded {
  //   [](int) {return returns(0); }, 
  //   [&](auto mptr){return returns(mptr->model.linear_predictor.predict_xb(newdata,newoffset));}
  // };
  // 
  // auto res = std::visit(functor,model.ptr);
  // auto xb = std::visit(functor2,model.ptr);
  // if(m>0){
  //   samps = glmmr::maths::sample_MVN(std::get<VectorMatrix>(res),m);
  // } else {
  //   samps.setZero();
  // }
  // return Rcpp::List::create(
  //   Rcpp::Named("linear_predictor") = wrap(std::get<Eigen::VectorXd>(xb)),
  //   Rcpp::Named("re_parameters") = wrap(std::get<VectorMatrix>(res)),
  //   Rcpp::Named("samples") = wrap(samps)
  // );
  
}