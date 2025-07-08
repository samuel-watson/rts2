#include "rtsheader.h"
#include <cmath>

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
void rtsModel__update_u(SEXP xp, SEXP u_,bool append,int covtype_, int lptype_){
  Eigen::MatrixXd u = as<Eigen::MatrixXd>(u_);
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&](auto mptr){mptr->update_u(u,append);}
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
SEXP rtsModel__u(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {  return returns(0);}, 
    [](auto mptr){return returns(mptr->re.zu_);}
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
SEXP rtsModel__ar_chol(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.ar_matrix(true));}
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
SEXP rtsModel__L(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->model.covariance.D(true,false));}
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

// [[Rcpp::export]]
SEXP rtsModel__P(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {return returns(0);},
    [](auto mptr){return returns(mptr->model.linear_predictor.P());}
  };
  auto p = std::visit(functor,model.ptr);
  return wrap(std::get<int>(p));
}

// [[Rcpp::export]]
void rtsModel__ml_theta(SEXP xp, int algo, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&algo](auto ptr){
      switch(algo){
      case 1:
        ptr->optim.template ml_theta<NEWUOA>();
        break;
      case 2:
        ptr->optim.template ml_theta<LBFGS>();
        break;
      case 3:
        ptr->optim.template ml_theta<DIRECT>();
        break;
      default:
        ptr->optim.template ml_theta<BOBYQA>();
      break;
      }
    }
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__ml_beta(SEXP xp, int algo, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&algo](auto ptr){
      switch(algo){
      case 1:
        ptr->optim.template ml_beta<NEWUOA>();
        break;
      case 2:
        ptr->optim.template ml_beta<LBFGS>();
        break;
      case 3:
        ptr->optim.template ml_beta<DIRECT>();
        break;
      default:
        ptr->optim.template ml_beta<BOBYQA>();
      break;
      }
    }
  };
  std::visit(functor,model.ptr);
}

// [[Rcpp::export]]
void rtsModel__ml_rho(SEXP xp, int algo, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) {}, 
    [&algo](auto ptr){
      switch(algo){
      case 1:
        ptr->optim.template ml_rho<NEWUOA>();
        break;
      case 2:
        ptr->optim.template ml_rho<LBFGS>();
        break;
      case 3:
        ptr->optim.template ml_rho<DIRECT>();
        break;
      default:
        ptr->optim.template ml_rho<BOBYQA>();
      break;
      }
    }
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

// // [[Rcpp::export]]
// void rtsModel__laplace_nr_beta_u(SEXP xp, int covtype_, int lptype_){
//   TypeSelector model(xp,covtype_,lptype_);
//   auto functor = overloaded {
//     [](int) {}, 
//     [](auto mptr){mptr->optim.laplace_nr_beta_u();}
//   };
//   std::visit(functor,model.ptr);
// }

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
SEXP rtsModel__log_likelihood(SEXP xp, int covtype_, int lptype_){
  TypeSelector model(xp,covtype_,lptype_);
  auto functor = overloaded {
    [](int) { return returns(0);}, 
    [](auto mptr){return returns(mptr->optim.log_likelihood(true));}
  };
  auto ll = std::visit(functor,model.ptr);
  return wrap(std::get<double>(ll));
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
  } else if(covtype == 3 && lptype == 2){
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
  } else if(covtype == 3 && lptype == 3){
    XPtr<ModelHSGPRegionG> ptr(xp);
    Eigen::ArrayXXd M = ptr->intersection_infomat();
    return wrap(M);
  } else {
    Rcpp::stop("Invalid type");
  }
}

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
  } else if(covtype == 3){
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
                                 SEXP theta_, 
                                 SEXP rptr_,
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
                                 SEXP theta_, 
                                 SEXP rptr_,
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
SEXP Model_region_lp__new(SEXP formula_region_, 
                                 SEXP formula_grid_,
                                 SEXP data_region_, 
                                 SEXP data_grid_, 
                                 SEXP colnames_region_, 
                                 SEXP colnames_grid_,
                                 SEXP rptr_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  XPtr<rts::RegionData> rptr(rptr_);
  glmmr::Formula form_r(formula_region);
  glmmr::Formula form_g(formula_grid);
  
  XPtr<rts::regionLinearPredictor> ptr(new rts::regionLinearPredictor(form_r,form_g,data_region,data_grid,colnames_region,colnames_grid,*rptr),true);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_region_lp__P(SEXP ptr_){
  XPtr<rts::regionLinearPredictor> ptr(ptr_);
  Rcpp::Rcout << "\nR P: " << ptr->region_predictor.P() << " G P: " << ptr->grid_predictor.P();
  return wrap(ptr->P());
}

// [[Rcpp::export]]
SEXP Model_hsgp_region_grid_bits__new(SEXP formula_region_, 
                                 SEXP formula_grid_,
                                 SEXP data_region_, 
                                 SEXP data_grid_, 
                                 SEXP colnames_region_, 
                                 SEXP colnames_grid_,
                                 SEXP rptr_,
                                 int T,
                                 int m,
                                 SEXP L_){
  std::string formula_region = as<std::string>(formula_region_);
  std::string formula_grid = as<std::string>(formula_grid_);
  Eigen::ArrayXXd data_region = as<Eigen::ArrayXXd>(data_region_);
  Eigen::ArrayXXd data_grid = as<Eigen::ArrayXXd>(data_grid_);
  std::vector<std::string> colnames_region = as<std::vector<std::string> >(colnames_region_);
  std::vector<std::string> colnames_grid = as<std::vector<std::string> >(colnames_grid_);
  Eigen::ArrayXd L = as<Eigen::ArrayXd>(L_);
  XPtr<rts::RegionData> rptr(rptr_);
  XPtr<BitsHSGPRegion> ptr(new BitsHSGPRegion(formula_region,formula_grid,data_region,data_grid,colnames_region,colnames_grid,T,m,L,*rptr),true);
  glmmr::RandomEffects<BitsHSGPRegion> re(*ptr,data_grid.rows()*T,ptr->covariance.Q());
  ArrayXd xb = ptr->linear_predictor.xb();//ptr->xb();
  Rcpp::Rcout << "\nXb: " << xb.head(10).transpose() << "\nsize: " << xb.size();
  Rcpp::Rcout << "\nOffset size: " << ptr->data.offset.size();
  //glmmr::ModelMatrix<BitsHSGPRegion> matrix(*ptr,re);
  return ptr;
}

// [[Rcpp::export]]
SEXP Model_region_lp__xb(SEXP ptr_, bool grid){
  XPtr<rts::regionLinearPredictor> ptr(ptr_);
  if(grid){
    Eigen::ArrayXd xbg = ptr->grid_predictor.xb();
    return wrap(xbg);
  } else {
    Eigen::ArrayXd xbg = ptr->region_predictor.xb();
    return wrap(xbg);
  }
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

// [[Rcpp::export]]
double max_dist(const Eigen::ArrayXXd &x){
  // this is a brute force algorithm for max distance
  // it can be improved by finding convex hull and then using rotating calipers method
  // but I haven't had the time to implement that!
  int n = x.rows();
  double maxdist = 0;
  double dist = 0;
  for(int i = 1; i < n; i++){
    for(int j = 0; j<(i-1); j++){
      dist = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
      if(dist > maxdist) maxdist = dist;
    }
  }
  return maxdist;
}

// [[Rcpp::export]]
Eigen::ArrayXXd semivariogram(const Eigen::ArrayXXd &x,
                              const Eigen::ArrayXd &offs,
                              const Eigen::ArrayXd &y,
                              int nbins,
                              int nT){
  double maxd = max_dist(x);
  int n = y.size()/nT;
  Eigen::ArrayXd denoms = Eigen::ArrayXd::Zero(nbins);
  Eigen::ArrayXd sums = Eigen::ArrayXd::Zero(nbins);
  double binw = maxd/nbins;
  Eigen::ArrayXd z(y);
  for(int t = 0; t< nT; t++){
    z.segment(t*n,n) *= offs.inverse();
  }
  // Eigen::ArrayXd z = offs.inverse();
  // z *= y;
  double dist;
  // int n = x.rows();
  int binc;
  for(int t = 0; t < nT; t++){
    for(int i = 1; i < n; i++){
      for(int j = 0; j<(i-1); j++){
        dist = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
        binc = static_cast<int>(std::floor(dist/binw));
        denoms(binc) += offs(i)*offs(j);
        sums(binc) += offs(i)*offs(j)*(z(i+n*t)-z(j+n*t))*(z(i+n*t)-z(j+n*t));
      }
    }
  }
  
  denoms *= 2;
  Eigen::ArrayXXd result(nbins,2);
  for(int i=0; i<nbins; i++)result(i,0) = i*binw + binw/2;
  result.col(1) = denoms.inverse()*sums;
  return result;
}