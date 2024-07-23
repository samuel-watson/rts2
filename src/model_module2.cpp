#include "rtsheader.h"

using namespace Rcpp;


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
