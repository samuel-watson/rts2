#pragma once

#include "rts/rtsmaths.h"
#include "rts/rtsmodel.h"
#include "rts/rtsregionmodel.h"
#include <variant>

typedef rts::rtsModel<rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> > ModelAR;
typedef rts::rtsModel<rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> > ModelNNGP;
typedef rts::rtsModel<rts::rtsModelBits<rts::hsgpCovariance, glmmr::LinearPredictor> > ModelHSGP;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> > ModelARRegionG;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> > ModelNNGPRegionG;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::hsgpCovariance, rts::regionLinearPredictor> > ModelHSGPRegionG;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> > ModelARRegion;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> > ModelNNGPRegion;
typedef rts::rtsRegionModel<rts::rtsModelBits<rts::hsgpCovariance, glmmr::LinearPredictor> > ModelHSGPRegion;

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

struct TypeSelector
{
  std::variant<int, Rcpp::XPtr<ModelAR>, Rcpp::XPtr<ModelNNGP>, Rcpp::XPtr<ModelHSGP>, Rcpp::XPtr<ModelARRegion>, Rcpp::XPtr<ModelNNGPRegion>, 
                  Rcpp::XPtr<ModelHSGPRegion>, Rcpp::XPtr<ModelARRegionG>, Rcpp::XPtr<ModelNNGPRegionG>, Rcpp::XPtr<ModelHSGPRegionG> > ptr; 
  TypeSelector(SEXP xp, int covtype, int lptype) : ptr(0) {
    if(covtype == 1 && lptype == 1){
      Rcpp::XPtr<ModelAR> newptr(xp);
      ptr = newptr;
    } else if(covtype == 2 && lptype == 1){
      Rcpp::XPtr<ModelNNGP> newptr(xp);
      ptr = newptr;
    } else if(covtype == 3 && lptype == 1){
      Rcpp::XPtr<ModelHSGP> newptr(xp);
      ptr = newptr;
    } else if(covtype == 1 && lptype == 2){
      Rcpp::XPtr<ModelARRegion> newptr(xp);
      ptr = newptr;
    } else if(covtype == 2 && lptype == 2){
      Rcpp::XPtr<ModelNNGPRegion> newptr(xp);
      ptr = newptr;
    } else if(covtype == 3 && lptype == 2){
      Rcpp::XPtr<ModelHSGPRegion> newptr(xp);
      ptr = newptr;
    } else if(covtype == 1 && lptype == 3){
      Rcpp::XPtr<ModelARRegionG> newptr(xp);
      ptr = newptr;
    } else if(covtype == 2 && lptype == 3){
      Rcpp::XPtr<ModelNNGPRegionG> newptr(xp);
      ptr = newptr;
    } else if(covtype == 3 && lptype == 3){
      Rcpp::XPtr<ModelHSGPRegionG> newptr(xp);
      ptr = newptr;
    }
  }
};

using returns = std::variant<int, double, Eigen::VectorXd, Eigen::MatrixXd, std::vector<double>, std::vector<std::string>, vector_matrix >;
 