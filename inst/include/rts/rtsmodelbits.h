#pragma once

#include <glmmr/general.h>
#include <glmmr/modelbits.hpp>
#include <glmmr/randomeffects.hpp>
#include <glmmr/modelmatrix.hpp>
#include <glmmr/modelmcmc.hpp>
#include <glmmr/family.hpp>
#include <glmmr/modelextradata.hpp>
#include <glmmr/calculator.hpp>
#include <glmmr/formula.hpp>
#include "ar1covariance.h"
#include "nngpcovariance.h"
#include "hsgpcovariance.h"
#include "regionlinearpredictor.h"
#include "regiondata.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;

class rtsModelBitsBase {
public:
  glmmr::Formula          formula;
  glmmr::ModelExtraData   data;
  glmmr::Family           family;
  glmmr::calculator       calc;
  glmmr::calculator       vcalc;
  bool                    weighted = false;
  
  rtsModelBitsBase(const glmmr::Formula& formula_, 
                 const glmmr::ModelExtraData& data_,
                 const glmmr::Family& family_) : formula(formula_),
    data(data_), family(family_) {};
  
  rtsModelBitsBase(const std::string& formula_, 
                   const ArrayXXd& data_) : formula(formula_),
                    data(data_.rows()), family("poisson","log") {};
  
  rtsModelBitsBase(const rts::rtsModelBitsBase& bits) : formula(bits.formula),
    data(bits.data), family(bits.family) {};
  
  virtual int       n(){ return 0; };
  virtual ArrayXd   xb(){return ArrayXd::Zero(1);};
  ~rtsModelBitsBase() = default;
};

template<typename cov, typename linpred>
class rtsModelBits : public rtsModelBitsBase {
  cov           covariance;
  linpred       linear_predictor;
  rtsModelBits(){};
  ~rtsModelBits() = default;
  int           n() override;
  ArrayXd       xb() override;
};

template<>
class rtsModelBits<rts::ar1Covariance, LinearPredictor> : public rtsModelBitsBase {
public:
  rts::ar1Covariance  covariance;
  LinearPredictor     linear_predictor;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_) : 
    rtsModelBitsBase(formula_,data_),
    covariance(formula_,data_,colnames_, 1),
    linear_predictor(formula,data_,colnames_) {};
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               int T,
               const ArrayXXd& grid_data_) : 
    rtsModelBitsBase(formula_,data_),
    covariance(formula_,grid_data_,std::vector<std::string>({"X","Y"}), T),
    linear_predictor(formula,data_,colnames_) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::ar1Covariance, LinearPredictor>& bits) :
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int       n(){return linear_predictor.n();};
  ArrayXd   xb(){return linear_predictor.xb() + data.offset;};
  
};

template<>
class rtsModelBits<rts::nngpCovariance, LinearPredictor> : public rtsModelBitsBase {
public:
  rts::nngpCovariance   covariance;
  LinearPredictor       linear_predictor;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               int m,
               const rts::griddata& grid_) : 
    rtsModelBitsBase(formula_,data_),
    covariance(formula_,data_,colnames_, 1, m, grid_),
    linear_predictor(formula,data_,colnames_) {};
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               int T, int m,
               const rts::griddata& grid_,
               const ArrayXXd& grid_data_) : 
    rtsModelBitsBase(formula_,data_),
    covariance(formula_,grid_data_,std::vector<std::string>({"X","Y"}), T, m, grid_),
    linear_predictor(formula,data_,colnames_) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::nngpCovariance, LinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int       n(){return linear_predictor.n();};
  ArrayXd   xb(){return linear_predictor.xb() + data.offset;};
};

template<>
class rtsModelBits<rts::hsgpCovariance, LinearPredictor> : public rtsModelBitsBase {
public:
  rts::hsgpCovariance   covariance;
  LinearPredictor       linear_predictor;
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               int m,
               const Array2d& L) : 
    rtsModelBitsBase(formula_,data_),
    covariance(formula_,data_,colnames_, 1, m, L),
    linear_predictor(formula,data_,colnames_) {};
  
  rtsModelBits(const std::string& formula_,
               const ArrayXXd& data_,
               const strvec& colnames_,
               int T,
               int m,
               const ArrayXd& L,
               const ArrayXXd& grid_data_) : 
    rtsModelBitsBase(formula_,data_),
    covariance(formula_,grid_data_,std::vector<std::string>({"X","Y"}), T, m, L),
    linear_predictor(formula,data_,colnames_) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::hsgpCovariance, LinearPredictor>& bits) :
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int       n(){return linear_predictor.n();};
  ArrayXd   xb(){return linear_predictor.xb() + data.offset;};
  
};

// need to implement a calculator for the region model
// to enable the gradient functions.

template<>
class rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> : public rtsModelBitsBase {
public:
  glmmr::Formula                formula_grid;
  rts::ar1Covariance            covariance;
  rts::regionLinearPredictor    linear_predictor;
  
  rtsModelBits(const std::string& form_region,
               const std::string& form_grid,
               const Eigen::ArrayXXd &data_region,
               const Eigen::ArrayXXd &data_grid,
               const strvec& colnames_region,
               const strvec& colnames_grid,
               int T,
               rts::RegionData& region) : 
    rtsModelBitsBase(form_region,data_region),
    formula_grid(form_grid),
    covariance(form_grid,data_grid,colnames_grid, T),
    linear_predictor(formula,formula_grid,data_region,data_grid,colnames_region,colnames_grid,region) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    formula_grid(bits.formula_grid),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int       n(){return linear_predictor.n();};
  ArrayXd   xb(){return linear_predictor.xb() + data.offset;};
};

template<>
class rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> : public rtsModelBitsBase {
public:
  glmmr::Formula              formula_grid;
  rts::nngpCovariance         covariance;
  rts::regionLinearPredictor  linear_predictor;
  
  rtsModelBits(const std::string& form_region,
               const std::string& form_grid,
               const Eigen::ArrayXXd &data_region,
               const Eigen::ArrayXXd &data_grid,
               const strvec& colnames_region,
               const strvec& colnames_grid,
               rts::RegionData& region,
               const rts::griddata& grid_,
               int T, int m) : 
  rtsModelBitsBase(form_region,data_region),
  formula_grid(form_grid),
  covariance(form_grid,data_grid,colnames_grid, T, m, grid_),
  linear_predictor(formula,formula_grid,data_region,data_grid,colnames_region,colnames_grid,region) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    formula_grid(bits.formula_grid),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int       n(){return linear_predictor.n();};
  ArrayXd   xb(){return linear_predictor.xb() + data.offset;};
  
};

template<>
class rtsModelBits<rts::hsgpCovariance, rts::regionLinearPredictor> : public rtsModelBitsBase {
public:
  glmmr::Formula              formula_grid;
  rts::hsgpCovariance         covariance;
  rts::regionLinearPredictor  linear_predictor;
  
  rtsModelBits(const std::string& form_region,
               const std::string& form_grid,
               const Eigen::ArrayXXd &data_region,
               const Eigen::ArrayXXd &data_grid,
               const strvec& colnames_region,
               const strvec& colnames_grid,
                 int T, int m,
                 const ArrayXd& L,
                 rts::RegionData& region) : 
    rtsModelBitsBase(form_region,data_region),
    formula_grid(form_grid),
    covariance(form_grid,data_grid,colnames_grid, T, m, L),
    linear_predictor(formula,formula_grid,data_region,data_grid,colnames_region,colnames_grid,region) {};
  
  rtsModelBits(const rts::rtsModelBits<rts::hsgpCovariance, rts::regionLinearPredictor>& bits) : 
    rtsModelBitsBase(bits.formula,bits.data,bits.family),
    formula_grid(bits.formula_grid),
    covariance(bits.covariance), linear_predictor(bits.linear_predictor) {};
  
  int       n(){return linear_predictor.n();};
  ArrayXd   xb(){return linear_predictor.xb() + data.offset;};
};

}

typedef rts::rtsModelBits<rts::ar1Covariance, glmmr::LinearPredictor> BitsAR;
typedef rts::rtsModelBits<rts::nngpCovariance, glmmr::LinearPredictor> BitsNNGP;
typedef rts::rtsModelBits<rts::hsgpCovariance, glmmr::LinearPredictor> BitsHSGP;
typedef rts::rtsModelBits<rts::ar1Covariance, rts::regionLinearPredictor> BitsARRegion;
typedef rts::rtsModelBits<rts::nngpCovariance, rts::regionLinearPredictor> BitsNNGPRegion;
typedef rts::rtsModelBits<rts::hsgpCovariance, rts::regionLinearPredictor> BitsHSGPRegion;

template<>
inline VectorMatrix glmmr::RandomEffects<BitsAR>::predict_re(const ArrayXXd& newdata_,
                                                           const ArrayXd& newoffset_){
#ifdef R_BUILD
  if(model.covariance.data_.cols()!=newdata_.cols())Rcpp::stop("Different numbers of columns in new data");
  Rcpp::stop("Predicting random effects is not yet implemented in this package for LGCPs");
#endif
}

template<>
inline VectorMatrix glmmr::RandomEffects<BitsARRegion>::predict_re(const ArrayXXd& newdata_,
                                                             const ArrayXd& newoffset_){
#ifdef R_BUILD
  if(model.covariance.data_.cols()!=newdata_.cols())Rcpp::stop("Different numbers of columns in new data");
  Rcpp::stop("Predicting random effects is not yet implemented in this package for LGCPs");
#endif
}

template<>
inline VectorMatrix glmmr::RandomEffects<BitsNNGP>::predict_re(const ArrayXXd& newdata_,
                                                             const ArrayXd& newoffset_){
#ifdef R_BUILD
  if(model.covariance.data_.cols()!=newdata_.cols())Rcpp::stop("Different numbers of columns in new data");
  Rcpp::stop("Predicting random effects is not yet implemented in this package for LGCPs");
#endif
}

template<>
inline VectorMatrix glmmr::RandomEffects<BitsNNGPRegion>::predict_re(const ArrayXXd& newdata_,
                                                             const ArrayXd& newoffset_){
#ifdef R_BUILD
  if(model.covariance.data_.cols()!=newdata_.cols())Rcpp::stop("Different numbers of columns in new data");
  Rcpp::stop("Predicting random effects is not yet implemented in this package for LGCPs");
#endif
}

template<>
inline VectorMatrix glmmr::RandomEffects<BitsHSGP>::predict_re(const ArrayXXd& newdata_,
                                                             const ArrayXd& newoffset_){
#ifdef R_BUILD
  if(model.covariance.data_.cols()!=newdata_.cols())Rcpp::stop("Different numbers of columns in new data");
  Rcpp::stop("Predicting random effects is not yet implemented in this package for LGCPs");
#endif
}

template<>
inline VectorMatrix glmmr::RandomEffects<BitsHSGPRegion>::predict_re(const ArrayXXd& newdata_,
                                                             const ArrayXd& newoffset_){
#ifdef R_BUILD
  if(model.covariance.data_.cols()!=newdata_.cols())Rcpp::stop("Different numbers of columns in new data");
  Rcpp::stop("Predicting random effects is not yet implemented in this package for LGCPs");
#endif
}




