#pragma once

#include "rtsmodelbits.h"

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