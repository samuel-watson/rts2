#pragma once

#include <glmmr/linearpredictor.hpp>
#include "regiondata.h"

namespace rts {

using namespace Eigen;
using namespace glmmr;


class regionLinearPredictor {
  public:
    rts::RegionData&  region;
    LinearPredictor   region_predictor;
    LinearPredictor   grid_predictor;
    MatrixXd*         u = nullptr;
    dblvec            parameters;
    calculator&       calc;
    
    regionLinearPredictor(glmmr::Formula& form_region,glmmr::Formula& form_grid,const Eigen::ArrayXXd &data_region,
      const Eigen::ArrayXXd &data_grid,const strvec& colnames_region,const strvec& colnames_grid,rts::RegionData& region_);    
    regionLinearPredictor(const rts::regionLinearPredictor& linpred);
    
    void        update_parameters(const dblvec& parameters_);
    void        update_parameters(const Eigen::ArrayXd& parameters_);
    int         P();
    int         n();
    strvec      colnames();
    VectorXd    xb();
    ArrayXXd    xb_region(const MatrixXd& u);
    MatrixXd    X();
    strvec      parameter_names();
    VectorXd    parameter_vector();
    bool        any_nonlinear();
    void        update_u(MatrixXd* u_);
    VectorXd    predict_xb(const ArrayXXd& newdata_,
                        const ArrayXd& newoffset_);
    
};

}



inline rts::regionLinearPredictor::regionLinearPredictor(
      glmmr::Formula& form_region,
      glmmr::Formula& form_grid,
      const Eigen::ArrayXXd &data_region,
      const Eigen::ArrayXXd &data_grid,
      const strvec& colnames_region,
      const strvec& colnames_grid,
      rts::RegionData& region_
    ) : region(region_), region_predictor(form_region,data_region,colnames_region),
        grid_predictor(form_grid,data_grid,colnames_grid), 
        parameters(region_predictor.P() + grid_predictor.P(), 0.0),
        calc(region_predictor.calc) {
  if(calc.any_nonlinear)throw std::runtime_error("Nonlinear functional forms not yet compatible with aggregated data models");
};
    
inline rts::regionLinearPredictor::regionLinearPredictor(const rts::regionLinearPredictor& linpred) : region(linpred.region), 
      region_predictor(linpred.region_predictor),
      grid_predictor(linpred.grid_predictor), 
      parameters(linpred.parameters),
      calc(region_predictor.calc) {
  if(calc.any_nonlinear)throw std::runtime_error("Nonlinear functional forms not yet compatible with aggregated data models");
};
    

inline void rts::regionLinearPredictor::update_u(MatrixXd* u_)
{
  u = u_;
}

inline void rts::regionLinearPredictor::update_parameters(const dblvec& parameters_)
{
  int dblP = region_predictor.P() + grid_predictor.P();
  dblvec par(dblP);
  if(dblP != parameters_.size()){
    Rcpp::warning("Supplied parameter vector not equal to number of parameters, generating random values\n");
    std::random_device rd{};
    std::mt19937 gen{ rd() };
    std::normal_distribution d{ 0.0, 0.5 };
    auto random_norm = [&d, &gen] { return d(gen); };
    for (int j = 0; j < par.size(); j++) par[j] = random_norm();
    for(const auto &i: par)Rcpp::Rcout << " " << i;
  } else {
    par = parameters_;
  }
  
  dblvec r_beta(region_predictor.P());
  dblvec g_beta(grid_predictor.P());
  for(int i = 0; i < region_predictor.P(); i++) r_beta[i] = par[i];
  for(int i = 0; i < grid_predictor.P(); i++) g_beta[i] = par[i+region_predictor.P()];
  region_predictor.update_parameters(r_beta);
  grid_predictor.update_parameters(g_beta);
  parameters = par;
}

inline void rts::regionLinearPredictor::update_parameters(const Eigen::ArrayXd& parameters_)
{
  dblvec new_parameters(parameters_.data(),parameters_.data()+parameters_.size());
  update_parameters(new_parameters);
};

inline int rts::regionLinearPredictor::P(){
  return region_predictor.P() + grid_predictor.P();
}

inline int rts::regionLinearPredictor::n(){
  return region_predictor.n();
}

inline strvec rts::regionLinearPredictor::colnames()
{
  strvec cnames = region_predictor.colnames();
  strvec g_cnames = grid_predictor.colnames();
  cnames.insert(cnames.end(),g_cnames.begin(),g_cnames.end());
  return cnames;
}

inline VectorXd rts::regionLinearPredictor::xb(){
  if(u!=nullptr){
    ArrayXXd xbarr = rts::regionLinearPredictor::xb_region(*u).array().exp();
    VectorXd xbvec = xbarr.rowwise().mean().log().matrix();
    return xbvec;
  } else {
    MatrixXd uzero = MatrixXd::Zero(grid_predictor.n(),1);
    ArrayXXd xbarr = rts::regionLinearPredictor::xb_region(uzero).array().exp();
    VectorXd xbvec = xbarr.rowwise().mean().log().matrix();
    return xbvec;
  }
}

inline ArrayXXd rts::regionLinearPredictor::xb_region(const MatrixXd& u){
  MatrixXd xbg = u;
  xbg.colwise() += grid_predictor.xb();
  MatrixXd xbr = region.grid_to_region(xbg);
  xbr = xbr.array().log().matrix();
  xbr.colwise() += region_predictor.xb();
  return xbr;
}

inline MatrixXd rts::regionLinearPredictor::X()
{
  MatrixXd Xg = grid_predictor.X();
  MatrixXd Xrg = region.grid_to_region(Xg);
  
  // MatrixXd xbg = MatrixXd::Zero(grid_predictor.n(),1);
  // if(u != nullptr){
  //   xbg.conservativeResize(NoChange,u->cols());
  //   xbg = *u;
  // } 
  // xbg.colwise() += grid_predictor.xb();
  // ArrayXXd gmu = xbg.array().exp().rowwise().mean();
  // MatrixXd rmu = region.grid_to_region(gmu.matrix());
  // for(int i = 0; i < Xg.rows(); i++){
  //   Xg.row(i) *= gmu(i,0);
  // }
  // MatrixXd Xrg = region.grid_to_region(Xg);
  // for(int i = 0; i < Xrg.rows(); i++){
  //   Xg.row(i) *= 1/rmu(i,0);
  // }
  MatrixXd Xr(region_predictor.n(),region_predictor.P() + Xg.cols());
  Xr.block(0,0,region_predictor.n(),region_predictor.P()) = region_predictor.X();
  Xr.block(0,region_predictor.P(),region_predictor.n(), Xg.cols()) = Xrg;
  return Xr;
}

inline strvec rts::regionLinearPredictor::parameter_names(){
  strvec pnames = region_predictor.parameter_names();
  strvec g_pnames = grid_predictor.parameter_names();
  pnames.insert(pnames.end(),g_pnames.begin(),g_pnames.end());
  return pnames;
}

inline VectorXd rts::regionLinearPredictor::parameter_vector(){
  VectorXd cnames = region_predictor.parameter_vector();
  VectorXd g_cnames = grid_predictor.parameter_vector();
  VectorXd pvec(cnames.size() + g_cnames.size());
  pvec.head(cnames.size()) = cnames;
  pvec.tail(g_cnames.size()) = g_cnames;
  return pvec;
}

inline bool rts::regionLinearPredictor::any_nonlinear(){
  return true;
}

// currently sets the random effects for the new prediction to zero, needs to sample from the random effects,
// or the random effects should be passed as an argument.
inline VectorXd rts::regionLinearPredictor::predict_xb(const ArrayXXd& newdata_,
                                                       const ArrayXd& newoffset_){
  // rts::regionLinearPredictor newlinpred(region_predictor.form,
  //                                       grid_predictor.form,
  //                                       newdata_region_,
  //                                       newdata_grid_,
  //                                       region_predictor.colnames(),
  //                                       grid_predictor.colnames(),
  //                                       region);
  // newlinpred.update_parameters(parameters);
  // 
  // MatrixXd xbg(grid_predictor.n(),1);
  // xbg.col(0) = grid_predictor.xb();
  // MatrixXd xbr = region.grid_to_region(xbg);
  // xbr = xbr.array().log().matrix();
  // xbr.colwise() += region_predictor.xb();
  // 
  // VectorXd xb = xbr.col(0) + newoffset_.matrix();
  // return xb;
  return region_predictor.predict_xb(newdata_,newoffset_);
}
