# pragma once

#include <glmmr/general.h>
#include <glmmr/openmpheader.h>

namespace rts {

using namespace Eigen;

class RegionData {
public: 
  ArrayXi n_cell;
  ArrayXi cell_id;
  ArrayXd q_weights;
  int gridT;
  int gridN;
  int nRegion;
  RegionData(const ArrayXi &n_cell_,
             const ArrayXi &cell_id_,
             const ArrayXd &q_weights_,
             int N_, int T_) : 
    n_cell(n_cell_), cell_id(cell_id_), q_weights(q_weights_), 
    gridT(T_), gridN(N_), 
    nRegion(n_cell.size()-1) {};
  RegionData(const rts::RegionData& region) : n_cell(region.n_cell), cell_id(region.cell_id),
    q_weights(region.q_weights), gridT(region.gridT), gridN(region.gridN), nRegion(region.nRegion) {};
  MatrixXd grid_to_region(const MatrixXd& u);
  MatrixXd region_design_matrix();
  MatrixXd grid_design_matrix();
};

}

inline MatrixXd rts::RegionData::region_design_matrix(){
  MatrixXd A = MatrixXd::Zero(q_weights.size(),nRegion*gridT);
  int nInter, r, t, l, idx1;
  for(r = 0; r < nRegion; r++){
    nInter = n_cell(r+1) - n_cell(r);
    for(l = 0; l < nInter; l++){
      idx1 = n_cell(r) + l;
      for(t = 0; t < gridT; t++){
        A(idx1,r+nRegion*t) = 1;
      }
    }
  }
  return A;
}

inline MatrixXd rts::RegionData::grid_design_matrix(){
  MatrixXd A = MatrixXd::Zero(q_weights.size(),gridN*gridT);
  for(int i = 0; i < q_weights.size(); i++){
    for(int t = 0; t < gridT; t++){
      A(i,cell_id(i)+t*gridN) = 1;
    }
  }
  
  return A;
}

inline MatrixXd rts::RegionData::grid_to_region(const MatrixXd& u){
  MatrixXd regionu = MatrixXd::Zero(nRegion*gridT, u.cols());
  int nInter, r, t, l, j, idx1;
  if(n_cell(0)!=0)Rcpp::stop("Indexing does not start from zero");
  
  for(r = 0; r < nRegion; r++){
    for(t = 0; t < gridT; t++){
      nInter = n_cell(r+1) - n_cell(r);
      for(l = 0; l < nInter; l++){
        idx1 = n_cell(r) + l;
        for(j = 0; j < u.cols(); j++){
          regionu(r + nRegion*t, j) += q_weights(idx1)*exp(u(cell_id(idx1) + t*gridN,j));
        }
      }
    }
  }
  
  return regionu;
}