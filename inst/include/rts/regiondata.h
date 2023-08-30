# pragma once

#include <glmmr/general.h>
#include <glmmr/openmpheader.h>
#include "griddata.h"

namespace rts {

using namespace Eigen;

class RegionData {
public: 
  const ArrayXi &n_cell;
  const ArrayXi &cell_id;
  const ArrayXd &q_weights;
  const rts::griddata& grid;
  int nRegion;
  RegionData(const ArrayXi &n_cell_,
             const ArrayXi &cell_id_,
             const ArrayXd &q_weights_,
             const rts::griddata& grid_) : 
    n_cell(n_cell_), cell_id(cell_id_), q_weights(q_weights_), grid(grid_), 
    nRegion(n_cell.size()-1) {};
  RegionData(const rts::RegionData& region) : n_cell(region.n_cell), cell_id(region.cell_id),
    q_weights(region.q_weights), grid(region.grid), nRegion(region.nRegion) {};
  MatrixXd grid_to_region(const MatrixXd& u);
};

}

inline MatrixXd rts::RegionData::grid_to_region(const MatrixXd& u){
  MatrixXd regionu = MatrixXd::Zero(nRegion*grid.T, u.cols());
  
#pragma omp parallel for
  for(int r=0; r<nRegion;r++){
    for(int t=0; t<grid.T; t++){
      int nInter = n_cell(r+1)-n_cell(r);
      for(int j=0; j<u.cols(); j++){
        double accum = 0;
        for(int l=0; l<nInter; l++){
          accum += q_weights(n_cell(r)-1+l)*exp(u(cell_id(n_cell(r)-1+l) + t*grid.N,j));
        }
        regionu(r + t*nRegion,j) = accum;
      }
    }
  }
  
  return regionu;
}