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
    nRegion(n_cell.size()-1) { setup_design_matrices(); };
  RegionData(const rts::RegionData& region) : n_cell(region.n_cell), cell_id(region.cell_id),
    q_weights(region.q_weights), gridT(region.gridT), gridN(region.gridN), nRegion(region.nRegion),
    region_to_intersection(region.region_to_intersection), grid_to_intersection(region.grid_to_intersection),
    grid_region(region.grid_region) {};
  MatrixXd grid_to_region(const MatrixXd& u);
  sparse region_design_matrix();
  sparse grid_design_matrix();
  sparse grid_to_region_matrix();

protected:
  void setup_design_matrices();
  sparse region_to_intersection;
  sparse grid_to_intersection;
  sparse grid_region;
};

}

inline void rts::RegionData::setup_design_matrices(){
  sparse A(q_weights.size(),nRegion,true);
  for(int r = 0; r < nRegion; r++){
    for(int l = n_cell(r); l < n_cell(r+1); l++){
      A.insert(l,r,1);
    }
  }
  region_to_intersection = A;
  
  sparse B(q_weights.size(),gridN,true);
  for(int i = 0; i < q_weights.size(); i++){
    B.insert(i,cell_id(i),1);
  }
  grid_to_intersection = B;
  
  sparse C(nRegion*gridT,gridN*gridT,true);
  int nInter, r, t, l, idx1;
  for(r = 0; r < nRegion; r++){
    for(l = n_cell(r); l < n_cell(r+1); l++){
      for(t = 0; t < gridT; t++){
        C.insert(r + nRegion*t,cell_id(l) + t*gridN,q_weights(l));
      }
    }
  }
  grid_region = C;
  
}

inline sparse rts::RegionData::region_design_matrix(){
  // sparse A(q_weights.size(),nRegion,true);
  // for(int r = 0; r < nRegion; r++){
  //   for(int l = n_cell(r); l < n_cell(r+1); l++){
  //     A.insert(l,r,1);
  //   }
  // }
  // return A;
  return region_to_intersection;
}

inline sparse rts::RegionData::grid_design_matrix(){
  // sparse A(q_weights.size(),gridN,true);
  // for(int i = 0; i < q_weights.size(); i++){
  //   A.insert(i,cell_id(i),1);
  // }
  // return A;
  return grid_to_intersection;
}

inline sparse rts::RegionData::grid_to_region_matrix(){
  // sparse A(nRegion*gridT,gridN*gridT,true);
  // int nInter, r, t, l, idx1;
  // for(r = 0; r < nRegion; r++){
  //   for(l = n_cell(r); l < n_cell(r+1); l++){
  //     for(t = 0; t < gridT; t++){
  //       A.insert(r + nRegion*t,cell_id(l) + t*gridN,q_weights(l));
  //     }
  //   }
  // }
  // return A;
  return grid_region;
}

inline MatrixXd rts::RegionData::grid_to_region(const MatrixXd& u){
  MatrixXd regionu = MatrixXd::Zero(nRegion*gridT, u.cols());
  int nInter, r, t, l, j, idx1;
  if(n_cell(0)!=0)Rcpp::stop("Indexing does not start from zero");
  
  for(r = 0; r < nRegion; r++){
    for(t = 0; t < gridT; t++){
      for(l = n_cell(r); l < n_cell(r+1); l++){
        for(j = 0; j < u.cols(); j++){
          regionu(r + nRegion*t, j) += q_weights(l)*exp(u(cell_id(l) + t*gridN,j));
        }
      }
    }
  }
  
  return regionu;
}