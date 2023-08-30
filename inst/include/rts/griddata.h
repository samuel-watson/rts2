#pragma once

#include <glmmr/general.h>

namespace rts {

using namespace Eigen;

class griddata {
public:
  int T = 1; // number of time periods
  int N = 1; // number of cells
  ArrayXXd X; // centroids
  ArrayXXi NN = ArrayXXi::Constant(1,1,1);
  int m = 10;
  
  griddata(const ArrayXXd& X_, int T_) : X(X_), T(T_) {};
  griddata(const ArrayXXd& X_, int T_, int m) : X(X_), T(T_) {genNN(m);};
  griddata(const rts::griddata& g) : X(g.X), N(g.N), T(g.T) {};

  void genNN(int m);
  void setup(const ArrayXXd& X_, int T_);
  void setup(const ArrayXXd& X_, int T_, int m);
};

}

inline void rts::griddata::genNN(int M){
  int n = X.rows();
  m = M;
  NN.conservativeResize(M,n);
  NN = ArrayXXi::Constant(M,n,n+1);
  for(int i=1; i<n; i++){
    ArrayXd dist(i);
    if(i > M){
      for(int j=0; j<(i-1); j++){
        dist(j) = sqrt((X(i,0) - X(j,0))*(X(i,0) - X(j,0))+(X(i,1) - X(j,1))*(X(i,1) - X(j,1)));
      }
      NN.col(i) = rts::top_i_pq(dist,M);
    } else {
      for(int j = 0; j<i; j++){
        NN(j,i) = j;
      }
    }
  }
}

inline void rts::griddata::setup(const ArrayXXd& X_, int T_){
  X.conservativeResize(X_.rows(),X_.cols());
  X = X_;
  T = T_;
  N = X_.rows();
}

inline void rts::griddata::setup(const ArrayXXd& X_, int T_, int m){
  X.conservativeResize(X_.rows(),X_.cols());
  X = X_;
  T = T_;
  N = X_.rows();
  genNN(m);
}
