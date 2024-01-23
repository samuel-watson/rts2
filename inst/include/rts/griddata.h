#pragma once

#include <glmmr/general.h>

namespace rts {

using namespace Eigen;

class griddata {
public:
  ArrayXXd  X; // centroids
  int       T; // number of time periods
  int       N; // number of cells
  ArrayXXi  NN = ArrayXXi::Constant(1,1,1);
  int       m = 10;
  
  griddata(const ArrayXXd& X_, int T_) : X(X_), T(T_), N(X_.rows()) {};
  griddata(const ArrayXXd& X_, int T_, int M) : X(X_), T(T_), N(X_.rows()) {genNN(M);};
  griddata(const rts::griddata& g) : X(g.X), T(g.T), N(g.N) {};

  void      genNN(int M);
  void      setup(const ArrayXXd& X_, int T_);
  void      setup(const ArrayXXd& X_, int T_, int M);
  ArrayXi   top_i_pq(const ArrayXd& v, int n);
};
}


inline void rts::griddata::genNN(int M){
  int n = X.rows();
  m = M;
  NN.conservativeResize(M,n);
  NN = ArrayXXi::Constant(M,n,n);
  for(int i=1; i<n; i++){
    if(i > M){
      ArrayXd dist(i);
      for(int j=0; j<i; j++){
        dist(j) = (X(i,0) - X(j,0))*(X(i,0) - X(j,0))+(X(i,1) - X(j,1))*(X(i,1) - X(j,1));
      }
      NN.col(i) = top_i_pq(dist,M);
    } else {
      for(int j = 0; j<i; j++){
        NN(j,i) = j;
      }
    }
  }
}

inline ArrayXi rts::griddata::top_i_pq(const ArrayXd& v, int n) {
  typedef std::pair<double, int> Elt;
  
  struct ComparePair {
    bool operator()(const Elt& elt1, const Elt& elt2) const {
      return elt1.first < elt2.first;
    }
  };
  
  std::priority_queue< Elt, std::vector<Elt>, ComparePair > pq;
  std::vector<int> result;
  
  for (int i = 0; i < v.size(); i++) {
    if (pq.size() < n)
      pq.push(Elt(v(i), i));
    else {
      Elt elt = Elt(v(i), i);
      if (pq.top().first > elt.first) {
        pq.pop();
        pq.push(elt);
      }
    }
  }
  
  ArrayXi res(pq.size());
  int iter = 0;
  while (!pq.empty()) {
    res(iter) = pq.top().second;
    pq.pop();
    iter++;
  }
  
  return res;
}

inline void rts::griddata::setup(const ArrayXXd& X_, int T_){
  X.conservativeResize(X_.rows(),X_.cols());
  X = X_;
  T = T_;
  N = X_.rows();
}

inline void rts::griddata::setup(const ArrayXXd& X_, int T_, int M){
  X.conservativeResize(X_.rows(),X_.cols());
  X = X_;
  T = T_;
  N = X_.rows();
  genNN(M);
}
