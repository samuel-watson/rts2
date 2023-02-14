#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::export]]
double max_dist(const Eigen::ArrayXXd &x){
  // this is a brute force algorithm for max distance
  // it can be improved by finding convex hull and then using rotating calipers method
  // but I haven't had the time to implement that!
  int n = x.rows();
  double maxdist = 0;
  double dist = 0;
  for(int i = 1; i < n; i++){
    for(int j = 0; j<(i-1); j++){
      dist = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
      if(dist > maxdist) maxdist = dist;
    }
  }
  return maxdist;
}

// [[Rcpp::export]]
Eigen::ArrayXXd semivariogram(const Eigen::ArrayXXd &x,
                              const Eigen::ArrayXd &offs,
                              const Eigen::ArrayXd &y,
                              int nbins){
  double maxd = max_dist(x);
  Eigen::ArrayXd denoms = Eigen::ArrayXd::Zero(nbins);
  Eigen::ArrayXd sums = Eigen::ArrayXd::Zero(nbins);
  double binw = maxd/nbins;
  Eigen::ArrayXd z = offs.inverse();
  z *= y;
  double dist;
  int n = x.rows();
  int binc;
  for(int i = 1; i < n; i++){
    for(int j = 0; j<(i-1); j++){
      dist = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
      binc = static_cast<int>(std::floor(dist/binw));
      denoms(binc) += offs(i)*offs(j);
      sums(binc) += offs(i)*offs(j)*(z(i)-z(j))*(z(i)-z(j));
    }
  }
  denoms *= 2;
  Eigen::ArrayXXd result(nbins,2);
  for(int i=0; i<nbins; i++)result(i,0) = i*binw + binw/2;
  result.col(1) = denoms.inverse()*sums;
  return result;
}