#include <Rcpp.h>
#include <queue>
#include "../inst/include/rtsheader.h"

using namespace Rcpp;
using namespace std;

//' Generate matrix of nearest neighbours
//' 
//' Generates a M by n matrix of nearest neighbours. For observations with index less than M,
//' only i-1 nearest neighbours are generated.
//' 
//' @param x Matrix with x and y coordinates of the observations
//' @param M Number of nearest neighbours
//' @return A M by n matrix of nearest neighbours. Note the indexing of the observations returned
//' by this function starts at zero.
// [[Rcpp::export]]
IntegerMatrix genNN(const NumericMatrix &x,
                          int M){
  int n = x.rows();
  IntegerMatrix NN(M,n);
  std::fill(NN.begin(),NN.end(),n+1);
  for(int i=1; i<n; i++){
    NumericVector dist(i);
    if(i > M){
      for(int j=0; j<(i-1); j++){
        dist(j) = sqrt((x(i,0) - x(j,0))*(x(i,0) - x(j,0))+(x(i,1) - x(j,1))*(x(i,1) - x(j,1)));
      }
      NN.column(i) = rts::top_i_pq(dist,M);
    } else {
      for(int j = 0; j<i; j++){
        NN(j,i) = j;
      }
    }
  }
  return NN;
}