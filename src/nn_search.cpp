#include <Rcpp.h>
#include <queue>

using namespace Rcpp;
using namespace std;

IntegerVector top_i_pq(NumericVector v, int n) {
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;
  
  for (int i = 0; i != v.size(); ++i) {
    if (pq.size() < n)
      pq.push(Elt(v[i], i));
    else {
      Elt elt = Elt(v[i], i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }
  
  result.reserve(pq.size());
  while (!pq.empty()) {
    result.push_back(pq.top().second + 1);
    pq.pop();
  }
  
  return wrap(result);
}

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
      NN.column(i) = top_i_pq(dist,M);
    } else {
      for(int j = 0; j<i; j++){
        NN(j,i) = j;
      }
    }
  }
  return NN;
}