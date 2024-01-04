functions {
  array[] real sparse_mult(array[] int ai, array[] int ap, array[] real ax, vector x){
    int n = size(ap) - 1;
    array[n] real result = to_array_1d(rep_vector(0,n));
    for(i in 1:n){
      for(j in ap[i]:(ap[i+1]-1)){
        result[i] += x[ai[j]] * ax[j];
      }
    }
    return result;
  }
  real partial_sum1_lpdf(array[] real y, int start, int end){
    return std_normal_lpdf(y[start:end]);
  }
  real partial_sum2_lpmf(array[] int y,int start, int end, vector mu){
    return poisson_lpmf(y[start:end]|mu[start:end]);
  }
}
data {
  int N; // number of cells
  int nT;
  int Q; // number of random effects, normally N=Q
  int nRegion; // number of regions
  int ssize;
  matrix[N,Q] ZL;
  array[nRegion*nT + 1] int  Ap;
  array[ssize] int  Ai;
  array[ssize] real  Ax;
  array[nRegion*nT] int y;
  vector[nRegion*nT] Xb;
  real rho;
  matrix[nT,nT] ar_chol;
}
parameters {
   matrix[Q,nT] gamma;
}
transformed parameters {
  vector[N*nT] zu = to_vector(ZL*gamma*ar_chol);
}
model {
  int grainsize = 1;
  array[nRegion*nT] real u = sparse_mult(Ai,Ap,Ax,exp(zu)); 
  target += reduce_sum(partial_sum1_lpdf,to_array_1d(gamma),grainsize);
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,Xb .* to_vector(u));
}

