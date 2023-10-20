functions {
  real partial_sum1_lpdf(array[] real y, int start, int end){
    return std_normal_lpdf(y[start:end]);
  }
  real partial_sum2_lpmf(array[] int y,int start, int end, vector mu){
    return poisson_log_lpmf(y[start:end]|mu[start:end]);
  }
}
data {
  int N; // sample size
  int nT;
  int Q; // number of random effects, normally N=Q
  vector[N*nT] Xb;
  matrix[N,Q] ZL;
  array[N*nT] int y;
  real rho;
  matrix[nT,nT] ar_chol;
}
parameters {
  matrix[Q,nT] gamma;
}
transformed parameters {
  vector[Q*nT] zu = to_vector(ZL*gamma*ar_chol);
}
model {
  int grainsize = 1;
  target += reduce_sum(partial_sum1_lpdf,to_array_1d(gamma),grainsize);
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,Xb + zu);
}

