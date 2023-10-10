functions {
  real partial_sum1_lpdf(array[] real y, int start, int end){
    return std_normal_lpdf(y[start:end]);
  }
  real partial_sum2_lpmf(array[] int y,int start, int end, vector mu){
    return poisson_lpmf(y[start:end]|mu[start:end]);
  }
}
data {
  int N; // number of cells
  int Q; // number of random effects, normally N=Q
  int nRegion; // number of regions
  matrix[N,Q] ZL;
  matrix[nRegion,N] P;
  array[nRegion] int y;
}
parameters {
  array[Q] real gamma;
}
transformed parameters {
  
}
model {
  int grainsize = 1;
  vector[N] u = exp(ZL*to_vector(gamma));
  target += reduce_sum(partial_sum1_lpdf,gamma,grainsize);
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,P*u);
}

