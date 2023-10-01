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
  int Q; // number of random effects, normally N=Q
  vector[N] Xb;
  matrix[N,Q] ZL;
  array[N] int y;
}
parameters {
  array[Q] real gamma;
}
model {
  int grainsize = 1;
  target += reduce_sum(partial_sum1_lpdf,gamma,grainsize);
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,Xb + ZL*to_vector(gamma));
}

