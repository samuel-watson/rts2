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
  to_array_1d(gamma) ~ std_normal();
  y ~ poisson_log(Xb + zu);
}