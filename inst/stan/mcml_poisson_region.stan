data {
  int N; // number of cells
  int Q; // number of random effects, normally N=Q
  int nRegion; // number of regions
  matrix[N,Q] ZL;
  matrix[nRegion,N] P;
  array[nRegion] int y;
}
parameters {
  vector[Q] gamma;
}
model {
  vector[N] u = exp(ZL*to_vector(gamma));
  gamma ~ std_normal();
  y ~ poisson(P*u);
}

