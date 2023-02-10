data {
  int N; // sample size
  vector[N] Xb;
  matrix[N,N] Z;
  int y[N];
}
parameters {
  vector[N] gamma;
}
model {
  gamma ~ std_normal();
  y ~ poisson_log(Xb + Z*gamma);
}

