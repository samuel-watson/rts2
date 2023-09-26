data {
  int N; // sample size
  vector[N] Xb;
  matrix[N,N] ZL;
  int y[N];
}
parameters {
  vector[N] gamma;
}
model {
  gamma ~ std_normal();
  y ~ poisson_log(Xb + ZL*gamma);
}

