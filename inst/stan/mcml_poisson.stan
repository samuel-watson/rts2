data {
  int N; // sample size
  int Q; // number of random effects, normally N=Q
  vector[N] Xb;
  matrix[N,Q] ZL;
  int y[N];
}
parameters {
  vector[Q] gamma;
}
model {
  gamma ~ std_normal();
  y ~ poisson_log(Xb + ZL*gamma);
}

