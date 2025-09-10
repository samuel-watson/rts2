// functions {
//    vector zero_sum_constrain(vector y) {
//     int N = num_elements(y);
//     vector[N + 1] z = zeros_vector(N + 1);
//     real sum_w = 0;
//     for (ii in 1:N) {
//       int i = N - ii + 1; 
//       real n = i;
//       real w = y[i] * inv_sqrt(n * (n + 1));
//       sum_w += w;
//       z[i] += sum_w;     
//       z[i + 1] -= w * n;    
//     }
//     return z;
//   }
// }
data {
  int N; // sample size
  int nT;
  int Q; // number of random effects, normally N=Q
  vector[N*nT] Xb;
  matrix[N,Q] ZL;
  array[N*nT] int y;
  real rho;
  matrix[nT,nT] ar_chol;
  real constr_zero;
}
parameters {
  //vector[Q*nT - 1] gamma_raw;
  matrix[Q,nT] gamma;
}
transformed parameters {
  //matrix[Q,nT] gamma = to_matrix(zero_sum_constrain(gamma_raw),Q,nT); // alternative sum to zero constraint + fn above
  vector[N*nT] zu = to_vector(ZL*gamma*ar_chol);
}
model {
  to_vector(gamma) ~ std_normal();
  sum(gamma) ~ normal(0, constr_zero*Q*nT); // soft sum to zero constraint
  y ~ poisson_log(Xb + zu);
}
