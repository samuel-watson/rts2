functions {
  real calc_lambda(real p, real xg, vector w, vector f){
    return p*exp(xg)*(w' * exp(f));
  }
  real poisson_log_block_lpmf(array[] int y, array[] int i, int nT, 
                              int nR, int nS, vector pop, matrix X, 
                vector w, vector f, array[] int cell_id, array[] int n_cell,
                vector gamma){ 
    // calculate the IDs of the observation
    int n = size(y);
    vector[n] lambda;
    for(j in 1:n){
      int r = i[j]%nR;
      if(r == 0)r = nR;
      int t = to_int(ceil(i[j]*1.0/nR));
      int lsize = n_cell[r+1] - n_cell[r];
      array[lsize] int idx;
      for(l in 1:lsize)idx[l] = cell_id[n_cell[r]+l-1]+(t-1)*nS;
      lambda[j] = calc_lambda(pop[r+(t-1)*nR],
                         X[r+(t-1)*nR,]*gamma,
                         w[(n_cell[r]):(n_cell[r+1]-1)],
                         f[idx]);
    }
    
    return poisson_lpmf(y | lambda);                  
                              
  }
  real partial_sum1_lpdf(array[] real y, int start, int end){
    return std_normal_lpdf(y[start:end]);
  }
  // real partial_sum2_lpmf(array[] int y,int start, int end, int nT,
  //                       int nR, vector pop, matrix X, 
  //                       vector w, vector f, array[] int cell_id, array[] int n_cell,
  //                       vector gamma){
  //   array[end-start+1] int idx = linspaced_int_array(end-start+1,start,end);
  //   return poisson_log_block_lpmf(y[start:end]|idx, nT, nR, pop,
  //                                 X, w, f, cell_id, n_cell, gamma);
  // }
}
data {
  // define the model and problem
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=0> Q_g; //number of covariates
  array[D] real L; // boundary condition
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  int n_region; // number of regions
  int n_Q; // number of intersections
  array[n_region+1] int<lower=1> n_cell; //number of cells intersecting region  
  array[n_Q] int<lower=1> cell_id; // IDs of the cells intersecting the region
  vector[n_Q] q_weights; // proportionate weights
  
  // outcomes data
  array[n_region*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  array[M_nD,D] int indices; //indices
  vector[n_region*nT] popdens; //population density
  matrix[n_region*nT,Q] X;
  matrix[Nsample*nT,Q_g == 0 ? 1 : Q_g] X_g;
  
  // distribution parameters
  array[2] real prior_lscale;
  array[2] real prior_var;
  array[Q] real prior_linpred_mean;
  array[Q] real prior_linpred_sd;
  int mod;
  int<lower = 0, upper = 1> known_cov;
  real<lower=0> sigma_data; 
  real<lower=0> phi_data; 
}
transformed data {
  
}

parameters {
  array[known_cov ? 0 : 1] real<lower=1e-05> phi_param; //length scale
  array[known_cov ? 0 : 1] real<lower=1e-05> sigma_param;
  vector[Q] gamma;
  vector[Q_g] gamma_g;
  real<lower=-1,upper=1> ar;
  array[Nsample*nT] real f_raw;
}

transformed parameters{
  matrix[Nsample,Nsample] L;
  vector[Nsample*nT] f;
  real<lower=1e-05> sigma;
  real<lower=1e-05> phi;
  if(known_cov){
    sigma = sigma_data;
    phi = phi_data;
  } else {
    sigma = sigma_param[1];
    phi = phi_param[1];
  }

  L = genChol(Nsample, sigma, phi, x_grid, mod);
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        f[1:Nsample] = (1/(1-ar^2))*L*to_vector(f_raw[1:Nsample]);
      } else {
        f[(Nsample*(t-1)+1):(t*Nsample)] = ar*L*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + L*to_vector(f_raw[(Nsample*(t-1)+1):(t*Nsample)]);
      }
    } else {
      f = L*to_vector(f_raw);
    }

  }
  
  if(Q_g > 0){
    f += X_g * gamma_g;
  }

}
model{
  vector[n_region*nT] lambda_r = rep_vector(0,n_region*nT);
  if(!known_cov){
    phi_param ~ normal(prior_lscale[1],prior_lscale[2]);
    sigma_param ~ normal(prior_var[1],prior_var[2]);
  }
  ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  if(Q_g > 0){
    gamma_g ~ normal(0,2);
  }
  array[n_region*nT] int idx = linspaced_int_array(n_region*nT,1,n_region*nT);
  target += reduce_sum(partial_sum1_lpdf,f_raw,grainsize);
  y ~ poisson_log_block(idx, nT, n_region, Nsample, popdens, X, q_weights, f, cell_id, n_cell, gamma);
}

generated quantities{
  vector[Nsample*nT] y_grid_predict;
  vector[n_region*nT] region_predict;

  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(f[i]);
  }
  
  region_predict = rep_vector(0,n_region*nT);
  for(r in 1:n_region){
    for(t in 1:nT){
      for(l in 1:(n_cell[r+1]-n_cell[r])){
        region_predict[r+(t-1)*n_region] += popdens[r+(t-1)*n_region]*exp(X[r+(t-1)*n_region,]*gamma)*
          q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*Nsample]);
      }
    }
  }
  
}
