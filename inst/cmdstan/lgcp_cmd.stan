functions {
  matrix genChol(int n, real alpha, real theta, matrix x, int mod){
    matrix[n,n] L = rep_matrix(0,n,n);
    real s;
    real dist;
    L[1,1] = alpha;
    
    for(i in 2:n){
      dist = sqrt((x[1,1] - x[i,1]) * (x[1,1] - x[i,1]) + (x[1,2] - x[i,2]) * (x[1,2] - x[i,2]));
      if(mod == 0){
        L[i,1] = alpha * exp(-(dist*dist)/(theta*theta));
      }  else {
        L[i,1] = alpha * exp(-dist/theta);
      }
    }
    
    for(j in 2:n){
      s = 0;
      for(k in 1:(j-1)){
        s = s + L[j,k]*L[j,k];
      }
      L[j,j] = sqrt(alpha - s);
      for(i in (j+1):n){
        dist = sqrt((x[j,1] - x[i,1]) * (x[j,1] - x[i,1]) + (x[j,2] - x[i,2]) * (x[j,2] - x[i,2]));
        s = 0;
        for(k in 1:(j-1)){
          s = s + L[j,k]*L[i,k];
        }
        if(mod == 0){
          L[i,j] = 1/L[j,j] * (alpha *  alpha * exp(-(dist*dist)/(theta*theta)) - s);
        }  else {
          L[i,j] = 1/L[j,j] * (alpha * alpha * exp(-dist/theta) - s);
        }
      }
    }
    return L;
  }
  real partial_sum2_lpmf(array[] int y,int start, int end, vector mu){
    return poisson_log_lpmf(y[start:end]|mu[start:end]);
  }
  real partial_sum1_lpdf(array[] real y, int start, int end){
    return std_normal_lpdf(y[start:end]);
  }
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  array[Nsample*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[Nsample*nT] popdens; //population density
  matrix[Nsample*nT,Q] X;
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
  vector[Nsample*nT] logpopdens = log(popdens);
}

parameters {
  array[known_cov ? 0 : 1] real<lower=1e-05> phi_param; //length scale
  array[known_cov ? 0 : 1] real<lower=1e-05> sigma_param;
  vector[Q] gamma;
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
}
model{
  if(!known_cov){
    to_vector(phi_param) ~ normal(prior_lscale[1],prior_lscale[2]);
    sigma_param ~ normal(prior_var[1],prior_var[2]);
  }
  ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  int grainsize = 1;
  
  target += reduce_sum(partial_sum1_lpdf,f_raw,grainsize);
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,X*gamma+logpopdens+f);

}

generated quantities{
  vector[Nsample*nT] y_grid_predict;

  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(X[i,]*gamma+logpopdens[i]+f[i]);
  }
}
