functions {
  matrix genChol(int n, real alpha, real theta, real[] dists, int mod){
    matrix[n,n] L = rep_matrix(0,n,n);
    real s;
    real dist;
    L[1,1] = alpha;
    int idx;
    
    for(i in 2:n){
      idx = i-1;
      dist = dists[idx];
      if(mod == 0){
        L[i,1] = alpha * exp(-1.0*(dist*dist)/(theta*theta));
      }  else {
        L[i,1] = alpha * exp(-1.0*dist/theta);
      }
    }
    for(j in 2:n){
      s = 0;
      for(k in 1:(j-1)){
        s = s + L[j,k]*L[j,k];
      }
      L[j,j] = sqrt(alpha - s);
      if(j < n){
        for(i in (j+1):n){
          idx = (n-1)*(j-1)-(((j-2)*(j-1))%/%2)+(i-j-1)+1;
          dist = dists[idx];
          s = 0;
          for(k in 1:(j-1)){
            s = s + L[j,k]*L[i,k];
          }
          if(mod == 0){
            L[i,j] = 1/L[j,j] * (alpha *  exp(-(dist*dist)/(theta*theta)) - s);
          }  else {
            L[i,j] = 1/L[j,j] * (alpha * exp(-dist/theta) - s);
          }
        }
      }
    }
    return L;
  }
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=0> Q_g; //number of covariates
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  int n_region; // number of regions
  int n_Q; // number of intersections
  int<lower=1> n_cell[n_region+1]; //number of cells intersecting region  
  int<lower=1> cell_id[n_Q]; // IDs of the cells intersecting the region
  vector[n_Q] q_weights; // proportionate weights
  
  int y[Nsample*nT]; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[Nsample*nT] popdens; //population density
  matrix[Nsample*nT,Q] X;
  matrix[Nsample*nT,Q_g == 0 ? 1 : Q_g] X_g;
  real prior_lscale[2];
  real prior_var[2];
  real prior_linpred_mean[Q];
  real prior_linpred_sd[Q];
  int mod;
  int<lower = 0, upper = 1> known_cov;
  real<lower=0> sigma_data; 
  real<lower=0> phi_data; 
}
transformed data {
  vector[Nsample*nT] logpopdens = log(popdens);
  matrix[known_cov ? Nsample : 0, known_cov ? Nsample : 0] L_data;
  real dists[(Nsample*(Nsample-1))%/%2];
  for(i in 1:(Nsample-1)){
    for(j in (i+1):Nsample){
      dists[(Nsample-1)*(i-1)-(((i-2)*(i-1))%/%2)+(j-i-1)+1] = sqrt((x_grid[i,1] - x_grid[j,1]) * (x_grid[i,1] - x_grid[j,1]) +
              (x_grid[i,2] - x_grid[j,2]) * (x_grid[i,2] - x_grid[j,2]));
    }
  }
  if(known_cov){
    L_data = genChol(Nsample, sigma_data, phi_data, dists, mod);
  }
}
parameters {
  real<lower=1e-05> phi_param[known_cov ? 0 : 1]; //length scale
  real<lower=1e-05> sigma_param[known_cov ? 0 : 1];
  vector[Q] gamma;
  vector[Q_g] gamma_g;
  real<lower=-1,upper=1> ar[nT > 1 ? 0 : 1];
  vector[Nsample*nT] f_raw;
}
transformed parameters {
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
  
  if(!known_cov) {
    L = genChol(Nsample, sigma, phi, dists, mod);
  } else {
    L = L_data;
  }
  
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        f[1:Nsample] = (1/(1-ar[1]^2))*f_raw[1:Nsample];
      } else {
        f[(Nsample*(t-1)+1):(t*Nsample)] = ar[1]*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + f_raw[(Nsample*(t-1)+1):(t*Nsample)];
      }
    } else {
      f = f_raw;
    }
  }
  if(Q_g > 0){
    f += X_g * gamma_g;
  }
}
model{
  vector[n_region*nT] lambda_r = rep_vector(0,n_region*nT);
  vector[Nsample] zeros = rep_vector(0,Nsample);
  if(!known_cov){
    phi ~ normal(prior_lscale[1],prior_lscale[2]);
    sigma ~ normal(prior_var[1],prior_var[2]);
  }
  if(nT > 1)ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  if(Q_g > 0){
    gamma_g ~ normal(0,2);
  }
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        f_raw[1:Nsample] ~ multi_normal_cholesky(zeros,L);
      } else {
        f_raw[(Nsample*(t-1)+1):(t*Nsample)] ~ multi_normal_cholesky(zeros,L);
      }
    } else {
      f_raw ~ multi_normal_cholesky(zeros,L);
    }
  }
  
  for(r in 1:n_region){
    for(t in 1:nT){
      for(l in 1:(n_cell[r+1]-n_cell[r])){
        lambda_r[r+(t-1)*n_region] += popdens[r+(t-1)*n_region]*exp(X[r+(t-1)*n_region,]*gamma)*
          q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*Nsample]);
      }
    }
  }
  
  y ~ poisson_log(X*gamma+logpopdens+f);
}
generated quantities{
  vector[Nsample*nT] y_grid_predict;
  vector[n_region*nT] region_predict;
  
  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(X[i,]*gamma+logpopdens[i]+f[i]);
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

