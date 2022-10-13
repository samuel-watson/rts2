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
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  int y[Nsample*nT]; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[Nsample*nT] popdens; //population density
  matrix[Nsample*nT,Q] X;
  real prior_lscale[2];
  real prior_var[2];
  real prior_linpred_mean[Q];
  real prior_linpred_sd[Q];
  int mod;
}

transformed data {
  
  vector[Nsample*nT] logpopdens = log(popdens);
}

parameters {
  real<lower=1e-05> phi; //length scale
  real<lower=1e-05> sigma;
  vector[Q] gamma;
  real<lower=-1,upper=1> ar;
  vector[Nsample*nT] f_raw;
}

transformed parameters {
  matrix[Nsample,Nsample] L;
  vector[Nsample*nT] f;
  L = genChol(Nsample, sigma, phi, x_grid, mod);
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        f[1:Nsample] = (1/(1-ar^2))*f_raw[1:Nsample];
      } else {
        f[(Nsample*(t-1)+1):(t*Nsample)] = ar*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + f_raw[(Nsample*(t-1)+1):(t*Nsample)];
      }
    } else {
      f = f_raw;
    }

  }
}

model{
  vector[Nsample] zeros = rep_vector(0,Nsample);
  phi ~ normal(prior_lscale[1],prior_lscale[2]);
  sigma ~ normal(prior_var[1],prior_var[2]);
  ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
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
  
  y ~ poisson_log(X*gamma+logpopdens+f);
}

generated quantities{
  vector[Nsample*nT] y_grid_predict;
  
  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(X[i,]*gamma+logpopdens[i]+f[i]);
  }
}

