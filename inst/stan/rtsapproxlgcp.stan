functions {
  vector lambda_nD(array[] real L, array[] int m, int D) {
    vector[D] lam;
    for(i in 1:D){
      lam[i] = ((m[i]*pi())/(2*L[i]))^2; }

    return lam;
  }
  real spd_nD(real sigma, real phi, vector w, int D, int mod) {
    real S;
    vector[2] wsq;
    wsq = w .* w;
    
    if(mod == 0){
      // squared exponential
      S = sigma * sqrt(2*pi())^D * phi * phi * exp(-0.5*(phi*phi*(wsq[1] + wsq[2])));
    } else {
      // exponential
      S = sigma * 4 * pi() * phi * phi * (1 + phi*phi*(wsq[1] + wsq[2]))^(-1.5);
    }

    return S;
  }
  vector phi_nD(array[] real L, array[] int m, matrix x) {
    int c = cols(x);
    int r = rows(x);

    matrix[r,c] fi;
    vector[r] fi1;
    for (i in 1:c){
      fi[,i] = 1/sqrt(L[i])*sin(m[i]*pi()*(x[,i]+L[i])/(2*L[i]));
    }
    fi1 = fi[,1];
    for (i in 2:c){
      fi1 = fi1 .* fi[,i];
    }
    return fi1;
  }
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  array[D] real L; // boundary condition
  int<lower=1> M; // number of basis functions (per dimension)
  int<lower=1> M_nD; //total basis functions m1*m2*...*mD
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  array[Nsample*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  array[M_nD,D] int indices; //indices
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
  matrix[Nsample,M_nD] PHI;
  vector[Nsample*nT] logpopdens = log(popdens);
  array[known_cov ? M_nD : 0] real diagSPD_data;

  for (m in 1:M_nD){
    PHI[,m] = phi_nD(L, indices[m,], x_grid);
  }
  
  if(known_cov){
    for(m in 1:M_nD){
      diagSPD_data[m] =  sqrt(spd_nD(sigma_data, phi_data, sqrt(lambda_nD(L, indices[m,], D)), D, mod));
    }
  }
}

parameters {
  matrix[M_nD,nT] beta;
  array[known_cov ? 0 : 1] real<lower=1e-05> phi_param; //length scale
  array[known_cov ? 0 : 1] real<lower=1e-05> sigma_param;
  vector[Q] gamma;
  array[known_cov ? 0 : 1] real<lower=-1,upper=1> ar;
}

transformed parameters{
  vector[Nsample*nT] f;
  vector[M_nD] diagSPD;
  vector[M_nD] SPD_beta;
  real<lower=1e-05> sigma;
  real<lower=1e-05> phi;
  if(known_cov){
    sigma = sigma_data;
    phi = phi_data;
    diagSPD = to_vector(diagSPD_data);
  } else {
    sigma = sigma_param[1];
    phi = phi_param[1];
    for(m in 1:M_nD){
      diagSPD[m] =  sqrt(spd_nD(sigma, phi, sqrt(lambda_nD(L, indices[m,], D)), D, mod));
    }
  }

  for(t in 1:nT){
    SPD_beta = diagSPD .* beta[,t];
    if(nT>1){
      if(t==1){
        f[1:Nsample] = (1/(1-ar[1]^2))*PHI * SPD_beta;
      } else {
        f[(Nsample*(t-1)+1):(t*Nsample)] = ar[1]*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + PHI * SPD_beta;
      }
    } else {
      f[1:Nsample] = PHI * SPD_beta;
    }

  }

}
model{
  to_vector(beta) ~ normal(0,1);
  if(!known_cov){
    phi_param ~ normal(prior_lscale[1],prior_lscale[2]);
    sigma_param ~ normal(prior_var[1],prior_var[2]);
  }
  if(nT > 1)ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  y ~ poisson_log(X*gamma+logpopdens+f);

}

generated quantities{
  vector[Nsample*nT] y_grid_predict;

  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(X[i,]*gamma+logpopdens[i]+f[i]);
  }
}
