functions {
  vector lambda_nD(array[] real L, array[] int m, int D) {
    vector[D] lam;
    for(i in 1:D){
      lam[i] = ((m[i]*pi())/(2*L[i]))^2; }

    return lam;
  }
  real spd_nD(real sigma, row_vector phi, vector w, int D, int mod) {
    real S;
    real S1;
    vector[2] phisq; 
    vector[2] wsq;
    phisq = (phi .* phi)';
    wsq = w .* w;
    
    if(mod == 0){
      // squared exponential
      S = sigma^2 * sqrt(2*pi())^D * prod(phi) * exp(-0.5*((phi .* phi) * (w .* w)));
    } else {
      // exponential
      S1 = sigma^2 * 4 * pi() * tgamma(1.5)/tgamma(0.5);
      S = S1 * prod(phi) * (1 + phisq' * wsq)^(-2);
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
}
transformed data {
  matrix[Nsample,M_nD] PHI;
  vector[Nsample*nT] logpopdens = log(popdens);

  for (m in 1:M_nD){
    PHI[,m] = phi_nD(L, indices[m,], x_grid);
  }
}

parameters {
  matrix[M_nD,nT] beta;
  row_vector<lower=1e-05>[D] phi; //length scale
  real<lower=1e-05> sigma;
  vector[Q] gamma;
  real<lower=-1,upper=1> ar;
}

transformed parameters{
  vector[Nsample*nT] f;
  //vector[Nsample] f_tilde;
  vector[M_nD] diagSPD;
  vector[M_nD] SPD_beta;

  for(m in 1:M_nD){
    diagSPD[m] =  sqrt(spd_nD(sigma, phi, sqrt(lambda_nD(L, indices[m,], D)), D, mod));
  }

  for(t in 1:nT){
    SPD_beta = diagSPD .* beta[,t];
    if(nT>1){
      if(t==1){
        f[1:Nsample] = (1/(1-ar^2))*PHI * SPD_beta;
      } else {
        f[(Nsample*(t-1)+1):(t*Nsample)] = ar*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + PHI * SPD_beta;
      }
    } else {
      f[1:Nsample] = PHI * SPD_beta;
    }

  }

}
model{
  to_vector(beta) ~ normal(0,1);
  phi ~ normal(prior_lscale[1],prior_lscale[2]);
  sigma ~ normal(prior_var[1],prior_var[2]);
  ar ~ normal(0,1);
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
