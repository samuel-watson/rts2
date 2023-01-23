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
  real calc_lambda(real p, real xg, vector w, vector f){
    return p*exp(xg)*(w' * exp(f));
  }
  real poisson_log_block_lpmf(array[] int y, array[] int i, int nT, 
                              int nR, array[] real pop, matrix X, 
                vector w, vector f, array[] int cell_id, array[] int n_cell){ 
    // calculate the IDs of the observation
    int n = size(y);
    vector[n] lambda;
    for(j in 1:n){
      int r = i%nR;
      if(r == 0)r = nR;
      int t = ceil(i/nR);
      lambda[j] = calc_lambda(pop[r+(t-1)*nT],
                         X[r+(t-1)*nT,]*gamma,
                         w[(n_cell[r]):(n_cell[r+1]-1)],
                         f[cell_id[(n_cell[r]):(n_cell[r+1]-1)] + (t-1)*nT]);
    }
    
    return poisson_lpmf(y | lambda);                  
                              
  }
  real partial_sum2_lpmf(array[] int y,int start, int end, int nT,
                        int nR, array[] real pop, matrix X, 
                        vector w, vector f, array[] int cell_id, array[] int n_cell){
    return poisson_log_block_lpmf(y[start:end],start:end, nT, nR, pop,
                                  X, w, f, cell_id, n_cell);
  }
}
data {
  // define the model and problem
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  array[D] real L; // boundary condition
  int<lower=1> M; // number of basis functions (per dimension)
  int<lower=1> M_nD; //total basis functions m1*m2*...*mD
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  int n_region; // number of regions
  int n_Q; // number of intersections
  array[n_region] int<lower=1> n_cell; //number of cells intersecting region  
  array[n_Q] int<lower=1> cell_id; // IDs of the cells intersecting the region
  vector[n_Q] q_weights; // proportionate weights
  
  // outcomes data
  array[n_region*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  array[M_nD,D] int indices; //indices
  vector[n_region*nT] popdens; //population density
  matrix[n_region*nT,Q] X;
  
  // distribution parameters
  real prior_lscale[2];
  real prior_var[2];
  real prior_linpred_mean[Q];
  real prior_linpred_sd[Q];
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
  vector[n_region*nT] lambda_r = rep_vector(0,n_region*nT);
  to_vector(beta) ~ normal(0,1);
  phi ~ normal(prior_lscale[1],prior_lscale[2]);
  sigma ~ normal(prior_var[1],prior_var[2]);
  ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  int grainsize = 1;
  // for(r in 1:n_region){
  //   for(t in 1:nT){
  //     for(l in 1:(n_cell[r+1]-n_cell[r])){
  //       lambda_r[r+(t-1)*nT] += popdens[r+(t-1)*nT]*exp(X[r+(t-1)*nT,]*gamma)*
  //         q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*nT]);
  //     }
  //   }
  // }
  // 
  // y ~ poisson(lambda_r); // we can parallelise this 
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,nT,n_region,popdens,
                        X, q_weights, f, cell_id, n_cell);
}

generated quantities{
  vector[Nsample*nT] y_grid_predict;
  vector[n_region*nT] region_predict;

  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(f[i]);
  }
  
  for(r in 1:n_region){
    for(t in 1:nT){
      for(l in 1:(n_cell[r+1]-n_cell[r])){
        region_predict[r+(t-1)*nT] += popdens[r+(t-1)*nT]*exp(X[r+(t-1)*nT,]*gamma)*
          q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*nT]);
      }
    }
  }
  
}
