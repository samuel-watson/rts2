functions {
  matrix getAD(real alpha, real theta, int M, int n,
                 array[] real dists, array[,] int NN, int mod){
    matrix[M+1,n] AD = rep_matrix(0,M+1,n);
    int idxlim;
    int idx1; 
    int idx2;
    real dist;
    AD[M+1,1] = alpha;
    matrix[M,M] smat = rep_matrix(0,M,M);
    vector[M] svec = rep_vector(0,M);
    int index;
    
    for(i in 2:n){
      idxlim = i<=(M) ? i-1 : M;
      for(j in 1:idxlim){
        smat[j,j] = alpha;
      }
      if(idxlim > 1){
        for(j in 1:(idxlim-1)){
          for(k in (j+1):idxlim){
            idx1 = NN[j,i] < NN[k,i] ? NN[j,i] : NN[k,i];
            idx2 = NN[j,i] < NN[k,i] ? NN[k,i] : NN[j,i];
            index = (n-1)*(idx1-1)-(((idx1-2)*(idx1-1))%/%2)+(idx2-idx1-1)+1;
            dist = dists[index];
            if(mod == 0){
              smat[j,k] = alpha * exp(-1.0*(dist*dist)/(theta*theta));
              smat[k,j] = smat[j,k];
            } else {
              smat[j,k] = alpha * exp(-1.0*dist/theta);
              smat[k,j] = smat[j,k];
            }
          }
        }
      }
      for(j in 1:idxlim){
        index = (n-1)*(NN[j,i]-1)-(((NN[j,i]-2)*(NN[j,i]-1))%/%2)+(i - NN[j,i] -1)+1;
        dist = dists[index];
        svec[j] = mod == 0 ? alpha * exp(-1.0*(dist*dist)/(theta*theta)) : alpha * exp(-1.0*dist/theta);
      }
      AD[1:idxlim,i] = mdivide_left_spd(smat[1:idxlim,1:idxlim] , svec[1:idxlim]);
      AD[M+1,i] = alpha - dot_product(AD[1:idxlim,i],svec[1:idxlim]);
    }
    return(AD);
   }
   
   real nngp_split_lpdf(array[] real u, matrix AD, array[,] int NN, int start){
    int n = cols(AD);
    int M = rows(AD) - 1;
    real logdetD;
    real qf;
    real au;
    real ll;
    int idxlim;
    if(start <= M){
      idxlim = start - 1;
    } else {
      idxlim = M;
    }
    logdetD = 0;
    for(i in 1:n){
      logdetD += log(AD[M+1,i]);
    }
    qf = u[1]*u[1]/AD[M+1,1];
    for(i in 2:n){
      au = u[i] - dot_product(AD[1:idxlim,i],to_vector(u[NN[1:idxlim,i]]));
      qf += au*au/AD[M+1,i];
      if(idxlim < M){
        idxlim+=1;
      }
    }
    
    ll = -0.5*logdetD - 0.5*qf - 0.5*n*log(2*pi());
    return ll;
   }
   
   real partial_sum_lpdf(array[] real u, int start, int end, matrix AD, array[,] int NN){
     return nngp_split_lpdf(u[start:end] | AD[,start:end], NN[,start:end], start);
   }
  
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
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=0> Q_g; //number of covariates
  int<lower=1> M; // number of nearest neighbours
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  array[M,Nsample] int NN;
  int n_region; // number of regions
  int n_Q; // number of intersections
  array[n_region+1] int<lower=1> n_cell; //number of cells intersecting region  
  array[n_Q] int<lower=1> cell_id; // IDs of the cells intersecting the region
  vector[n_Q] q_weights; // proportionate weights
  
  array[n_region*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[n_region*nT] popdens; //population density
  matrix[n_region*nT,Q] X;
  matrix[Nsample*nT,Q_g == 0 ? 1 : Q_g] X_g;
  
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
  matrix[known_cov ? M+1 : 0,known_cov ? Nsample : 0] AD_data;
  array[(Nsample*(Nsample-1))%/%2] real dists;
  
  for(i in 1:(Nsample-1)){
    for(j in (i+1):Nsample){
      dists[(Nsample-1)*(i-1)-(((i-2)*(i-1))%/%2)+(j-i-1)+1] = sqrt((x_grid[i,1] - x_grid[j,1]) * (x_grid[i,1] - x_grid[j,1]) +
              (x_grid[i,2] - x_grid[j,2]) * (x_grid[i,2] - x_grid[j,2]));
    }
  }
  
  if(known_cov){
    AD_data = getAD(sigma_data, phi_data, M, Nsample, dists, NN, mod);
  }
}

parameters {
  array[known_cov ? 0 : 1] real<lower=1e-05> phi_param; //length scale
  array[known_cov ? 0 : 1] real<lower=1e-05> sigma_param;
  vector[Q] gamma;
  vector[Q_g] gamma_g;
  array[nT > 1 ? 1 : 0] real<lower=-1,upper=1> ar;
  vector[Nsample*nT] f_raw;
}

transformed parameters {
  matrix[M +1,Nsample] AD;
  real<lower=1e-05> phi; //length scale
  real<lower=1e-05> sigma;
  vector[Nsample*nT] f;
  if(known_cov){
    sigma = sigma_data;
    phi = phi_data;
    AD = AD_data;
  } else {
    sigma = sigma_param[1];
    phi = phi_param[1];
    AD = getAD(sigma, phi, M, Nsample, dists, NN, mod);
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
  
}

model{
  vector[n_region*nT] lambda_r = rep_vector(0,n_region*nT);
  vector[Nsample*nT] region_mean = rep_vector(0,n_region*nT);
  if(!known_cov){
    phi_param ~ normal(prior_lscale[1],prior_lscale[2]);
    sigma_param ~ normal(prior_var[1],prior_var[2]);
  }
  if(nT>1)ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  
  int grainsize = 1;
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        target += reduce_sum(partial_sum_lpdf,to_array_1d(f_raw[1:Nsample]),grainsize,AD,NN);
      } else {
        target += reduce_sum(partial_sum_lpdf,to_array_1d(f_raw[(Nsample*(t-1)+1):(t*Nsample)]),grainsize,AD,NN);
      }
    } else {
      target += reduce_sum(partial_sum_lpdf,to_array_1d(f_raw),grainsize,AD,NN);
    }
  }
  
  if(Q_g > 0){
    gamma_g ~ normal(0,2);
    region_mean = X_g * gamma_g;
  }
  
  array[n_region*nT] int idx = linspaced_int_array(n_region*nT,1,n_region*nT);
  y ~ poisson_log_block(idx, nT, n_region, Nsample, popdens, X, q_weights, f+region_mean, cell_id, n_cell, gamma);
}

generated quantities{
  vector[Nsample*nT] y_grid_predict;
  vector[n_region*nT] region_predict;
  vector[Nsample*nT] region_mean_predict = rep_vector(0,n_region*nT);
  if(Q_g > 0){
    region_mean_predict = X_g * gamma_g;
  }

  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(f[i]);
  }
  
  region_predict = rep_vector(0,n_region*nT);
  for(r in 1:n_region){
    for(t in 1:nT){
      for(l in 1:(n_cell[r+1]-n_cell[r])){
        region_predict[r+(t-1)*n_region] += popdens[r+(t-1)*n_region]*exp(X[r+(t-1)*n_region,]*gamma)*
          q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*Nsample] + region_mean_predict[cell_id[n_cell[r]+l-1] + (t-1)*Nsample]);
      }
    }
  }
}
