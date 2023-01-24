functions {
  matrix getAD(real alpha, real theta, 
                 matrix x, array[,] int NN, int mod){
                   int n = rows(x);
    int M = size(NN);
    matrix[M,n] A = rep_matrix(0,M,n);
    row_vector[n] D = rep_row_vector(0,n);
    matrix[M+1,n] AD;
    int idxlim;
    matrix[M,M] smat;
    vector[M] svec;
    real dist;
    
    D[1] = alpha;
    
    for(i in 2:n){
      idxlim = i<=(M) ? i-1 : M;
      smat[1:idxlim,1:idxlim] = rep_matrix(0,idxlim,idxlim);
      svec[1:idxlim] = rep_vector(0,idxlim);
      for(j in 1:idxlim){
        smat[j,j] = alpha;
      }
      if(idxlim > 1){
        for(j in 1:(idxlim-1)){
          for(k in (j+1):idxlim){
            dist = sqrt((x[NN[j,i],1] - x[NN[k,i],1]) * (x[NN[j,i],1] - x[NN[k,i],1]) +
              (x[NN[j,i],2] - x[NN[k,i],2]) * (x[NN[j,i],2] - x[NN[k,i],2]));
            if(mod == 0){
              smat[j,k] = alpha * alpha * exp(-(dist*dist)/(theta*theta));
              smat[k,j] = alpha * alpha * exp(-(dist*dist)/(theta*theta));
            } else {
              smat[j,k] = alpha * alpha * exp(-dist/theta);
              smat[k,j] = alpha * alpha * exp(-dist/theta);
            }
            
          }        
        }
      }
      for(j in 1:idxlim){
        dist = sqrt((x[NN[j,i],1] - x[i,1])^2 + (x[NN[j,i],2] - x[i,2])^2);
        svec[j] = alpha * alpha * exp(-dist/theta);
      }
      A[1:idxlim,i] = mdivide_left_spd(smat[1:idxlim,1:idxlim] , svec[1:idxlim]);
      D[i] = alpha - dot_product(A[1:idxlim,i],svec[1:idxlim]);
    }
    
    AD = append_row(A,D);
    
    return(AD);
    
   }
   
   real nngp_lpdf(vector u, matrix AD, array[,] int NN){
    int n = cols(AD);
    int M = rows(AD) - 1;
    real logdetD;
    real qf;
    real au;
    real ll;
    int idxlim;
    matrix[M,n] A = AD[1:M,];
    vector[n] D = AD[M+1,]';
    
     logdetD = 0;
    for(i in 1:n){
      logdetD += log(D[i]);
    }
    qf = u[1]*u[1]/D[1];
    for(i in 2:n){
      idxlim = i<=(M) ? i-1 : M;
      au = u[i] - dot_product(A[1:idxlim,i],to_vector(u[NN[1:idxlim,i]]));
      qf += au*au/D[i];
    }
    ll = -0.5*logdetD - 0.5*qf - 0.5*n*pi();
    return ll;
   }
  
  real calc_lambda(real p, real xg, vector w, vector f){
    return p*exp(xg)*(w' * exp(f));
  }
  real poisson_log_block_lpmf(array[] int y, array[] int i, int nT, 
                              int nR, vector pop, matrix X, 
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
      for(l in 1:lsize)idx[l] = cell_id[n_cell[r]+l-1]+(t-1)*nR;
      lambda[j] = calc_lambda(pop[r+(t-1)*nR],
                         X[r+(t-1)*nR,]*gamma,
                         w[(n_cell[r]):(n_cell[r+1]-1)],
                         f[idx]);
    }
    
    return poisson_lpmf(y | lambda);                  
                              
  }
  // real partial_sum2_lpmf(array[] int y,int start, int end, int nT,
  //                       int nR, array[] real pop, matrix X, 
  //                       vector w, vector f, array[] int cell_id, array[] int n_cell){
  //   return poisson_log_block_lpmf(y[start:end],start:end, nT, nR, pop,
  //                                 X, w, f, cell_id, n_cell);
  // }
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
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
  
  array[2] real prior_lscale;
  array[2] real prior_var;
  array[Q] real prior_linpred_mean;
  array[Q] real prior_linpred_sd;
  int mod;
}

parameters {
  real<lower=1e-05> phi; //length scale
  real<lower=1e-05> sigma;
  vector[Q] gamma;
  real<lower=-1,upper=1> ar;
  vector[Nsample*nT] f_raw;
}

transformed parameters {
  matrix[M +1,Nsample] AD;
  vector[Nsample*nT] f;
  AD = getAD(sigma, phi, x_grid, NN, mod);
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
  vector[n_region*nT] lambda_r = rep_vector(0,n_region*nT);
  phi ~ normal(prior_lscale[1],prior_lscale[2]);
  sigma ~ normal(prior_var[1],prior_var[2]);
  ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  
  
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        f_raw[1:Nsample] ~ nngp(AD, NN);
      } else {
        f_raw[(Nsample*(t-1)+1):(t*Nsample)] ~ nngp(AD, NN);
      }
    } else {
      f_raw ~ nngp(AD, NN);
    }
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
  // target += reduce_sum(partial_sum2_lpmf,y,grainsize,nT,n_region,popdens,
  //                       X, q_weights, f, cell_id, n_cell);
  array[n_region*nT] int idx = linspaced_int_array(n_region*nT,1,n_region*nT);
  y ~ poisson_log_block(idx, nT, n_region, popdens, X, q_weights, f, cell_id, n_cell, gamma);
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
          q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*n_region]);
      }
    }
  }
}
