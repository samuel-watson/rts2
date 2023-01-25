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
   
   real nngp_split_lpdf(array[] real u, matrix AD, array[,] int NN, int start){
    int n = cols(AD);
    int M = rows(AD) - 1;
    real logdetD;
    real qf;
    real au;
    real ll;
    int idxlim;
    matrix[M,n] A = AD[1:M,];
    vector[n] D = AD[M+1,]';
    if(start <= M){
      idxlim = start - 1;
    } else {
      idxlim = M;
    }
    logdetD = 0;
    for(i in 1:n){
      logdetD += log(D[i]);
    }
    qf = u[1]*u[1]/D[1];
    for(i in 2:n){
      au = u[i] - dot_product(A[1:idxlim,i],to_vector(u[NN[1:idxlim,i]]));
      qf += au*au/D[i];
      if(idxlim < M){
        idxlim+=1;
      }
    }
    ll = -0.5*logdetD - 0.5*qf - 0.5*n*pi();
    return ll;
   }
   
   real partial_sum_lpdf(array[] real u, int start, int end, matrix AD, array[,] int NN){
     return nngp_split_lpdf(u | AD[,start:end], NN[,start:end], start);
   }
   
   real partial_sum2_lpmf(array[] int y,int start, int end, vector mu){
    return poisson_log_lpmf(y[start:end]|mu[start:end]);
  }
   
   // real nngp_lpdf(vector u, matrix AD, int[,] NN){
   //  int n = cols(AD);
   //  int M = rows(AD) - 1;
   //  real logdetD;
   //  real qf;
   //  real au;
   //  real ll;
   //  int idxlim;
   //  matrix[M,n] A = AD[1:M,];
   //  vector[n] D = AD[M+1,]';
   //  
   //   logdetD = 0;
   //  for(i in 1:n){
   //    logdetD += log(D[i]);
   //  }
   //  qf = u[1]*u[1]/D[1];
   //  for(i in 2:n){
   //    idxlim = i<=(M) ? i-1 : M;
   //    au = u[i] - dot_product(A[1:idxlim,i],to_vector(u[NN[1:idxlim,i]]));
   //    qf += au*au/D[i];
   //  }
   //  ll = -0.5*logdetD - 0.5*qf - 0.5*n*pi();
   //  return ll;
   // }
  
  
}
data {
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=1> M; // number of nearest neighbours
  int<lower=1> Nsample; //number of observations per time period
  int nT; //number of time periods
  array[M,Nsample] int NN;
  array[Nsample*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[Nsample*nT] popdens; //population density
  matrix[Nsample*nT,Q] X;
  array[2] real prior_lscale;
  array[2] real prior_var;
  array[Q] real prior_linpred_mean;
  array[Q] real prior_linpred_sd;
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
  
  phi ~ normal(prior_lscale[1],prior_lscale[2]);
  sigma ~ normal(prior_var[1],prior_var[2]);
  ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  
  int grainsize = 1;
  for(t in 1:nT){
    if(nT>1){
      if(t==1){
        //f_raw[1:Nsample] ~ nngp(AD, NN);
        // f_raw[1:Nsample] ~ nngp(sigma, phi, x_grid, NN, mod);
        target += reduce_sum(partial_sum_lpdf,to_array_1d(f_raw[1:Nsample]),grainsize,AD,NN);
      } else {
        target += reduce_sum(partial_sum_lpdf,to_array_1d(f_raw[(Nsample*(t-1)+1):(t*Nsample)]),grainsize,AD,NN);
        //f_raw[(Nsample*(t-1)+1):(t*Nsample)] ~ nngp(AD, NN);
        // f_raw[(Nsample*(t-1)+1):(t*Nsample)] ~ nngp(sigma, phi, x_grid, NN, mod);
      }
    } else {
      target += reduce_sum(partial_sum_lpdf,to_array_1d(f_raw),grainsize,AD,NN);
      //f_raw ~ nngp(AD, NN);
      // f_raw ~ nngp(sigma, phi, x_grid, NN, mod);
    }

  }
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,X*gamma+logpopdens+f);
  //y ~ poisson_log(X*gamma+logpopdens+f);
}

generated quantities{
  vector[Nsample*nT] y_grid_predict;
  
  for(i in 1:(Nsample*nT)){
    y_grid_predict[i] = exp(X[i,]*gamma+logpopdens[i]+f[i]);
  }
}
