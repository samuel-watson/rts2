functions {
  matrix genChol(int n, real alpha, real theta, array[] real dists, int mod){
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
          //dist = sqrt((x[j,1] - x[i,1]) * (x[j,1] - x[i,1]) + (x[j,2] - x[i,2]) * (x[j,2] - x[i,2]));
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
        index = (n-1)*(NN[j,i]-1)-(((NN[j,i]-2)*(NN[j,i]-1))%/%2)+(i - NN[j,i]-1)+1;
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
}
data {
  // type of model
  int mod;
  int<lower = 0, upper = 2> approx; // 0 = full, 1 = HSGP, 2 = NNGP
  int is_region; // 0 = no, 1 = yes
  
  // common data
  int<lower=1> D; //number of dimensions
  int<lower=1> Q; //number of covariates
  int<lower=1> Nsample; //number of observations per time period
  int<lower=1> n_region; // number of regions
  int nT; //number of time periods
  array[is_region ? n_region*nT : Nsample*nT] int y; //outcome
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[is_region ? n_region*nT : Nsample*nT] popdens; //population density
  matrix[is_region ? n_region*nT : Nsample*nT,Q] X;
  
  // priors
  array[2] real prior_lscale;
  array[2] real prior_var;
  array[Q] real prior_linpred_mean;
  array[Q] real prior_linpred_sd;
  
  // known covariates?
  int<lower = 0, upper = 1> known_cov;
  real<lower=0> sigma_data; 
  real<lower=0> phi_data; 
  
  // regional data
  int n_Q; // number of intersections
  array[is_region ? n_region+1 : 1] int<lower=1> n_cell; //number of cells intersecting region  
  array[is_region ? n_Q : 1] int<lower=1> cell_id; // IDs of the cells intersecting the region
  vector[is_region ? n_Q : 1] q_weights; // proportionate weights
  int<lower=0> Q_g; //number of covariates
  matrix[is_region ? Nsample*nT : 1,Q_g == 0 ? 1 : Q_g] X_g;
  
  // HSGP
  int<lower=1> M; // number of basis functions (per dimension) or nearest neighbours
  int<lower=1> M_nD; //total basis functions m1*m2*...*mD
  array[approx == 1 ? M_nD : 1,approx == 1 ? D : 1] int indices; //indices
  array[D] real L; // boundary condition
  
  // NNGP
  array[approx == 2 ? M : 1, approx == 2 ? Nsample : 1] int NN;
}
transformed data {
  vector[is_region ? n_region*nT : Nsample*nT] logpopdens = log(popdens);
  // full GP
  matrix[known_cov && approx == 0 ? Nsample : 0, known_cov  && approx == 0 ? Nsample : 0] L_data;
  array[approx != 1 ? (Nsample*(Nsample-1))%/%2 : 0] real dists;
  // HSGP
  matrix[approx == 1 ? Nsample : 0, approx == 1 ? M_nD : 0] PHI;
  array[known_cov && approx == 1 ? M_nD : 0] real diagSPD_data;
  // NNGP 
  matrix[known_cov && approx == 2 ? M+1 : 0,known_cov && approx == 2 ? Nsample : 0] AD_data;
  
  if(approx == 0 || approx == 2){
    for(i in 1:(Nsample-1)){
      for(j in (i+1):Nsample){
        dists[(Nsample-1)*(i-1)-(((i-2)*(i-1))%/%2)+(j-i-1)+1] = sqrt((x_grid[i,1] - x_grid[j,1]) * (x_grid[i,1] - x_grid[j,1]) +
                (x_grid[i,2] - x_grid[j,2]) * (x_grid[i,2] - x_grid[j,2]));
      }
    }
  }
  
  if(approx == 0){
    if(known_cov){
      L_data = genChol(Nsample, sigma_data, phi_data, dists, mod);
    }
  } else if(approx == 1){
    for (m in 1:M_nD){
      PHI[,m] = phi_nD(L, indices[m,], x_grid);
    }
    if(known_cov){
      for(m in 1:M_nD){
        diagSPD_data[m] =  sqrt(spd_nD(sigma_data, phi_data, sqrt(lambda_nD(L, indices[m,], D)), D, mod));
      }
    }
  } else if(approx == 2){
    if(known_cov){
      AD_data = getAD(sigma_data, phi_data, M, Nsample, dists, NN, mod);
    }
  }
}

parameters {
  array[known_cov ? 0 : 1] real<lower=1e-05> phi_param; //length scale
  array[known_cov ? 0 : 1] real<lower=1e-05> sigma_param;
  vector[Q] gamma;
  vector[Q_g] gamma_g;
  array[nT > 1 ? 1 : 0] real<lower=-1,upper=1> ar;
  array[(approx == 0 || approx == 2) ? Nsample*nT : 2] real f_raw;
  matrix[approx == 1 ? M_nD : 0, approx == 1 ? nT : 0] beta;
}

transformed parameters{
  vector[Nsample*nT] f;
  real<lower=1e-05> sigma;
  real<lower=1e-05> phi;
  // full GP
  matrix[approx == 0 ? Nsample : 0,approx == 0 ? Nsample : 0] Lmat;
  // HSGP
  vector[approx == 1 ? M_nD : 0] diagSPD;
  vector[approx == 1 ? M_nD : 0] SPD_beta;
  // NNGP 
  matrix[approx == 2 ? M +1 : 0,approx == 2 ? Nsample : 0] AD;
  
  if(known_cov){
    sigma = sigma_data;
    phi = phi_data;
    if(approx == 1)diagSPD = to_vector(diagSPD_data);
    if(approx == 2)AD = AD_data;
  } else {
    sigma = sigma_param[1];
    phi = phi_param[1];
    if(approx == 1){
      for(m in 1:M_nD){
        diagSPD[m] =  sqrt(spd_nD(sigma, phi, sqrt(lambda_nD(L, indices[m,], D)), D, mod));
      }
    }
    if(approx == 2){
      AD = getAD(sigma, phi, M, Nsample, dists, NN, mod);
    }
  }
  
  if(approx == 0){
    if(!known_cov) {
      Lmat = genChol(Nsample, sigma, phi, dists, mod);
    } else {
      Lmat = L_data;
    }
    for(t in 1:nT){
      if(nT>1){
        if(t==1){
          f[1:Nsample] = (1/(1-ar[1]^2))*Lmat*to_vector(f_raw[1:Nsample]);
        } else {
          f[(Nsample*(t-1)+1):(t*Nsample)] = ar[1]*Lmat*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + Lmat*to_vector(f_raw[(Nsample*(t-1)+1):(t*Nsample)]);
        }
      } else {
        f = Lmat*to_vector(f_raw);
      }
    }
  } else if(approx == 1){
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
  } else if(approx == 2){
    for(t in 1:nT){
      if(nT>1){
        if(t==1){
          f[1:Nsample] = (1/(1-ar[1]^2))*to_vector(f_raw[1:Nsample]);
        } else {
          f[(Nsample*(t-1)+1):(t*Nsample)] = ar[1]*f[(Nsample*(t-2)+1):((t-1)*Nsample)] + to_vector(f_raw[(Nsample*(t-1)+1):(t*Nsample)]);
        }
      } else {
        f = to_vector(f_raw);
      }
    }
  }
}
model{
  vector[is_region ? n_region*nT : 0] lambda_r;
  vector[is_region ? Nsample*nT : 0] region_mean;
  real accum = 0;
  
  if(is_region){
    lambda_r = rep_vector(0,n_region*nT);
    region_mean = rep_vector(0,Nsample*nT);
  }
  
  if(approx == 1)to_vector(beta) ~ normal(0,1);
  if(!known_cov){
    phi_param ~ normal(prior_lscale[1],prior_lscale[2]);
    sigma_param ~ normal(prior_var[1],prior_var[2]);
  }
  if(nT > 1)ar ~ normal(0,1);
  for(q in 1:Q){
    gamma[q] ~ normal(prior_linpred_mean[q],prior_linpred_sd[q]);
  }
  
  if(approx == 2){
      for(t in 1:nT){
        if(nT>1){
          if(t==1){
            target += nngp_split_lpdf(to_array_1d(f_raw[1:Nsample])|AD,NN,1);
          } else {
            target += nngp_split_lpdf(to_array_1d(f_raw[(Nsample*(t-1)+1):(t*Nsample)])|AD,NN,1);
          }
        } else {
          target += nngp_split_lpdf(to_array_1d(f_raw)|AD,NN,1);
        }
      }
    }
  
  if(is_region){
    if(Q_g > 0){
      gamma_g ~ normal(0,2);
      region_mean = X_g * gamma_g;
    }
    
    if(approx == 0){
      f_raw ~ std_normal();
    } 
    
    for(r in 1:n_region)
    {
      for(t in 1:nT)
      {
        accum = 0;
        lambda_r[r+(t-1)*n_region] = popdens[r+(t-1)*n_region]*exp(X[r+(t-1)*n_region,]*gamma);
        for(l in 1:(n_cell[r+1]-n_cell[r]))
        {
          accum += q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*Nsample] + region_mean[cell_id[n_cell[r]+l-1] + (t-1)*Nsample]);
        }
        lambda_r[r+(t-1)*n_region] *= accum;
      }
    }
    y ~ poisson(lambda_r);
    
  } else {
    y ~ poisson_log(X*gamma+logpopdens+f);
  } 
  
}

generated quantities{
  vector[Nsample*nT] y_grid_predict;
  vector[is_region ? n_region*nT : 0] region_predict;
  vector[is_region ? Nsample*nT : 0] region_mean_predict;
  
  if(!is_region){
    for(i in 1:(Nsample*nT)){
      y_grid_predict[i] = exp(X[i,]*gamma+logpopdens[i]+f[i]);
    }
  } else {
    region_predict = rep_vector(0,n_region*nT);
    region_mean_predict = rep_vector(0,Nsample*nT);
    if(Q_g > 0){
      region_mean_predict = X_g * gamma_g;
    }
    for(i in 1:(Nsample*nT))
    {
      y_grid_predict[i] = exp(f[i] + region_mean_predict[i]);
    }
    for(r in 1:n_region){
      for(t in 1:nT){
        for(l in 1:(n_cell[r+1]-n_cell[r])){
          region_predict[r+(t-1)*n_region] += popdens[r+(t-1)*n_region]*exp(X[r+(t-1)*n_region,]*gamma)*
            q_weights[n_cell[r]+l-1]*exp(f[cell_id[n_cell[r]+l-1] + (t-1)*Nsample]+region_mean_predict[cell_id[n_cell[r]+l-1] + (t-1)*Nsample]);
        }
      }
    }
  }
}
