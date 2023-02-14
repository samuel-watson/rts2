functions {
  real partial_sum1_lpdf(array[] real y, int start, int end){
    return std_normal_lpdf(y[start:end]);
  }
  real partial_sum2_lpmf(array[] int y,int start, int end, vector mu){
    return poisson_lpmf(y[start:end]|mu[start:end]);
  }
  vector gen_lambda(int nT,int nR, int nS, vector xb, 
                vector w, vector gamma, array[] int cell_id, 
                array[] int n_cell){ 
    vector[nT*nR] lambda;
    vector[nS] u;
    for(t in 1:nT){
      u = gamma[((t-1)*nS+1):(t*nS)];
      for(r in 1:nR){
        int lsize = n_cell[r+1] - n_cell[r];
        array[lsize] int idx;
        vector[lsize] f = rep_vector(0,lsize);
        for(l in 1:lsize){
          f[l] = gamma[cell_id[n_cell[r]+l-1] +(t-1)*nS];
        }
        lambda[(t-1)*nR + r] = exp(xb[r+(t-1)*nR])* (w[(n_cell[r]):(n_cell[r+1]-1)]' * exp(f));
      }
    }
    return lambda;                
  }
}
data {
  int N; // number of cells
  int nT; // number of time periods 
  int nRegion; // number of regions
  int n_Q;
  vector[nRegion*nT] Xb;
  matrix[N*nT,N*nT] ZL;
  array[nRegion*nT] int y;
  array[nRegion+1] int<lower=1> n_cell; //number of cells intersecting region  
  array[n_Q] int<lower=1> cell_id; // IDs of the cells intersecting the region
  vector[n_Q] q_weights; 
}
parameters {
  array[N*nT] real gamma;
}
transformed parameters {
  
}
model {
  int grainsize = 1;
  vector[N*nT] u = ZL*to_vector(gamma);
  vector[nRegion*nT] mu = gen_lambda(nT,nRegion,N,Xb,q_weights,u,cell_id,n_cell);
  target += reduce_sum(partial_sum1_lpdf,gamma,grainsize);
  target += reduce_sum(partial_sum2_lpmf,y,grainsize,mu);
}

