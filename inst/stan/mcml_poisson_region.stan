functions {
  vector gen_lambda(int nT,int nR, int nS, vector xb, 
                vector w, vector gamma, int[] cell_id, 
                int[] n_cell){ 
    vector[nT*nR] lambda;
    vector[nS] u;
    for(t in 1:nT){
      u = gamma[((t-1)*nS+1):(t*nS)];
      for(r in 1:nR){
        int lsize = n_cell[r+1] - n_cell[r];
        int idx[lsize];
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
  vector[N*nT] Xb_cell;
  matrix[N*nT,N*nT] ZL;
  int y[nRegion*nT];
  int<lower=1> n_cell[nRegion+1]; //number of cells intersecting region  
  int<lower=1> cell_id[n_Q]; // IDs of the cells intersecting the region
  vector[n_Q] q_weights; 
}
parameters {
  vector[N*nT] gamma;
}
model {
  vector[N*nT] u = Xb_cell + ZL*gamma;
  vector[nRegion*nT] mu = gen_lambda(nT,nRegion,N,Xb,q_weights,u,cell_id,n_cell);
  gamma ~ std_normal();
  y ~ poisson(mu);
}

