data {
  // dimensions
  int<lower=1> N;
  int<lower=1, upper=N> N_y;
  int<lower=1> p;
  int<lower=1> q;
  // covariates
  matrix[N_y, p] X;
  matrix[N, q] Z;
  // responses
  int<lower=0, upper=1> D[N];
  vector[N_y] y;
}
parameters {
  vector[p] beta;
  vector[q] gamma;
  real<lower=-1,upper=1> rho;
  real<lower=0> sigma_e;
}

transformed parameters{
  real<lower=0> sigma2=sigma_e^2;
}
model {
  // naive (truncated) priors
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  rho ~ uniform(-1, 1);
  sigma_e ~ cauchy(0, 4);
  {
    // log-likelihood
    vector[N_y] Xb = X * beta;
    vector[N] Zg = Z * gamma;
    real sigmat=(1 - rho^2);
    int ny = 1;
    for(n in 1:N) {
      if(D[n] > 0) {
        real mut= -(Zg[n] + rho / sigma_e * (y[ny] - Xb[ny]));
        target += normal_lpdf(y[ny] | Xb[ny], sigma_e) + normal_lcdf(0 | mut, sqrt(sigmat));
        ny += 1;
      }
      else {
        target += normal_lcdf(0 | Zg[n],1);
      }
    }
  }
}

generated quantities {
	vector[N] log_lik;
	real likel;
	vector[N_y] Xb = X * beta;
    vector[N] Zg = Z * gamma;
	    int ny = 1;
	    real sigmat=(1 - rho^2);
    for(n in 1:N) {
      if(D[n] > 0) {
        real mut= -(Zg[n] + rho / sigma_e * (y[ny] - Xb[ny]));
        log_lik[n] = normal_lpdf(y[ny] | Xb[ny], sigma_e) + normal_lcdf(0 | mut, sqrt(sigmat));
        ny += 1;
      }
      else {
        log_lik[n]= normal_lcdf(0 | Zg[n],1);
      }
    }
    likel=sum(log_lik);
    real EAIC=-2*likel+2*(p+q+2);
    real EBIC=-2*likel+log(N)*(p+q+2);
}