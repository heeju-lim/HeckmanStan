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
  beta ~ normal(0, 10);
  gamma ~ normal(0, 10);
  rho ~ uniform(-1, 1);
  sigma_e ~ cauchy(0, 4);
  {
    // log-likelihood
    vector[N_y] Xb = X * beta;
    vector[N] Zg = Z * gamma;
    int ny = 1;
    for(n in 1:N) {
      if(D[n] > 0) {
        target += normal_lpdf(y[ny] | Xb[ny], sigma_e) + log(Phi((Zg[n] + rho / sigma_e * (y[ny] - Xb[ny])) / sqrt(1 - rho^2)));
        ny += 1;
      }
      else {
        target += log(Phi(-Zg[n]));
      }
    }
  }
}

generated quantities {
	vector[N] log_lik;
	vector[N] log_invlik;
	vector[N_y] Xb = X * beta;
    vector[N] Zg = Z * gamma;
	    int ny = 1;
    for(n in 1:N) {
      if(D[n] > 0) {
        log_lik[n] = normal_lpdf(y[ny] | Xb[ny], sigma_e) + log(Phi((Zg[n] + rho / sigma_e * (y[ny] - Xb[ny])) / sqrt(1 - rho^2)));
        ny += 1;
      }
      else {
        log_lik[n]= log(Phi(-Zg[n]));
      }
      log_invlik[n]= exp(-log_lik[n]);
    }
    real likel=sum(log_lik);
    real EAIC= -2*likel+2*(p+q+2);
    real EBIC= -2*likel+log(N)*(p+q+2);
}
