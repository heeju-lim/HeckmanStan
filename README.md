
<!--  README.md is generated from README.Rmd. Please edit that file -->

# HeckmanStan

[![CRAN_Status_Badge]()]()

The goal of skewlmm is to fit skew robust linear mixed models, using
scale mixture of skew-normal linear mixed models with possible
within-subject dependence structure, using an EM-type algorithm. In
addition, some tools for model adequacy evaluation are available.

For more information about the model formulation and estimation, please
see Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale
mixture of skew‐normal linear mixed models with within‐subject serial
dependence. *Statistics in Medicine*. DOI:
[10.1002/sim.8870](https://doi.org/10.1002/sim.8870).

## Installation

<!-- You can install the released version of lmmsmsn from [CRAN](https://CRAN.R-project.org) with: -->

You can install skewlmm from GitHub with:

``` r
remotes::install_github("heeju-lim/HeckmanStan")
```

Or you can install the released version of skewlmm from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("HeckmanStan")
```

## Example

This is a basic example which shows you how to fit the Heckman Selection
model using bayesian approach:

``` r
library(HeckmanStan)
n<- 100
w<- cbind(1,rnorm(n),rnorm(n))
x<- cbind(w[,1:2])
type="CN"
sigma2<- 1
rho<-0.7
beta<- c(1,0.5)
gama<- c(1,0.3,-.5)
nu=c(0.1,0.1)
data<-geraHeckman(x,w,beta,gama,sigma2,rho,nu,type=type)
y<-data$y
cc<-data$cc
# Fit Heckman Normal Stan model
fit.n_stan <- HeckmanStan(y, x, w, cc, family="Normal", thin = 5, chains = 1, iter = 1000, warmup = 100)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8.3e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.83 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 101 / 1000 [ 10%]  (Sampling)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Sampling)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Sampling)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Sampling)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.052 seconds (Warm-up)
#> Chain 1:                0.424 seconds (Sampling)
#> Chain 1:                0.476 seconds (Total)
#> Chain 1:
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> 
#> -------------------------------------------------------------
#> Posterior mean(Mean), standard deviation(Sd) and HPD interval
#> -------------------------------------------------------------
#>              Mean      Sd  HPD(95%) Lower Upper Bound
#> beta[1]   1.27885 0.27396         0.79400     1.75725
#> beta[2]   0.49754 0.14095         0.29748     0.78284
#> gamma[1]  0.82469 0.20499         0.56272     1.23160
#> gamma[2]  0.23784 0.22621        -0.00464     0.51124
#> gamma[3] -0.58268 0.25725        -0.93358    -0.29479
#> rho      -0.00194 0.40346        -0.82199     0.68585
#> sigma_e   1.12353 0.11410         0.94474     1.36140
#> sigma2    1.27526 0.26408         0.89254     1.85342
#> -------------------------------------------------------------
#> Warning: 
#> 1 (1.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.
#> Model selection criteria
#> ----------------------------------------
#>          Looic     WAIC       CPO
#> Value 330.8016 330.6705 -165.3123
#> ----------------------------------------
#> 
```

``` r
qoi=c("beta","gamma","sigma_e","sigma2", "rho","EAIC","EBIC")
print(fit.n_stan[[1]],par=qoi)
#> Inference for Stan model: anon_model.
#> 1 chains, each with iter=1000; warmup=100; thin=5; 
#> post-warmup draws per chain=180, total post-warmup draws=180.
#> 
#>            mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#> beta[1]    1.28    0.02 0.24   0.85   1.13   1.26   1.43   1.80   136 1.01
#> beta[2]    0.50    0.01 0.13   0.27   0.40   0.51   0.60   0.76   270 1.00
#> gamma[1]   0.84    0.01 0.16   0.55   0.74   0.82   0.95   1.20   143 1.00
#> gamma[2]   0.26    0.01 0.13   0.03   0.18   0.26   0.35   0.53   153 1.00
#> gamma[3]  -0.61    0.01 0.16  -0.93  -0.71  -0.60  -0.49  -0.34   132 0.99
#> sigma_e    1.13    0.01 0.11   0.95   1.06   1.10   1.18   1.39   196 1.00
#> sigma2     1.28    0.02 0.26   0.91   1.12   1.22   1.39   1.92   192 1.00
#> rho        0.01    0.04 0.41  -0.81  -0.25   0.03   0.31   0.70   134 1.02
#> EAIC     337.49    0.27 3.48 332.21 334.72 336.98 339.68 344.54   166 1.02
#> EBIC     355.73    0.27 3.48 350.45 352.96 355.22 357.92 362.78   166 1.02
#> 
#> Samples were drawn using NUTS(diag_e) at Tue Oct 22 18:23:57 2024.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

``` r
print(fit.n_stan[[2]])
#> [[1]]
#>                  Mean        Sd  HPD(95%) Lower Upper Bound
#> beta[1]   1.278845783 0.2739615     0.794002748   1.7572510
#> beta[2]   0.497538298 0.1409466     0.297481796   0.7828414
#> gamma[1]  0.824691394 0.2049882     0.562717297   1.2315960
#> gamma[2]  0.237836930 0.2262101    -0.004637697   0.5112411
#> gamma[3] -0.582682024 0.2572459    -0.933581977  -0.2947870
#> rho      -0.001937362 0.4034588    -0.821986105   0.6858459
#> sigma_e   1.123525929 0.1140995     0.944742267   1.3614039
#> sigma2    1.275264107 0.2640751     0.892537950   1.8534205
#> 
#> [[2]]
#>          Looic     WAIC       CPO
#> Value 330.8016 330.6705 -165.3123
```

``` r

# Plots for stanfit objects : 
library(rstan)
#> 필요한 패키지를 로딩중입니다: StanHeaders
#> 
#> rstan version 2.32.6 (Stan version 2.32.2)
#> For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores()).
#> To avoid recompilation of unchanged Stan programs, we recommend calling
#> rstan_options(auto_write = TRUE)
#> For within-chain threading using `reduce_sum()` or `map_rect()` Stan functions,
#> change `threads_per_chain` option:
#> rstan_options(threads_per_chain = 1)
#> Do not specify '-march=native' in 'LOCAL_CPPFLAGS' or a Makevars file
```

``` r
plot(fit.n_stan[[1]], pars=c("beta[1]","beta[2]", "gamma[1]", "gamma[2]", "gamma[3]", "rho", "sigma_e"))
#> ci_level: 0.8 (80% intervals)
#> outer_level: 0.95 (95% intervals)
```

<img src="man/figures/README-example1-1.png" width="70%" style="display: block; margin: auto;" />

``` r
plot(fit.n_stan[[1]], plotfun="hist", pars=c("beta[1]","beta[2]", "gamma[1]", "gamma[2]", "gamma[3]", "rho", "sigma_e"))
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="man/figures/README-example1-2.png" width="70%" style="display: block; margin: auto;" />

``` r
plot(fit.n_stan[[1]], plotfun="trace", pars=c("beta[1]","beta[2]", "gamma[1]", "gamma[2]", "gamma[3]", "rho", "sigma_e"))
```

<img src="man/figures/README-example1-3.png" width="70%" style="display: block; margin: auto;" />
