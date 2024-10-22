#' Fit the Heckman Selection Stan model using the Normal, Student-t or Contaminated Normal distributions.
#'
#' `HeckmanStan()` fits the Heckman selection model using a Bayesian approach to address sample selection bias.
#'

#' @param y A response vector.
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param cc A missing indicator vector (1=observed, 0=missing) .
#' When using the CN distribution, the initial values for the proportion of bad observations and the degree of contamination.
#' @param family The family to be used (Normal, T or CN).
#' @param init Parameters specifies the initial values for model parameters.
#' @param thin An Interval at which samples are retained from the MCMC process to reduce autocorrelation.
#' @param chains The number of chains to run during the MCMC sampling. Running multiple chains is useful for checking convergence.
#' @param iter The total number of iterations for the MCMC sampling, determining how many samples will be drawn.
#' @param warmup The number of initial iterations that will be discarded as the algorithm stabilizes before collecting samples.

#' @return An object of the class HeckmanStan with all the outputs provided from the function.
#'
#' @examples
#' \donttest{
#' ################################################################################
#' # Simulation
#' ################################################################################
#' library(mvtnorm)
#' n<- 100
#' w<- cbind(1,rnorm(n),rnorm(n))
#' x<- cbind(w[,1:2])
#' type="CN"
#' sigma2<- 1
#' rho<-0.7
#' beta<- c(1,0.5)
#' gama<- c(1,0.3,-.5)
#' nu=c(0.1,0.1)
#' data<-geraHeckman(x,w,beta,gama,sigma2,rho,nu,type=type)
#' y<-data$y
#' cc<-data$cc

#' # Fit Heckman Normal Stan model
#' fit.n_stan <- HeckmanStan(y, x, w, cc, family="Normal", thin = 5, chains = 1, iter = 10000, warmup = 1000)
#'qoi=c("beta","gamma","sigma_e","sigma2", "rho","EAIC","EBIC")
#'print(fit.n_stan[[1]],par=qoi)
#'print(fit.n_stan[[2]])
#'
#' # Plots for stanfit objects : https://mc-stan.org/rstan/reference/stanfit-method-plot.html#:~:text=plotfun.%20A%20character%20string%20naming%20the%20plotting%20function%20to%20apply
#' library(rstan)
#'plot(fit.n_stan[[1]], pars=c("beta[1]","beta[2]", "gamma[1]", "gamma[2]", "gamma[3]", "rho", "sigma_e"))
#'plot(fit.n_stan[[1]], plotfun="hist", pars=c("beta[1]","beta[2]", "gamma[1]", "gamma[2]", "gamma[3]", "rho", "sigma_e"))
#'plot(fit.n_stan[[1]], plotfun="trace", pars=c("beta[1]","beta[2]", "gamma[1]", "gamma[2]", "gamma[3]", "rho", "sigma_e"))
#'plot(fit.n_stan[[1]], plotfun = "rhat")
#'
#' # Fit Heckman Student-t Stan model
#' fit.t_stan <- HeckmanStan(y, x, w, cc, family="T", thin = 5, chains = 1, iter = 10000, warmup = 1000)
#'qoi=c("beta","gamma","sigma_e","sigma2", "rho","nu")
#'print(fit.t_stan[[1]],par=qoi)
#'print(fit.t_stan[[2]])
#'
#'
#'################################################################################
#'## Application: Wage data ##Data used in the OGUNDIMU paper.
#'################################################################################
#'library(AER)
#'data("PSID1976")
#'cc <- ifelse(PSID1976$participation=="yes", 1, 0) # cc == 0  missing data
#'attach(PSID1976)
#'
#'# Transformed the factors variables into quantitative variables
#'city     <- ifelse(city=="yes", 1, 0)     # yes == 1 and no == 0
#'college  <- ifelse(college=="yes", 1, 0)  # yes == 1 and no == 0
#'hcollege <- ifelse(hcollege=="yes", 1, 0) # yes == 1 and no == 0
#'
#'x1= hwage
#'x2= youngkids
#'x3= tax
#'x4= feducation
#'x5= education
#'x6= city

#library(HeckmanEM)
#'x <- cbind(1,x5,x6)
#'w <- cbind(1,x1,x2,x3,x4,x5,x6)
#'y <- ifelse(wage>0, log(wage), 0)
#'n<-length(cc)
#' detach(PSID1976)
#'
#' library("HeckmanEM")
#' #Normal: STAN
#' ResN<- HeckmanEM(y, x, w, cc, nu = NULL, family="Normal", error = 1e-05,iter.max = 200, im=TRUE, criteria = TRUE)
#' initf1 <- function() {list(beta=ResN$beta, gamma=ResN$gamma,sigma_e=ResN$sigma,rho=ResN$rho)}
#' fit.n_stan <-  HeckmanStan( y, x, w, cc, family="Normal", init=initf1, thin = 5, chains = 1, iter = 10000, warmup = 1000)
#' print(fit.cn_stan[[1]],par=c("beta","gamma","sigma_e","sigma2", "rho","nu1","nu2","likel","EAIC","EBIC"))
#' print(fit.cn_stan[[2]])
#'
#' initf2 <- function() {list(beta=ResN$beta, gamma=ResN$gamma,sigma_e=ResN$sigma,rho=ResN$rho, nu1=0.1,nu2=0.1)}
#' fit.cn_stan <-  HeckmanStan( y, x, w, cc, family="CN", init=initf2, thin = 5, chains = 1, iter = 10000, warmup = 1000)
#'
#' initf3 <- function() {list(beta=ResN$beta, gamma=ResN$gamma,sigma_e=ResN$sigma,rho=ResN$rho, nu=3)}
#' fit.t_stan <-  HeckmanStan( y, x, w, cc, family="T", init=initf3, thin = 5, chains = 1, iter = 10000, warmup = 1000)
#'
#' }
#'
#' @import loo
#' @importFrom rstan stan
#' @export
HeckmanStan <- function(y, x, w, cc, family="CN", init="random", thin = 5, chains = 1, iter = 10, warmup = 5){

  n<-length(cc)
  data = list(N = n, N_y = sum(cc==1), p = ncol(x), q = ncol(w), X = x[cc > 0, ], Z = w, D = cc, y = y[cc > 0])

  if (family != "Normal" && family !="normal" && family !="T" && family !="t" && family !="CN" && family !="cn" ) stop("Family not recognized! Obly families allowed are: \"Normal\", \"T\" and \"CN\".")
  if(!is.vector(y)) stop("y must be a vector!")
  if(!is.vector(cc)) stop("y must be a vector!")

  if(is.vector(x)) x <- as.matrix(x)
  if(is.vector(w)) w <- as.matrix(w)
  if(!is.matrix(x)) stop("y must be a matrix!")
  if(!is.matrix(w)) stop("y must be a matrix!")



  if(family == "Normal" || family == "normal"){

    stan_file <- system.file("stan", "HeckmanNormal.stan",package = "HeckmanStan")
    out <- rstan::stan(file=stan_file,
                       data =data , init=init,
                       thin = thin, chains = chains, iter = iter, warmup = warmup)

      }

  if((family == "T" || family == "t")){

    stan_file <- system.file("stan", "HeckmanT.stan",package = "HeckmanStan")
    out <- rstan::stan(file=stan_file,
                data =data , init=init,
                thin = thin, chains = chains, iter = iter, warmup = warmup)

  }


  if((family == "CN" || family == "cn")){

    stan_file <- system.file("stan", "HeckmanCNormal.stan",package = "HeckmanStan")
    out <- rstan::stan(file=stan_file,
                data =data , init=init,
                thin = thin, chains = chains, iter = iter, warmup = warmup)

  }

  cat('\n')
  cat('-------------------------------------------------------------\n')
  cat('Posterior mean(Mean), standard deviation(Sd) and HPD interval\n')
  cat('-------------------------------------------------------------\n')

  paramT <- matrix(nrow = 0, ncol = 4)

  nsize<-out@par_dims$beta+out@par_dims$gamma+ifelse(is.null(out@par_dims$rho),0,1)+ ifelse(is.null(out@par_dims$sigma_e),0,1)+ ifelse(is.null(out@par_dims$sigma2),0,1)+ifelse(is.null(out@par_dims$nu),0,1)+ifelse(is.null(out@par_dims$nu1),0,1)+ifelse(is.null(out@par_dims$nu2),0,1)
  for(i in 1:nsize){
  pname<-names(out@sim$samples[[1]])
  param <- mean(out@sim$samples[[1]][[pname[i]]])
  se <- sd(out@sim$samples[[1]][[pname[i]]])
  HPDTot<-hpd(out@sim$samples[[1]][[pname[i]]], alpha=0.05)

  paramT <- rbind(paramT, c(param, se, HPDTot))
  }

  dimnames(paramT) <- list(c(pname[1:nsize]),c("Mean", "Sd", " HPD(95%) Lower","Upper Bound"))
  print(round(paramT,digits=5))

  cat('-------------------------------------------------------------\n')
# fit_n <- extract(out, par=c("EAIC","EBIC"))
  log_lik_n <- loo::extract_log_lik(out, merge_chains = FALSE)
  loo_n <- loo::loo(log_lik_n)
  WAIC_n = loo::waic(log_lik_n)
  critFin<-cbind(loo_n$estimates[[3]],  WAIC_n$estimates[[3]], CPO(out))
# critFin<-cbind(mean(fit_n$EAIC), mean(fit_n$EBIC))
# c("EAIC","EBIC")
  dimnames(critFin)<-list(c("Value"),c("Looic", "WAIC","CPO"))
  cat('Model selection criteria\n')
  cat('----------------------------------------\n')
  print(critFin)

  cat('----------------------------------------\n')
  cat('\r \n')
  output<-list(paramT, critFin)

  return(list(out, output ))
}
