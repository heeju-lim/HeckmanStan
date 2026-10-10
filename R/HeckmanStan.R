#' Fit the Heckman Selection Stan model using the Normal, Student-t or Contaminated Normal distributions.
#'
#' `HeckmanStan()` fits the Heckman selection model using a Bayesian approach to address sample selection bias.
#'

#' @param y A response vector.
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param cc A missing indicator vector (1=observed, 0=missing) .
#' @param family The distribution family to be used (Normal, T, or CN).
#' @param init Parameters specifies the initial values for model parameters.
#' @param thin An Interval at which samples are retained from the MCMC process to reduce autocorrelation.
#' @param chains The number of chains to run during the MCMC sampling. Running multiple chains is useful for checking convergence.
#' @param iter The total number of iterations for the MCMC sampling, determining how many samples will be drawn.
#' @param warmup The number of initial iterations that will be discarded as the algorithm stabilizes before collecting samples.

#' @return An object of class \code{HeckmanStan}, which is a list containing two elements:
#' \itemize{
#'   \item \code{list[[1]]}: Includes inference results from the Stan model, along with EAIC and EBIC.
#'   \item \code{list[[2]]}: Includes the HPC confidence intervals, along with LOOIC, WAIC, and CPO.
#' }
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
#' family="CN"
#' sigma2<- 1
#' rho<-0.7
#' beta<- c(1,0.5)
#' gamma<- c(1,0.3,-.5)
#' nu=c(0.1,0.1)
#' data<-geraHeckman(x,w,beta,gamma,sigma2,rho,nu,family=family)
#' y<-data$y
#' cc<-data$cc

#' # Fit Heckman Normal Stan model
#' fit.n_stan <- HeckmanStan(y, x, w, cc, family="Normal"
#'                          , thin = 5, chains = 1, iter = 10000, warmup = 1000)
#'qoi=c("beta","gamma","sigma_e","sigma2", "rho","EAIC","EBIC")
#'print(fit.n_stan[[1]],par=qoi)
#'print(fit.n_stan[[2]])
#'
#' require(rstan)
#'plot(fit.n_stan[[1]], pars=qoi)
#'plot(fit.n_stan[[1]], plotfun="hist", pars=qoi)
#'plot(fit.n_stan[[1]], plotfun="trace", pars=qoi)
#'plot(fit.n_stan[[1]], plotfun = "rhat")
#'
#' }
#'
#' @import loo
#' @importFrom rstan stan
#' @export
HeckmanStan <- function(y, x, w, cc, family="CN", init="random", thin = 5, chains = 1, iter = 10, warmup = 5){

  n<-length(cc)
  data = list(N = n, N_y = sum(cc==1), p = ncol(x), q = ncol(w), X = x[cc > 0, ], Z = w, D = cc, y = y[cc > 0])

  if (!family %in% c("Normal", "normal",
                     "T", "t",
                     "CN", "cn",
                     "SN", "sn")) {
    stop('Family not recognized! Only families allowed are: "Normal", "T", "CN" and "SN".')
  }
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


  if((family == "SN" || family == "sn")){

    stan_file <- system.file("stan", "HeckmanSkewNormal.stan",package = "HeckmanStan")
    out <- rstan::stan(file=stan_file,
                       data =data , init=init,
                       thin = thin, chains = chains, iter = iter, warmup = warmup)

  }

  paramT <- matrix(nrow = 0, ncol = 4)

  nsize<-out@par_dims$beta+out@par_dims$gamma+ifelse(is.null(out@par_dims$rho),0,1)+ ifelse(is.null(out@par_dims$sigma_e),0,1)+ ifelse(is.null(out@par_dims$sigma2),0,1)+ifelse(is.null(out@par_dims$nu),0,1)+ifelse(is.null(out@par_dims$nu1),0,1)+ifelse(is.null(out@par_dims$nu2),0,1)
  for(i in 1:nsize){
    pname<-names(out@sim$samples[[1]])
    param <- mean(out@sim$samples[[1]][[pname[i]]])
    se <- sd(out@sim$samples[[1]][[pname[i]]])
    HPDTot<-hpd(out@sim$samples[[1]][[pname[i]]], alpha=0.05)

    paramT <- rbind(paramT, c(param, se, HPDTot))
  }

##  dimnames(paramT) <- list(c(pname[1:nsize]),c("Mean", "Sd", " HPD(95%) Lower","Upper Bound"))

# fit_n <- extract(out, par=c("EAIC","EBIC"))
##  log_lik_n <- loo::extract_log_lik(out, merge_chains = FALSE)
##  loo_n <- loo::loo(log_lik_n)
##  WAIC_n = loo::waic(log_lik_n)
##  critFin<-cbind(loo_n$estimates[[3]],  WAIC_n$estimates[[3]], CPO(out))
# critFin<-cbind(mean(fit_n$EAIC), mean(fit_n$EBIC))
#c("EAIC","EBIC")
##  ll_total <- rowSums(loo::extract_log_lik(out))      # draw별 전체 로그우도

  #EAIC <- -2 * mean(ll_total) + 2 * nsize
  #EBIC <- -2 * mean(ll_total) + nsize * log()

##  dimnames(critFin)<-list(c("Value"),c("Looic", "WAIC","CPO"))
##  output<-list(paramT, critFin)

##  return(list(out, output ))
  return(out)
  }
