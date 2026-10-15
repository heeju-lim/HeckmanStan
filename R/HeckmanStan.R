#' Fit the Heckman Selection Stan model using the Normal, Student-t, Contaminated Normal or Skew-Normal distributions.
#'
#' `HeckmanStan()` fits the Heckman selection model using a Bayesian approach to address sample selection bias.
#'
#' @param y A response vector.
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param cc A missing indicator vector (1=observed, 0=missing) .
#' @param family The distribution family to be used ("Normal", "T", "CN" or "SN").
#' @param init Parameters specifies the initial values for model parameters.
#' @param thin An Interval at which samples are retained from the MCMC process to reduce autocorrelation.
#' @param chains The number of chains to run during the MCMC sampling. Running multiple chains is useful for checking convergence.
#' @param iter The total number of iterations for the MCMC sampling, determining how many samples will be drawn.
#' @param warmup The number of initial iterations that will be discarded as the algorithm stabilizes before collecting samples.
#' @param ... Further arguments passed to \code{rstan::stan} (e.g. \code{seed}, \code{cores}, \code{control}).
#'
#' @return An object of class \code{HeckmanStan}, which is a list containing two elements:
#' \itemize{
#'   \item \code{list[[1]]}: The \code{stanfit} object (draws of all parameters, including EAIC and EBIC).
#'   \item \code{list[[2]]}: A list with \code{estimates} (posterior mean, sd and 95\% HPD interval) and
#'   \code{criteria} (EAIC, EBIC, LOOIC, WAIC and CPO).
#' }
#' EAIC, EBIC, LOOIC and WAIC: smaller is better. CPO is reported as the log pseudo-marginal
#' likelihood, LPML = sum of log(CPO_i): larger is better.
#' The full \code{loo} and \code{waic} objects are stored in \code{attr(fit, "loo")} and \code{attr(fit, "waic")}.
#'
#' @examples
#' \donttest{
#' ################################################################################
#' # Simulation
#' ################################################################################
#' library(mvtnorm)
#' set.seed(1)
#' n<- 300
#' w<- cbind(1,rnorm(n),rnorm(n))
#' x<- cbind(w[,1:2])
#' sigma2<- 1
#' rho<-0.7
#' beta<- c(1,0.5)
#' gamma<- c(1,0.3,-.5)
#' lambda<- 2
#' data<-geraHeckman(x,w,beta,gamma,sigma2,rho,nu=lambda,family="SN")
#' y<-data$y
#' cc<-data$cc
#'
#' # Fit Heckman Skew-Normal Stan model
#' fit.sn <- HeckmanStan(y, x, w, cc, family="SN",
#'                       chains = 2, iter = 2000, warmup = 1000)
#' fit.sn                 # estimates + EAIC, EBIC, LOOIC, WAIC, CPO
#' fit.sn[[2]]$criteria
#'
#' qoi=c("beta","gamma","sigma","lambda","rho","EAIC","EBIC")
#' print(fit.sn[[1]],par=qoi)
#'
#' require(rstan)
#'plot(fit.sn[[1]], pars=qoi[1:5])
#'plot(fit.sn[[1]], plotfun="trace", pars=qoi[1:5])
#'
#' }
#'
#' @import loo
#' @importFrom rstan stan
#' @export
HeckmanStan <- function(y, x, w, cc, family="CN", init="random", thin = 1, chains = 1, iter = 2000, warmup = 1000, ...){

  fam <- switch(tolower(family), normal = "Normal", t = "T", cn = "CN", sn = "SN", NA)
  if (is.na(fam))
    stop('Family not recognized! Only families allowed are: "Normal", "T", "CN" and "SN".')

  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.data.frame(w)) w <- as.matrix(w)
  if(is.vector(x)) x <- as.matrix(x)
  if(is.vector(w)) w <- as.matrix(w)
  if(!is.matrix(x)) stop("x must be a matrix!")
  if(!is.matrix(w)) stop("w must be a matrix!")
  y  <- as.vector(y)
  cc <- as.vector(cc)
  if(!all(cc %in% c(0, 1))) stop("cc must contain only 0 and 1!")
  if(length(y) != length(cc) || nrow(x) != length(cc) || nrow(w) != length(cc))
    stop("y, x, w and cc must have the same number of observations!")
  if(anyNA(y[cc == 1]) || anyNA(x[cc == 1, ]) || anyNA(w))
    stop("Missing values (NA) are only allowed in y where cc == 0!")

  n <- length(cc)
  data <- list(N = n, N_y = sum(cc == 1), p = ncol(x), q = ncol(w),
               X = x[cc == 1, , drop = FALSE], Z = w, D = as.integer(cc),
               y = y[cc == 1])

  stan_name <- switch(fam,
                      Normal = "HeckmanNormal.stan",
                      T      = "HeckmanT.stan",
                      CN     = "HeckmanCNormal.stan",
                      SN     = "HeckmanSkewNormal_centered.stan")
  stan_file <- system.file("stan", stan_name, package = "HeckmanStan")
  out <- rstan::stan(file = stan_file, data = data, init = init,
                     thin = thin, chains = chains, iter = iter, warmup = warmup, ...)
  if (out@mode != 0L || length(out@sim) == 0L)
    stop("Stan sampling failed; see the messages above.")

  ## ---- posterior summary: mean, sd, 95% HPD (all chains, warmup excluded) ----
  scale_par <- if (fam == "SN") "sigma" else "sigma_e"
  extra_par <- switch(fam, Normal = NULL, T = "nu", CN = c("nu1", "nu2"), SN = "lambda")
  pars <- c("beta", "gamma", "rho", scale_par, "sigma2", extra_par)

  draws <- as.matrix(out, pars = pars)
  paramT <- t(apply(draws, 2, function(d) c(mean(d), sd(d), hpd(d, alpha = 0.05))))
  colnames(paramT) <- c("Mean", "Sd", "HPD(95%) Lower", "HPD(95%) Upper")

  ## ---- model-selection criteria: EAIC, EBIC, LOOIC, WAIC, CPO ----
  k <- npar_heckman(fam, data$p, data$q)
  crit <- criteria_heckman(out, k)

  output <- list(paramT, crit$criteria)
  names(output) <- c("estimates", "criteria")

  res <- list(out, output)
  attr(res, "loo")  <- crit$loo
  attr(res, "waic") <- crit$waic
  class(res) <- c("HeckmanStan", "list")
  return(res)
}

#' @method print HeckmanStan
#' @export
print.HeckmanStan <- function(x, digits = 4, ...) {
  cat("Heckman selection model (Stan)\n\n")
  cat("Posterior estimates:\n")
  print(round(x[[2]]$estimates, digits))
  cat("\nModel selection criteria:\n")
  cat("(EAIC, EBIC, LOOIC, WAIC: smaller is better; CPO = LPML: larger is better)\n")
  print(round(x[[2]]$criteria, digits))
  invisible(x)
}
