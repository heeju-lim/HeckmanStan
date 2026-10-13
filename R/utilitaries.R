#' @import mvtnorm

########################################
### log-sum-exp (numerically stable)
########################################
logSumExp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

########################################
### CALCULO DE CPO
########################################
# CPO_i = 1 / E_post[ 1 / f(y_i | theta) ]  (harmonic mean of the likelihood)
# Returns LPML = sum_i log(CPO_i)  (larger is better).
# Computed on the log scale from the S x N log_lik matrix, so exp(-log_lik)
# never overflows (the old version used log_invlik = exp(-log_lik), which is
# Inf for observations with very small likelihood).
#' @importFrom rstan extract
CPO <- function(result) {
  ll <- if (is.matrix(result)) result else loo::extract_log_lik(result, parameter_name = "log_lik")
  S <- nrow(ll)
  logCPO <- apply(ll, 2, function(l) log(S) - logSumExp(-l))
  sum(logCPO)
}

########################################
### Number of parameters of each family
########################################
# Used for EAIC / EBIC. Must match the penalty in the .stan files.
npar_heckman <- function(family, p, q) {
  switch(family,
         Normal = p + q + 2,   # beta, gamma, rho, sigma_e
         T      = p + q + 3,   # + nu
         CN     = p + q + 4,   # + nu1, nu2
         SN     = p + q + 3)   # beta, gamma, rho, sigma, lambda
}

########################################
### Model-selection criteria
########################################
# EAIC = -2 * E_post[log L] + 2k       (smaller is better)
# EBIC = -2 * E_post[log L] + k log(N) (smaller is better)
# LOOIC, WAIC from the 'loo' package    (smaller is better)
# CPO  = LPML = sum_i log CPO_i         (larger is better)
criteria_heckman <- function(out, k) {
  ll_arr <- loo::extract_log_lik(out, parameter_name = "log_lik", merge_chains = FALSE)
  ll <- loo::extract_log_lik(out, parameter_name = "log_lik", merge_chains = TRUE)
  N <- ncol(ll)
  ll_total <- rowSums(ll)

  EAIC <- -2 * mean(ll_total) + 2 * k
  EBIC <- -2 * mean(ll_total) + log(N) * k

  r_eff <- loo::relative_eff(exp(ll_arr))
  loo_fit  <- suppressWarnings(loo::loo(ll_arr, r_eff = r_eff))
  waic_fit <- suppressWarnings(loo::waic(ll))
  pk <- loo::pareto_k_values(loo_fit)
  nbad <- sum(pk > 0.7)
  if (nbad > 0)
    warning(nbad, " observation(s) with Pareto k > 0.7: LOOIC may be unreliable.", call. = FALSE)

  crit <- cbind(EAIC  = EAIC,
                EBIC  = EBIC,
                LOOIC = loo_fit$estimates["looic", "Estimate"],
                WAIC  = waic_fit$estimates["waic", "Estimate"],
                CPO   = CPO(ll))
  rownames(crit) <- "Value"
  list(criteria = crit, loo = loo_fit, waic = waic_fit)
}

########################################
### CALCULO hpd
########################################
hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]))
}
