#' @import mvtnorm

########################################
### CALCULO DE CPO
########################################
#' @importFrom rstan extract
CPO<-function(result)
{
  ll = data.frame(rstan::extract(result, pars = "log_invlik"))
  CPOb<-sum(log(1/apply(ll,2,mean)))
  return(CPOb)
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
