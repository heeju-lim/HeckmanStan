#' Generating Heckman data  : Normal, Student-t, Slash and Laplace
#'
#' `geraHeckman()` generates a random sample from the Heckman selection model (Normal, Student-t or CN).
#'
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param beta Values for the beta vector.
#' @param gamma Values for the gamma vector.
#' @param sigma2 Value for the variance.
#' @param rho Value for the dependence between the response and missing value.
#' @param nu When using the t- distribution, the initial value for the degrees of freedom.
#' @param family The distribution family to be used (Normal, T, or CN).
#' @return Return an object with the response (y) and missing values (cc).
#'
#' @examples
#' \donttest{
#' n <- 100
#' rho <- .6
#' cens <- 0.25
#' nu <- 4
#' set.seed(20200527)
#' w <- cbind(1,runif(n,-1,1),rnorm(n))
#' x <- cbind(w[,1:2])
#'
#' family <- "T"
#' c <- qt(cens, df=nu)
#'
#' sigma2 <- 1
#' beta <- c(1,0.5)
#' gamma<- c(1,0.3,-.5)
#' gamma[1] <- -c*sqrt(sigma2)
#'
#' data <- geraHeckman(x,w,beta,gamma,sigma2,rho,nu,family=family)
#' }
#' @import stats
#' @export
geraHeckman<-function(x,w,beta,gamma,sigma2,rho,nu,family="T"){
n<-nrow(x)
rhoa<- rho*sqrt(sigma2)

if(family=="Normal"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)
resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}

if(family=="T"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
errorTerms<- mvtnorm::rmvt(n, Sigma,df=nu)
resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}


if(family=="SL"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)*1/sqrt(rbeta(n,nu,1))
resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}


if(family=="La"){#Laplace
  Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
  errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)*sqrt(rgamma(n,1.5,0.125))
  resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}

if(family=="CN"){#Contaminated Normal
  p <- runif(n)
  u <- rep(1,n)
  u[p<nu[1]] <- nu[2]
  Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
  errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)/sqrt(u)
  resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}

cc<-(resp[,2]> 0)+0
resp[cc==0,1]<-0

return=list(y=resp[,1],cc=cc)
}
#'
#' @references
#'  Lachos, V. H., Prates, M. O., & Dey, D. K. (2021).
#'  Heckman selection-t model: Parameter estimation via the EM-algorithm.
#'  Journal of Multivariate Analysis, 184, 104737.
#' @name geraHeckman
NULL
