#' @import mvtnorm
################################################################################
### likelihood function  : Normal
################################################################################
likeL<-function(y, x, w, cc, beta, gama, Sigma, rho){
sigma2<- Sigma[1,1]
sigma<- sqrt(sigma2)
rho<- Sigma[1,2]/sigma
med<-(y-x%*%beta)/sigma
#fd<- cc*log(dnorm(med)*pnorm((rho*med+w%*%gama)/sqrt(1-rho^2))/sigma)+(1-cc)*log(pnorm(-w%*%gama))
fd <- ifelse(cc==1, dnorm(med,log=T)+pnorm((rho*med+w%*%gama)/sqrt(1-rho^2), log=T)-log(sigma), pnorm(-w%*%gama,log=T))
return(sum(fd))
}

################################################################################
### likelihood function  : Normal   used to direct maximization
################################################################################
likeLmax<-function(theta, y, x, w, cc){
p<-ncol(x)
qq<-ncol(w)
beta<-as.matrix(theta[1:p],p,1)
gama<-as.matrix(theta[(p+1):(p+qq)],qq,1)
sigma<-theta[p+qq+1]
rho<- theta[p+qq+2]
med<-(y-x%*%beta)/sigma
#fd<- cc*log(dnorm(med)*pnorm((rho*med+w%*%gama)/sqrt(1-rho^2))/sigma)+(1-cc)*log(pnorm(-w%*%gama))
fd <- ifelse(cc==1, dnorm(med,log=T)+pnorm((rho*med+w%*%gama)/sqrt(1-rho^2), log=T)-log(sigma), pnorm(-w%*%gama,log=T))
return(-sum(fd))
}

################################################################################
### likelihood function  : Student--t
################################################################################
likeLt<-function(nu, y, x, w, cc, beta, gama, Sigma, rho){
sigma2<- Sigma[1,1]
sigma<- sqrt(sigma2)
rho<- Sigma[1,2]/sigma
med<-(y-x%*%beta)/sigma
aux1<-sqrt((nu+1)/(nu+med^2))*((rho*med+w%*%gama)/(sqrt(1-rho^2)))
#fd<-cc*log(dt(med,nu)*pt(aux1,nu+1)/sigma)+(1-cc)*log(pt(-w%*%gama,nu))
fd<-ifelse(cc==1, dt(med,nu,log=T)+pt(aux1,nu+1,log=T)-log(sigma), pt(-w%*%gama,nu,log=T))
return(sum(fd))
}

################################################################################
### likelihood function  : Student--t
################################################################################
likeLt2<-function(nu, y, x, w, cc, beta, gama, Sigma, rho){
nu2 <- exp(nu) + 3.001
sigma2<- Sigma[1,1]
sigma<- sqrt(sigma2)
rho<- Sigma[1,2]/sigma
med<-(y-x%*%beta)/sigma
aux1<-sqrt((nu2+1)/(nu2+med^2))*((rho*med+w%*%gama)/(sqrt(1-rho^2)))
#fd<-cc*log(dt(med,nu2)*pt(aux1,nu2+1)/sigma)+(1-cc)*log(pt(-w%*%gama,nu2))
fd<-ifelse(cc==1, dt(med,nu2,log=T)+pt(aux1,nu2+1,log=T)-log(sigma), pt(-w%*%gama,nu2,log=T))
return(-sum(fd))
}

################################################################################
### likelihood function  : Student--t   used to direct maximization
################################################################################
liketLmaxnu<-function(theta, y, x, w, cc){
p<-ncol(x)
qq<-ncol(w)
beta<-as.matrix(theta[1:p],p,1)
gama<-as.matrix(theta[(p+1):(p+qq)],qq,1)
sigma<-theta[p+qq+1]
rho<- theta[p+qq+2]
nu<- theta[p+qq+3]
med<-(y-x%*%beta)/sigma
aux1<-sqrt((nu+1)/(nu+med^2))*((rho*med+w%*%gama)/(sqrt(1-rho^2)))
#fd<-cc*log(dt(med,nu)*pt(aux1,nu+1)/sigma)+(1-cc)*log(pt(-w%*%gama,nu))
fd<-ifelse(cc==1, dt(med,nu,log=T)+pt(aux1,nu+1,log=T)-log(sigma), pt(-w%*%gama,nu,log=T))
return(-sum(fd))
}

################################################################################
### likelihood function  : Laplace   used to direct maximization
################################################################################

likeLapmax<-function(theta, y, x, w, cc){
n<-nrow(x)
p<-ncol(x)
qq<-ncol(w)
beta<-as.matrix(theta[1:p],p,1)
gama<-as.matrix(theta[(p+1):(p+qq)],qq,1)
sigma<-theta[p+qq+1]
rho<- theta[p+qq+2]

xi<-0.00000001
psi<-0.25
mux<-x%*%beta
muw<-w%*%gama
mu1.2<-  muw+rho/sigma*(y-mux)
sigma1.2<-1-rho^2
delta1<- (y-mux)^2/sigma^2

fd<-matrix(0,n,1)

for (j in 1:n){
paraX<- ghyp(lambda = 1.5, chi = xi, psi = psi, mu = mux[j], sigma = sigma^2,
gamma = 0)
paraW<- ghyp(lambda = 1.5, chi = xi, psi = psi, mu = muw[j], sigma = 1 ,
gamma = 0)
paraC<-ghyp(lambda = 1.5, chi = xi, psi = psi, mu = mu1.2[j], sigma = sigma1.2,
gamma = 0)

aux1<- pghyp(0, object = paraC,  n.sim = 1000, subdivisions =200 , lower.tail = FALSE)
 if(aux1 == 0) aux1 <- .Machine$double.xmin
aux2<- pghyp(0, object = paraW,  n.sim = 1000, subdivisions = 200)
 if(aux2 == 0) aux2 <- .Machine$double.xmin
fd[j]<-ifelse(cc[j]==1, dghyp(y[j], object = paraX, logvalue = TRUE)
+log(aux1),
 log(aux2))
}
return(-sum(fd))
}

################################################################################
## likeLapmax2: Corrigido 01/26/2024
################################################################################
likeLapmax2<-function(theta, y, x, w, cc){
n<-nrow(x)
p<-ncol(x)
qq<-ncol(w)
beta<-as.matrix(theta[1:p],p,1)
gama<-as.matrix(theta[(p+1):(p+qq)],qq,1)
sigma<-theta[p+qq+1]
rho<- theta[p+qq+2]
#rho<-(exp(2*arcrho)-1)/(exp(2*arcrho)+1)

xi<- 0.00000001
psi<- 0.25
mux<- x%*%beta
muw<- w%*%gama
mu1.2<-  muw+rho/sigma*(y-mux)
sigma1.2<- 1-rho^2

paraX<- ghyp(lambda = 1, chi = xi, psi = psi, mu = 0, sigma = sigma^2,
gamma = 0)

Auxdp<-dghyp(y-mux, object = paraX, logvalue = TRUE)

paraW<- ghyp(lambda = 1, chi = xi, psi = psi, mu = 0, sigma = 1 ,
gamma = 0)

aux2<- pghyp(rep(0,n)-muw, object = paraW,  n.sim = 10000, subdivisions = 100)
 if(length(which(aux2 == 0)) > 0) aux2[which(aux2 == 0)] <- .Machine$double.xmin

fd<- matrix(0,n,1)
aux1<- matrix(0,n,1)

for (j in 1:n){
if(cc[j]==1){
paraC<-ghyp(lambda = 1, chi = xi, psi = psi, mu = mu1.2[j], sigma = sigma1.2,
gamma = 0)
aux1[j]<- pghyp(0, object = paraC,  n.sim = 10000, subdivisions =100 , lower.tail = FALSE)
 if(aux1[j] == 0) aux1[j] <- .Machine$double.xmin
 }
fd[j]<-ifelse(cc[j]==1, Auxdp[j]+log(aux1[j]), log(aux2[j]))
}
return(-sum(fd))
}

##test
#n<-1000
#sigma2<- 1
#rho=0.8
#beta<- c(1,0.5)
#gama<- c(1,0.3,-.5)
#w<- cbind(1,runif(n,-1,1),rnorm(n))
#x<-cbind(w[,1:2])
#type <- "T"
#datas<-geraHeckman(x,w,beta,gama,sigma2,rho,nu,type=type)
#y<-datas$y
#cc<-datas$cc
#theta<-c(beta,gama,sigma2,rho)
#likeLapmax(theta, y, x, w, cc)
#likeLapmax2(theta, y, x, w, cc)


### Sem nu
liketLmax<-function(theta, y, x, w, cc, nu){
p<-ncol(x)
qq<-ncol(w)
beta<-as.matrix(theta[1:p],p,1)
gama<-as.matrix(theta[(p+1):(p+qq)],qq,1)
sigma<-theta[p+qq+1]
rho<- theta[p+qq+2]
#nu<- theta[p+qq+3]
med<-(y-x%*%beta)/sigma
aux1<-sqrt((nu+1)/(nu+med^2))*((rho*med+w%*%gama)/(sqrt(1-rho^2)))
fd<-cc*log(dt(med,nu)*pt(aux1,nu+1)/sigma)+(1-cc)*log(pt(-w%*%gama,nu))
return(-sum(fd))
}


################################################################################
#### Generating Heckman data  : Normal, Student-t, Slash and Laplace
################################################################################

geraHeckman<-function(x,w,beta,gama,sigma2,rho,nu,type="T"){
n<-nrow(x)
rhoa<- rho*sqrt(sigma2)

if(type=="Normal"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)
resp<- cbind(x%*%beta,w%*%gama)+ errorTerms
}

if(type=="T"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
errorTerms<- mvtnorm::rmvt(n, Sigma,df=nu)
resp<- cbind(x%*%beta,w%*%gama)+ errorTerms
}


if(type=="SL"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)*1/sqrt(rbeta(n,nu,1))
resp<- cbind(x%*%beta,w%*%gama)+ errorTerms
}


if(type=="La"){#Laplace
  Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
  errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)*sqrt(rgamma(n,1.5,0.125))
  resp<- cbind(x%*%beta,w%*%gama)+ errorTerms
}

if(type=="CN"){#Contaminated Normal
  p <- runif(n)
  u <- rep(1,n)
  u[p<nu[1]] <- nu[2]
  Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
  errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)/sqrt(u)
  resp<- cbind(x%*%beta,w%*%gama)+ errorTerms
}

cc<-(resp[,2]> 0)+0
resp[cc==0,1]<-0

return=list(y=resp[,1],cc=cc)
}





################################################################################
# Algorithm Laplace Optim Parallel
################################################################################

LaplaceDM.alg<-function(y, x, w, cc){

n<-nrow(x)
y<-matrix(y,n,1)
p<-ncol(x)
q<-ncol(w)

ResN<- HeckmanEM(c(y), x, w, cc, nu = NULL, family="Normal", error = 1e-05,iter.max = 100, im=FALSE, criteria = TRUE)

beta<-ResN$beta
gama<-ResN$gamma
sigma<- ResN$sigma
rho<- ResN$rho

#cl <- makeCluster(2) # set the number of processor cores
cl <- parallel::makeCluster(parallel::detectCores()-1)
setDefaultCluster(cl=cl) # set 'cl' as default cluster

clusterExport(cl,"ghyp")
clusterExport(cl,"dghyp")
clusterExport(cl,"pghyp")



OpL<-try(optimParallel(c(beta,gama,sigma,rho), likeLapmax2,
      lower = c(rep(-Inf, ncol(x)),rep(-Inf, ncol(w)),0.1,-0.99), upper = c(rep(Inf, ncol(x)),rep(Inf, ncol(w)),Inf,0.99),
      y=y, x=x, w=w, cc=cc,hessian=TRUE))

setDefaultCluster(cl=NULL)
stopCluster(cl)

beta<-OpL$par[1:p]
gama<-OpL$par[(p+1):(p+q)]
sigma<-OpL$par[p+q+1]
rho<-OpL$par[p+q+2]

lk<- -OpL$value
aic <- -2*lk + 2*(p+q+2)
bic <- -2*lk + log(n)*(p+q+2)
desvios <- sqrt(diag(solve(OpL$hessian)))
 return(list(beta=beta, gamma=gama, rho=rho,sigma=sigma, sd=desvios, logL=lk, AIC=aic,BIC=bic))
}

################################################################################
# Algorithm Laplace without Optim Parallel
################################################################################

LaplaceDM.algN<-function(y, x, w, cc){

n<-nrow(x)
y<-matrix(y,n,1)
p<-ncol(x)
q<-ncol(w)

ResN<- HeckmanEM(c(y), x, w, cc, nu = 3, family="T", error = 1e-05,iter.max = 1000, im=FALSE, criteria = TRUE)

beta<-ResN$beta
gama<-ResN$gamma
sigma<- ResN$sigma
arcrho<- atanh(ResN$rho)


OpL<-try(optim(c(beta,gama,sigma,arcrho), likeLapmax2, method = "L-BFGS-B",
      lower = c(rep(-Inf, ncol(x)),rep(-Inf, ncol(w)),0.00001,-8), upper = c(rep(Inf, ncol(x)),rep(Inf, ncol(w)),Inf,8),
      y=y, x=x, w=w, cc=cc,hessian=TRUE,control = list(maxit = 30000, temp = 2000, trace = TRUE,
                            REPORT = 500)))

beta<-OpL$par[1:p]
gama<-OpL$par[(p+1):(p+q)]
sigma<-OpL$par[p+q+1]
arcrho<-OpL$par[p+q+2]
rho<-(exp(2*arcrho)-1)/(exp(2*arcrho)+1)
lk<- -OpL$value
aic <- -2*lk + 2*(p+q+2)
bic <- -2*lk + log(n)*(p+q+2)
desvios <- sqrt(diag(solve(OpL$hessian)))
return(list(beta=beta, gamma=gama, rho=rho,sigma=sigma,arcrho=arcrho, sd=desvios, logL=lk, AIC=aic,BIC=bic))
}

# atanh(0.9999999999999999)

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
