# Data generation

rlogit  <- function(nobs,beta) {
  
  nvars = length(beta)
  
  covs1 = matrix(rnorm(n = nobs*round(nvars/2),mean = 1,sd = 2),nrow = nobs,ncol = round(nvars/2)) 
  covs2 = matrix(rnorm(n = nobs*round(nvars/2),mean = 0,sd = 1),nrow = nobs,ncol =  nvars - round(nvars/2)) 
  
  x = cbind(covs1,covs2)
  
  u = exp(x %*% beta) / (1 + exp(x %*% beta))
  
  y = as.numeric(runif(nobs) < u)
  return(list('x'=x,'y'=y,'u'=u))
  
}


