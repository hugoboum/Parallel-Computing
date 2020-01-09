# Auxiliary functions

loglikl <- function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  
  l=0
  for (i in 1:nrow(x_sub)) {
    l = l - (y[i]*log(u[i]) + (1-y[i])*log(1-u[i]))
  }
  
  return(l)
}


gradient = function(x_sub,y,beta) {
  g= rep(0,ncol(x_sub))
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  
  for (i in 1:nrow(x_sub)) {
    g = g + x_sub[i,]*(u[i]-y[i])
  }
  return(as.matrix(g))
}

hessian = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  h = matrix(0,ncol(x_sub),ncol(x_sub))
  xtx = t(x_sub) %*% x_sub
  
  for (i in 1:nrow(x_sub)) {
    h = h + u[i]*(1-u[i]) * xtx 
  }
  return(h)
}

norm_vec <- function(x) sqrt(sum(x^2))


rhapson_newton <- function(x,y,betaO,its) {
  beta = betaO
  it=0
  while(it<its)  {
    lglik = loglikl(x_sub = x,y = y,beta = beta)
    grd= gradient(x_sub = x,y = y,beta = beta)
    hes= hessian(x_sub = x,y = y,beta = beta)
    dir = chol2inv(chol(hes)) %*% grd
    beta = beta - dir
    it = it + 1
  }
  return(beta)
}
