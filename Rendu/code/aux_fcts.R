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
  return(g)
}

hessian_old = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  h = matrix(0,ncol(x_sub),ncol(x_sub))
  
  for (i in 1:nrow(x_sub)) {
    h = h + u[i]*(1-u[i]) * t(x_sub) %*% x_sub
  }
  return(h)
}


hessian = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  h = matrix(0,ncol(x_sub),ncol(x_sub))
  txsub = t(x_sub)
  
  for (i in 1:nrow(x_sub)) {
    h = h + u[i]*(1-u[i]) * txsub %*% x_sub
  }
  return(h)
}

norm_vec <- function(x) sqrt(sum(x^2))


rhapson_newton_old <- function(x,y,betaO,eps) {
  beta = betaO
  dir = rep(eps+1,length(beta))
  
  while(norm_vec(dir)>eps)  {
    lglik = loglikl(x_sub = x,y = y,beta = beta)
    grd= gradient(x_sub = x,y = y,beta = beta)
    hes= hessian_old(x_sub = x,y = y,beta = beta)
    dir = solve(hes) %*% grd
    beta = beta - dir
  }
  return(beta)
}

rhapson_newton <- function(x,y,betaO,eps) {
  beta = betaO
  dir = rep(eps+1,length(beta))
  
  while(norm_vec(dir)>eps)  {
    lglik = loglikl(x_sub = x,y = y,beta = beta)
    grd= gradient(x_sub = x,y = y,beta = beta)
    hes= hessian(x_sub = x,y = y,beta = beta)
    dir = solve(hes) %*% grd
    beta = beta - dir
  }
  return(beta)
}
