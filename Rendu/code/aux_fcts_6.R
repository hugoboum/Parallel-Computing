# Auxiliary functions

loglikl <- function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  vec = mclapply(1:nrow(x_sub),function(i){y[i]*log(u[i]) + (1-y[i])*log(1-u[i])},mc.cores = 12)
  l = do.call(sum,vec)
  return(-l)
}

gradient = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  vec = mclapply(1:nrow(x_sub),function(i){x_sub[i,]*(u[i]-y[i])},mc.cores = 12)
  vec = rowSums(matrix(unlist(vec),ncol(x_sub)))
  return(as.matrix(vec))
}

hessian = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  xtx = t(x_sub) %*% x_sub
  vec = mclapply(1:nrow(x_sub),function(i){u[i]*(1-u[i]) * xtx},mc.cores = 12)
  vec = rowSums(matrix(unlist(vec),ncol(x_sub)**2))
  return(matrix(vec,ncol(x_sub),ncol(x_sub)))
}

norm_vec <- function(x) sqrt(sum(x^2))
  
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
