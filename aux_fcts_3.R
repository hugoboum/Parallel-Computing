# Auxiliary functions

library(Rcpp)
library(inline)

matmut <- cxxfunction(signature(tm="NumericMatrix",
                               tm2="NumericMatrix"),
                     plugin="RcppEigen",
                     body="
NumericMatrix tm22(tm2);
NumericMatrix tmm(tm);

const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));

Eigen::MatrixXd prod = ttm*ttm2;
return(wrap(prod));
                 ")


loglikl <- function(x_sub,y,beta) {
  u = exp(matmut(x_sub,beta)) / (1 + exp(matmut(x_sub,beta)))
  vec = sapply(1:nrow(x_sub),function(i){y[i]*log(u[i]) + (1-y[i])*log(1-u[i])})
  return(-sum(vec))
}


gradient = function(x_sub,y,beta) {
  u = exp(matmut(x_sub,beta)) / (1 + exp(matmut(x_sub,beta)))
  vec = sapply(1:nrow(x_sub),function(i){x_sub[i,]*(u[i]-y[i])})
  return(as.matrix(rowSums(vec)))
}

hessian = function(x_sub,y,beta) {
  u = exp(matmut(x_sub,beta)) / (1 + exp(matmut(x_sub,beta)))
  xtx = matmut(t(x_sub),x_sub)
  vec = sapply(1:nrow(x_sub),function(i){u[i]*(1-u[i]) * xtx})
  h = matrix(rowSums(vec),ncol(x_sub),ncol(x_sub))
  return(h)
}

norm_vec <- function(x) sqrt(sum(x^2))
  

rhapson_newton <- function(x,y,betaO,eps) {
  beta = betaO
  dir = rep(eps+1,length(beta))
  
  while(norm_vec(dir)>eps)  {
    lglik = loglikl(x_sub = x,y = y,beta = beta)
    grd= gradient(x_sub = x,y = y,beta = beta)
    hes= hessian(x_sub = x,y = y,beta = beta)
    dir = matmut(solve(hes),grd)
    beta = beta - dir
  }
  return(beta)
}
