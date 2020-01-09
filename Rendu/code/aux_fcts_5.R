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

src <- '
   Rcpp::List input(data);
   Rcpp::Function f(fun);
   Rcpp::List output(input.size());
   std::transform(input.begin(), input.end(), output.begin(), f);
   output.names() = input.names();
   return output;
   '
cpp_lapply <- cxxfunction(signature(data = "list", fun = "function"),src, plugin = "Rcpp")


require(parallel)
detectCores()

#system.time({means = do.call(rbind,mclapply(data,colMeans,mc.cores = 12))})



loglikl <- function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  vec = cpp_lapply(1:nrow(x_sub),function(i){y[i]*log(u[i]) + (1-y[i])*log(1-u[i])})
  l = do.call(sum,vec)
  return(-l)
}


gradient = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  vec = cpp_lapply(1:nrow(x_sub),function(i){x_sub[i,]*(u[i]-y[i])})
  vec = rowSums(matrix(unlist(vec),ncol(x_sub)))
  return(as.matrix(vec))
}

hessian = function(x_sub,y,beta) {
  u = exp(x_sub %*% beta) / (1 + exp(x_sub %*% beta))
  xtx = t(x_sub) %*% x_sub
  vec = cpp_lapply(1:nrow(x_sub),function(i){u[i]*(1-u[i]) * xtx})
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
