set.seed(123)

library(profvis)
library(MASS)

source('./datagen.R')
source('./aux_fcts.R')

n = 100
p=6
beta = c(10,5,1,0,0,0)
data = rlogit(nobs = n,beta = beta)
k=2

loglikl(x_sub = data$x[,1:2],y = data$y,beta = rep(1,2))
gradient(x_sub = data$x[,1:2],y = data$y,beta = rep(1,2))
hessian(x_sub = data$x[,1:2],y = data$y,beta = rep(1,2))

system.time({beta_s =  rhapson_newton(x = data$x[,1:k],y = data$y,betaO = rep(1,k),eps = 1e-6)})
beta_s

dat = data.frame(cbind(data$x[,1:k],data$y))
names = c(paste('x',1:k,sep = ''),'y')
colnames(dat) = names
system.time({g=glm(formula = y~.-1 ,family = binomial(link="logit") ,data = dat,start = rep(1,k))})
g$coefficients

p= profvis({
  beta_s =  rhapson_newton_old(x = data$x[,1:k],y = data$y,betaO = rep(1,k),eps = 1e-3)
})

htmltools::save_html(p,'run1a.html')

p= profvis({
  beta_s =  rhapson_newton(x = data$x[,1:k],y = data$y,betaO = rep(1,k),eps = 1e-3)
})

htmltools::save_html(p,'run1b.html')
