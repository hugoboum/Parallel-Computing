#Init

set.seed(123)

library(profvis)
library(MASS)

source('./datagen.R')
source('./aux_fcts_2b.R')

n = 100
beta = as.matrix(c(10,5,1,0,0,0))
data = rlogit(nobs = n,beta = beta)

#Rhapson-Newton algorithm with number of iterations as stopping criterion
k=2

system.time({beta_s =  rhapson_newton(x = data$x[,1:k],y = data$y,betaO = as.matrix(rep(1,k)),its = 4e4)})
beta_s

dat = data.frame(cbind(data$x[,1:k],data$y))
names = c(paste('x',1:k,sep = ''),'y')
colnames(dat) = names
system.time({g=glm(formula = y~.-1 ,family = binomial(link="logit") ,data = dat,start = rep(1,k))})
g$coefficients

#Cross-Validation

source('./basics.R')

p= profvis({
  cv.err = basic.cv(x = data$x[,1:k],y = data$y,betaO = as.matrix(rep(1,k)),its = 4e4,kf = 3)
})
htmltools::save_html(p,'run7a.html')
cv.err

#Model Comparison
##Create subsets of independant variables

s = cut(seq(ncol(data$x)),breaks = ncol(data$x) / k,labels = FALSE)
subsets = vector('list',length = (ncol(data$x) / k))

for (i in 1:length(subsets)) {
  subsets[[i]] = which(s==i)
}
subsets 

##Comparison

p= profvis({
  selected.model = basic.modelcomparison(x = data$x,models = subsets,y = data$y,betaO = as.matrix(rep(1,k)), its = 4e4 ,kf = 3)
})
htmltools::save_html(p,'run7b.html')
selected.model

#Model Selection

p= profvis({
  forward.model =  basic.modelselection(x = data$x,y = data$y,forward = TRUE,betaO = as.matrix(rep(1,ncol(data$x))), its = 4e4 ,kf = 3)
})
htmltools::save_html(p,'run7c.html')
forward.model

p= profvis({
  backward.model =  basic.modelselection(x = data$x,y = data$y,forward = FALSE,betaO = as.matrix(rep(1,ncol(data$x))), its = 4e4 ,kf = 3)
  
  })
htmltools::save_html(p,'run7d.html')
backward.model


#EXCEPTIONS
#Check the number of folds in basic.cv
cv.err = basic.cv(x = data$x[,1:k],y = data$y,betaO = as.matrix(rep(1,k)),its = 4e4,kf = n+1)

#Warn of LOOCV
cv.err = basic.cv(x = data$x[,1:k],y = data$y,betaO = as.matrix(rep(1,k)),its = 4e2,kf = n)

#We can check the type of models in basic.modelcomparison
selected.model = basic.modelcomparison(x = data$x,models = c(1:10),y = data$y,betaO = as.matrix(rep(1,k)), its = 4e4 ,kf = 3)

