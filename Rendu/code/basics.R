source('./aux_fcts_2b.R')

basic.cv = function(x,y,betaO,its,kf){
  #kf is the number of cv folds
  if (kf>nrow(x)){stop('Cannot partition the data with this much folds')}
  #Here we could shuffle the data if necessary
  #
  #LOOCV
  if (kf == nrow(x)){message('Performing LOOCV')}
  
  folds = cut(seq(length(y)),breaks = kf,labels = FALSE)
  errs = numeric(kf)
  
  for (i in 1:kf) {
    train.ind = which(folds != i,arr.ind = TRUE)
    beta_s =  rhapson_newton(x = as.matrix(x[train.ind,]),y = y[train.ind],betaO = betaO, its = its)
    preds = exp(x[-train.ind,] %*% beta_s) / (1 + exp(x[-train.ind,] %*% beta_s))
    preds = as.numeric(preds >= 0.5)
    errs[i] = sum((preds - y[-train.ind])**2) / length(preds)
  }
  return(mean(errs))
}

basic.modelcomparison <- function(x,models,y,betaO,its,kf) {
  #x should be a matrix of independant variables
  #models should be a list, of which each element denotes the subset of selected variables (column indices of x)
  
  if(typeof(models) != "list"){stop('models should be a list, of which each element denotes the subset of selected variables (column indices of x)')}
  
  models.errs = numeric(length(models))
  
  for (i in 1:length(models)) {
    models.errs[i] = basic.cv(x = as.matrix(x[,models[[i]]]),y = y,betaO = betaO,its = its,kf = kf)
  }
  return(list('vars'=models[[which.min(models.errs)]],'cv.err'= models.errs[which.min(models.errs)]))
}

basic.modelselection <- function(x,y,forward,betaO,its,kf) {
  # forward should be a boolean
  
  d_cv.err = 1
  t=0
  subsets = vector('list',length = ncol(x))
  
  if (forward==TRUE) {
    selected.model = 0
    cv.err = basic.cv(x =as.matrix(rep(1,length(y))) ,y = y,betaO = betaO[1] ,its = its,kf=kf)
    for (i in 1:length(subsets)) {
      subsets[[i]] = i
    }
    while ((d_cv.err > 0) & (t < ncol(x))) {
      #Save previous model
      m.selected.model = selected.model
      if (t>0) {
        m.selected.model$vars = sort(m.selected.model$vars)
      }
      #Compare models
      selected.model = basic.modelcomparison(x = x,models = subsets,y = y,betaO = as.matrix(betaO[1:(t+1)]),its = its,kf = kf) 
      #Update errors
      d_cv.err = cv.err - selected.model$cv.err
      cv.err = selected.model$cv.err
      #Update models
      for (i in 1:length(subsets)) {
        subsets[[i]] = c(subsets[[i]],selected.model$vars[1])
      }
      #Delete subsets with duplicates variables
      for (i in length(subsets):1) {
        if (length(subsets[[i]]) > length(unique(subsets[[i]]))) {
          subsets[[i]]=NULL
        }
      }
      t=t+1
    }
  }else{
    #BACKWARD
    selected.model = c(1:ncol(x))
    cv.err = basic.cv(x =x ,y = y,betaO = betaO ,its = its,kf=kf)
    for (i in 1:length(subsets)) {
      subsets[[i]] = c(1:length(subsets))[-i]
    }
    while ((d_cv.err >= 0) & (t < ncol(x))) {
      #Save previous model
      m.selected.model = selected.model
      if (t>0) {
        m.selected.model$vars = sort(m.selected.model$vars)
      }
      #Compare models
      selected.model = basic.modelcomparison(x = x,models = subsets,y = y,betaO = as.matrix(betaO[1:(length(betaO)-(t+1))]),its = its,kf = kf) 
      #Update errors
      d_cv.err = cv.err - selected.model$cv.err
      cv.err = selected.model$cv.err
      #Update models
      for (i in 1:(length(subsets)-1)) {
        subsets[[i]] = selected.model$vars[-i]
      }
      #Delete residual
      subsets[[length(subsets)]]=NULL 
      t=t+1
    }
  }
  #Remove last added,not useful,variable.Return instead the previous to last model
  return(m.selected.model)
}