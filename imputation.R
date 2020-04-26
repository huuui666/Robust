
## Functions for multiple imputation via iterative model-based imputation


# Multiple imputation
# 
# Input:
# xy    data set with missing values
# m     number of imputations
# DDC   a logical indicating whether to run the DetectDeviatingCells algorithm
#       before performing multiple imputation (such that the flagged outlying 
#       cells are also imputed)
# ...   additional arguments to be passed to function irmi() from package VIM
#       (for example, whether to use robust models)
# 
# Output:
# A list with the following components:
# imputed   this should again be a list in which each list element is one 
#           imputed data set
# m         the number of imputations


multimp <- function(xy, m, DDC = FALSE, ...) {
  # You should set a sensible default for the number of imputations m.
  #
  # You can use function DDC() from package cellWise, as well as function 
  # irmi() from package VIM for imputations.
  library(VIM)
  library(cellWise)
  xy = as.matrix(xy)
  if (DDC == TRUE){
    ddc_result = DDC(xy)
    indice = ddc_result$indcells
    xy[indice] = NA
  }
  if (missing(m)){
    missing_ob = sum(!complete.cases(xy))
    missing_percentage = missing_ob / nrow(xy)
    m = ceiling(missing_percentage*100)
  }
  
  imputed= irmi(xy,mi=m, robust=TRUE)
  output = list('imputed'= imputed, 'm' = m,'data'=xy)
}


# Fit regression models
# 
# Input:
# xyList   list of imputed data sets as returned by function multimp()
# ...      additional arguments to be passed to modeling function (for example,
#          control parameters for the MM-algorithm)
# 
# Output:
# A list with the following components:
# models   this should again be a list in which each list element is a 
#          regression model (fitted to the corresponding imputed data set)
# m        the number of imputations

fit <- function(xyList, ...) {
  #  use function lmrob() from package robustbase for the MM-estimator.
  library(robustbase)
  m = xyList$m
  data = xyList$imputed
  models = list()
  for (i in 1:m){
    models[[i]] = lmrob(data[[i]],control=lmrob.control(max.it=1000,k.max=1000,maxit.scale = 1000))
  
  }
  output = list('models' = models, 'm' = m)
}


# Pool point estimates and standard errors
#
# Input:
# fitList  a list as returned by function fit() containting the regression 
#          models fitted to the imputed data sets
# ...      additional arguments you may need to pass down to other functions
#
# Output:
# A matrix that contains the pooled coefficients in the first column, the
# pooled standard errors in the second column, the t-statistic in the third
# column, the estimated degrees of freedom in the fourth column, and the
# p-value in the fifth column 

pool <- function(fitList, ...) {
  m = fitList$m
  models = fitList$models
  p = length(models[[1]]$coefficients)
  coe = matrix(NA, p,m)
  for (i in 1:m){
    coe[,i] = models[[i]]$coefficients
  }
  pool_coe = apply(coe,1,mean)
  
  between_variance = apply(coe,1,var)
  
  sd =matrix(NA, p,m)
  for (i in 1:m){
    su = summary(models[[i]])
    sd[,i] = (su$coefficients[,2])^2
  }
  within_variance = apply(sd,1,mean)
  
  pool_var = within_variance + (m+1)/m *between_variance
  pool_sd = sqrt(pool_var)
  
  pool_t = pool_coe/ pool_sd
  
  gamma_m = (m+1)/m * between_variance/pool_var
  v_m = (m-1)/gamma_m^2
  vcomp = nrow(models[[1]]$model)-p
  v_obs = (vcomp+1)/(vcomp+3) * vcomp * (1-gamma_m)
  pool_freedom = v_m*v_obs / (v_m+v_obs)
  
  pool_p_value = matrix(NA,p,1)
  for (i in 1:p){
    pool_p_value[i,] = 2*pt(abs(pool_t[i]), pool_freedom[i], lower.tail = FALSE)
  }
  output_matrix = cbind(pool_coe,pool_sd,pool_t,pool_freedom,pool_p_value)
  colnames(output_matrix) = c('Estimate','Std.Error','t.values','df','Pr(>|t|)')
  return(output_matrix)
  
  
}



## Function for the bootstrap with kNN imputation and linear regression


# Input:
#
# x     data set with missing values
# R     number of bootstrap replications
# k     number of neighbors for kNN imputation
# DDC   a logical indicating whether to run the DetectDeviatingCells algorithm
#       before performing imputation (such that the flagged outlying cells are 
#       also imputed)
# ...   additional arguments to be passed to modeling function (for example,
#       control parameters for the MM-algorithm)
#
# Output:
# A list with the following components:
# replicates   a matrix with all coefficient estimates from all replications
# summary      a matrix that contains the point estimates of the coefficients 
#              in the first column, the standard errors in the second column, 
#              the z-statistic in the third column, and the p-value in the 
#              fourth column 

bootstrap <- function(x, R=500, k=5, DDC = FALSE, ...) {

  #  use function DDC() from package cellWise, function kNN() from 
  # package VIM for imputations, and function lmrob() from package robustbase 
  # for the MM-estimator.
  library(cellWise)
  library(VIM)
  x = as.matrix(x)
  if (DDC == TRUE){
    ddc_result = DDC(x)
    indice = ddc_result$indcells
    x[indice] = NA
  }
  replicates = matrix(NA,ncol(x),R)
  for (i in 1:R){
    x1 = x[sample(nrow(x),nrow(x),replace=TRUE),]
    data = kNN(as.data.frame(x1),k=k)
    data = data[,1:(ncol(data)/2)]
    fit = lmrob(data,control=lmrob.control(max.it=1000,k.max=1000,maxit.scale = 1000))
    replicates[,i] = fit$coefficients
  }
  pool_estimate = apply(replicates,1,mean)
  pool_var = apply(replicates,1,var)
  pool_sd = sqrt(pool_var)
  pool_z = pool_estimate/pool_sd
  pool_p = matrix(NA,ncol(x),1)
  for (i in 1:ncol(x)){
    pool_p[i,] = 2*pnorm(abs(pool_z[i]),lower.tail = FALSE)
  }
  su=cbind(pool_estimate,pool_sd,pool_z,pool_p)
  colnames(su) = c('Estimate','Std.Error','z.values','Pr(>|z|)')
  output = list('replicates' = replicates,'summary'=su)
}
