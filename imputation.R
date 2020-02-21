# -----------------------------------------
# Author:
# *Enter your name and student number here*
# -----------------------------------------

## Use this code skeleton to implement the procedures for obtaining point
## estimates and valid standard errors via multiple imputation or the 
## bootstrap.  Please submit your implementation together with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each student's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.



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
# any other information you want to keep

multimp <- function(xy, m, DDC = FALSE, ...) {
  # You should set a sensible default for the number of imputations m.
  #
  # You can use function DDC() from package cellWise, as well as function 
  # irmi() from package VIM for imputations.
  
  if (DDC == TRUE){
    ddc_result = DDC(xy)
    indice = ddc_result$indall
    xy[indice] = NA
  }
  missing_ob = nrow(x[!complete.cases(xy),])
  missing_percentage = missing_ob / nrow(xy)
  if (missing(m)) m = ceiling(missing_percentage*100)
  imputed= irmi(xy,mi=m)
  output = list('imputed'= imputed, 'm' = m)
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
# any other information you want to keep

fit <- function(xyList, ...) {
  # You can use function lmrob() from package robustbase for the MM-estimator.
  m = xyList$m
  data = xyList$imputed
  model = list()
  for (i in 1:m){
    models[[i]] = lmrob(data[[i]])
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
# p-value in the fifth column (see slide 50 of Lecture 5 for an example)

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
  v_m = (m-1)/gamma_m
  vcomp = nrow(models[[1]]$model)-p-1
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
#              fourth column (see slide 29 of Lecture 5 for an example)

bootstrap <- function(x, R, k, DDC = FALSE, ...) {
  # You should set a sensible default for the number of bootstrap replicates R 
  # and the number of neighbors k.
  #
  # You can use function DDC() from package cellWise, function kNN() from 
  # package VIM for imputations, and function lmrob() from package robustbase 
  # for the MM-estimator.
  
  
}
