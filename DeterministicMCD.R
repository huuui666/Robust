# --------------------------------------------------------------------
# Author:
# *Enter your group number here, as well as names and student numbers*
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the plug-in
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each group's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.



## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
  # *enter your code here*
  z1 = apply(z,2,tanh)
  return (cor(z1))
}

# spearman correlation matrix
corSpearman <- function(z) {
  # *enter your code here*
  return (cor(z, method='spearman'))
  
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
  # *enter your code here*
  n = nrow(z)
  t = qnorm ((apply(z,2,rank) - 1/3) / (n + 1/3) )
  return (cor(t))
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  # *enter your code here*
  z_norm = sqrt(apply(z^2,1,sum))
  k = z / z_norm
  return (crossprod(k))
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  # *enter your code here*
  z_norm = sqrt(apply(z^2,1,sum))
  Order = order(z_norm)
  first_half = Order[1:round(nrow(z) / 2)]
  return (cov(z[first_half,]))
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  # *enter your code here*
  # Hint: have a look at function covOGK() in package robustbase
  ogk = covOGK(z,sigmamu = s_Qn)
  return (ogk$cov)
}



## Main function for deterministic MCD algorithm

# Input:
# x ....... data matrix
# alpha ... proportion of observations to be used for subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return

covDetMCD <- function(x, alpha, ...) {
  # *enter your code here*
  #
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function h.alpha.n() from package robustbase to compute the 
  # subset size.
  n = nrow(x)
  p = ncol(x)
  alpha=0.8
  h = h.alpha.n(alpha, n, p)
  cor_choice <- function(x,choice){
    if (choice == 1){
      return (corHT(x))
    }
    else if (choice ==2 ){
      return (corSpearman(x))
    }
    else if (choice == 3){
      return (corNSR(x))
    }
    else if (choice == 4){
      return (covMSS(x))
    }
    else if (choice == 5){
      return (covBACON1(x))
    }
    else {
      return (rawCovOGK(x))
    }
  }
  raw_center = matrix(0,6,p)
  raw_cov = list()
  for (i in 1:6){
    subset = x[sample(n,h),]
    mu = apply(subset,2,mean)
    s = cor_choice(subset,i)
    initial_det = det(s)
    new_det = initial_det - 10^-6
    iter = 0
    while (abs(new_det - initial_det)> 10^-8){
      initial_det = new_det
      distance = mahalanobis(x,mu,s)
      Order = order(distance)
      new_set = Order[1:h]
      new_subset = x[new_set,]
      new_mu = apply(new_subset,2,mean)
      new_s = cor_choice(new_subset,i)
      new_det = det(new_s)
      iter =+ 1
      
    }
    raw_center[i,]=new_mu
    raw_cov[[i]]= new_s
    }
    
    output=list('raw_center'=raw_center, 'raw_cov' = raw_cov,'iter'=iter)
    output
  }
  
  
  
  
  
  




## Function for regression based on the deterministic MCD

# Input:
# x ........ matrix of explanatory variables
# y ........ response variable
# alpha .... proportion of observations to be used for the subset size in the 
#            MCD estimator
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, alpha, ...) {
  # *enter your code here*
}
