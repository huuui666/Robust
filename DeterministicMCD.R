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
  k = as.matrix(k)
  k[is.nan(k)]=0
  return (crossprod(k)/nrow(z))
}
#   cross = apply(k,1,function(y) (y%*%t(y)))
#   cross[is.nan(cross)]=0
#   output = matrix(apply(cross,1,sum),ncol(z),ncol(z))
#   return (output/nrow(z))
# }

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
  # D=diag(apply(z,2,qn))
  # Z= apply(z,1,function(y)(solve(D)%*%y))
  # Z=t(Z)
  # U = matrix(NA,ncol(Z),ncol(Z))
  # for (j in 1:ncol(Z)){
  #   for (k in 1:ncol(Z)){
  #     U[j,k] = 0.25 * (qn((Z[,j]+Z[,k])^2)-qn((Z[,j]-Z[,k])^2))
  #   }
  # }
  # E = eigen(U)$vectors
  # V = Z %*% E
  # l = apply(V,2,qn)
  # L = diag(l^2)
  # m = matrix(apply(V,2,median),ncol(z),1)
  # mu = E %*% m
  # sigma = E %*% L %*% t(E)
  # raw_mu = D %*% mu
  # raw_sigma = D %*% sigma %*% t(D)
  # return (raw_sigma)
  ogk = covOGK(z, n.iter=1, sigmamu = s_Qn,rcov=covGK)
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

mycovDetMCD <- function(x, alpha, ...) {
  # *enter your code here*
  #
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function h.alpha.n() from package robustbase to compute the 
  # subset size.
  mysqrtm = function(m){
    E = eigen(m)
    V = E$values
    Q = E$vectors
    Y = Q%*% diag(1/sqrt(V))%*%t(Q)
    return(solve(Y))
  }
  Qnman <- function(var){
    n=length(var)   
    h_k <- floor((n/2) + 1)
    k <- (h_k*(h_k-1))/2
    
    d <- 2.21914
    start <- 1
    values <- matrix(Inf, n*n, 1)
    for(i in 1:n){
      for(j in 1:n){
        if(i < j){
          values[start] <- abs(var[i] - var[j])
          start = start + 1
        }
      }
    }
    sorted <- sort(values, decreasing = FALSE)
    stat <- d*sorted[k]
    return(stat)
  }
  x_initial = x
  x = apply(x,2,function(y) (y-median(y))/Qnman(y))
  
     
  n = nrow(x)
  p = ncol(x)
  h = h.alpha.n(alpha, n, p)
  c_alpha = alpha/ pgamma(qchisq(alpha,p)/2,1+p/2,1)
  c_alpha_reweight = (1-0.025)/(pgamma(qchisq(1-0.025,p)/2,1+p/2,1))
  correct_eigen <- function(x,s){
    eigen_list = eigen(s)
    E = as.matrix(eigen_list$vectors)
    V = x %*% E
    l = apply(V,2,Qnman)
    L = diag(l^2)
    sigma = E %*% L %*% t(E)
    mu = mysqrtm(sigma) %*% apply(x %*% solve(mysqrtm(sigma)),2,median)
    output = list('center' = mu,'cov' = sigma)
    output
  }
  
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
  best_set = list()
  six_det = rep(0,6)
  h0 = ceiling(n/2)
  iter_list = rep(0,6)
  for (i in 1:6){
    #subset = x[sample(n,h),]
    #mu = apply(subset,2,mean)
    initial_s = cor_choice(x,i)
    correction = correct_eigen(x,initial_s)
    mu = correction$center
    s = correction$cov
   
    initial_distance = mahalanobis(x,mu,s)
    Order = order(initial_distance)
    initial_set = Order[1:h0]
    initial_subset = x_initial[initial_set,]
    s=cov(initial_subset)
    mu = apply(initial_subset,2,mean)
    initial_det = det(s)
    #new_det = initial_det + 10^-6
    new_det = Inf
    iter = 0
    while (new_det !=initial_det){
      if (iter>0){
        initial_det = new_det
      }
      
      distance = mahalanobis(x_initial,mu,s)
      
      
      Order = order(distance)
      new_set = Order[1:h]
      new_subset = x_initial[new_set,]
      mu = apply(new_subset,2,mean)
      s = cov(new_subset)
      new_det = det(s)

      iter = iter+1
      
    }
    iter_list[i]=iter
    raw_center[i,]=mu
    raw_cov[[i]]= s
    best_set[[i]] = new_subset
    six_det[i]=det(s)
    }
    min_index = which(six_det == min(six_det))
    if (length(min_index)>1){
      min_index = min_index[1]
    }
    distance_of_best = mahalanobis(x_initial,raw_center[min_index,],raw_cov[[min_index]]*c_alpha)
    best = which(distance_of_best <= qchisq(0.975,p))
    last_best = x_initial[best,]
    re_center = apply(last_best,2,mean)
    re_cov = cov(last_best)
    weights = rep(0,n)
    weights = replace(weights, best, 1)
    
      
    output=list('raw.center'=raw_center[min_index,], 'raw.cov' = raw_cov[[min_index]]*c_alpha,'iter'=iter_list,'center'=re_center, 'cov'=re_cov*c_alpha_reweight,
                'best'=best,'weights'=weights,'det' = six_det)
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
  z = as.matrix(cbind(y,x))
  x = as.matrix(x)
  fit = mycovDetMCD(z,alpha)
  mu = fit$center
  sigma = fit$cov
  beta = solve(sigma[2:nrow(sigma),2:ncol(sigma)]) %*% sigma[2:nrow(sigma),1]
  intercept =   mu[1]- t(mu[2:length(mu)]) %*% beta
  fitted.values = rep(intercept,nrow(x)) + x%*% beta
  residuals = fitted.values - y
  output = list('coefficients'=rbind(intercept,beta), 'MCD' = fit, 'fitted.values' = fitted.values,'residuals' = residuals,'weight'=fit$weights)  
}
