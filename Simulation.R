# -----------------------------------------
# Author:
# *Enter your name and student number here*
# zhaohui wang 466503
# -----------------------------------------

## Use this empty file to implement your simulation study.  Please submit your 
## code together with your report.


# Load R script with your implementation of the methods
source("Imputation.R")


# Put the code for your simulations below here



logit = function(x) {exp(x)/(1+exp(x))}

generate_cellwise_outlier = function (data, prob){ # add cellwise outlier into data
  p = ncol(data)
  data_cellwise = matrix(NA, nrow(data),p)
  for (i in 1:nrow(data)){
    I_p = diag(1,p)
    B = diag(rbinom(p,1,prob=prob),p)
    ss1 = matrix(0.2,p,p)#change correlation as you want
    diag(ss1)=1
    Z= rmvnorm(1,runif(p,15,20),sigma = ss1)
    data_cellwise[i,] = data[i,]%*%(I_p - B) + Z %*% B
  }
  return (data_cellwise)
}



add_missing_value = function(data,pattern_frequency, type = c('mar','mcar')){#add two types of missing mechanisms
  p = ncol(data)
  n = nrow(data)

  if (type == 'mar'){
    weight = matrix(runif(p,-5,5),p,p)
    diag(weight) = 0
    index = sample(seq(1,length(pattern_frequency)),size = n , replace=TRUE, prob = pattern_frequency)
    pattern_data = list()
    for (i in 1:length(pattern_frequency)){
      data_per_pattern = data[index==i,]
      score = data_per_pattern %*% weight[i,]
      right_prob = logit(-mean(score)+score)
      set_miss = 1- rbinom(nrow(data_per_pattern),1,right_prob)
      for (j in 1:length(set_miss)){
        if (set_miss[j] == 0){
          data_per_pattern[j,i] = NA
        }
      }
      pattern_data[[i]] = data_per_pattern
    }
    
    return (do.call(rbind,pattern_data))
  }else if (type=='mcar'){
    # all score will be 0.5, since random
    set_miss = 1-rbinom(n,1,0.5)
    for (j in 1:n){
      if(set_miss[j]==0){
        data [j, sample(1:p,1)] = NA
      }
    }
    return (data)
  }
}



mi_imputation = function(data, DDC=c('TRUE','FALSE')){ #make MI all functions into one
  if (DDC==TRUE){
    imput = multimp(data,DDC=TRUE)
  }else if (DDC==FALSE){
    imput = multimp(data,DDC=FALSE)
  }
  myfit = fit(imput)
  return(pool(myfit))
}


library(progress)
#data generation defined parameter
set.seed(1234)
simulation_number = 100
n = 300
p = 9
beta = matrix(runif(p,0,5),p,1)


run_simulation = function(simulation_number,beta,n,p,eplison=0.1, cellwise=TRUE){#simulation function 
  
  ss = matrix(0.7,p,p)#change correlation as you want
  diag(ss)=1
  pattern_frequency= rep(0.25,4)
  MI = list()
  MI_DDC = list()
  BOOT =list()
  BOOT_DDC = list()
  
  MI_mar = list()
  MI_DDC_mar = list()
  BOOT_mar =list()
  BOOT_DDC_mar = list()
  pb <- progress_bar$new(total = simulation_number)
  
  for (i in 1: simulation_number){
    pb$tick()
    Sys.sleep(1 / simulation_number)
    x =  rmvnorm(n,rep(0,p),sigma = ss)
    y= x%*%beta + matrix(rnorm(n,0,1),n,1)
    data = cbind(y,x)  
    if (cellwise==TRUE){
      data =generate_cellwise_outlier(data,eplison)
    }
    data_with_mcar = add_missing_value(data,pattern_frequency,type='mcar')
    data_with_mar = add_missing_value(data,pattern_frequency,type='mar')
    MI[[i]] = mi_imputation(data_with_mcar,DDC = FALSE)
    MI_DDC[[i]] = mi_imputation(data_with_mcar,DDC = TRUE)
    BOOT[[i]] = bootstrap(data_with_mcar,DDC=FALSE)$summary
    BOOT_DDC[[i]] = bootstrap(data_with_mcar,DDC=TRUE)$summary
    MI_mar[[i]] = mi_imputation(data_with_mar,DDC = FALSE)
    MI_DDC_mar[[i]] = mi_imputation(data_with_mar,DDC = TRUE)
    BOOT_mar[[i]] = bootstrap(data_with_mar,DDC=FALSE)$summary
    BOOT_DDC_mar[[i]] = bootstrap(data_with_mar,DDC=TRUE)$summary
    
  }
  output= list('MI_mcar' = MI, 'MI_DDC_mcar' = MI_DDC,'BOOT_mcar' = BOOT, 'BOOT_DDC_mcar' = BOOT_DDC,
          'MI_mar' = MI_mar,'MI_DDC_mar' = MI_DDC_mar, 'BOOT_mar' = BOOT_mar,'BOOT_DDC_mar' = BOOT_DDC_mar)
  
}

simulation_result_with_cellwise = run_simulation(simulation_number,beta,n,p)
simulation_result_without_cellwise = run_simulation(simulation_number,beta,n,p,cellwise = FALSE)


check_interval = function(coeff,beta){ # check if true coefficient fall into 95% confidence interval
  simulation_number = ncol(coeff)
  fall_in = matrix(NA, length(beta),simulation_number)
  std_error = apply(coeff,1,sd)
  for (i in 1:simulation_number){
    #std_error = data_list[[i]][-1,2]
    interval = 1.96* std_error
    lower = coeff[,i] - interval
    upper =coeff[,i] + interval
    fall_in[,i] = ifelse(sapply(beta,function(p) any(lower<p & upper>p)), 1,0)
  }
  return(apply(fall_in,1,mean))
}

check_point = function(coeff,beta){ # output MSE for each variable
  sq_matrix = matrix(NA,nrow(coeff),ncol(coeff))
  for (i in 1: ncol(coeff)){
    sq_matrix[,i]  = (coeff[,i]-beta)^2
  }
  return (apply(sq_matrix,1,mean))
}

check_performance = function(simulaion_result,simulation_number,p,beta){ # check perfromance from the result of simulation
  MI_coeff_mcar = matrix(NA,p,simulation_number)
  MI_DDC_coeff_mcar = matrix(NA,p,simulation_number)
  BOOT_coeff_mcar = matrix(NA,p,simulation_number)
  BOOT_DDC_coeff_mcar = matrix(NA,p,simulation_number)
  for (i in 1:simulation_number){
    MI_coeff_mcar[,i] = simulaion_result$MI_mcar[[i]][-1,1]
    MI_DDC_coeff_mcar[,i] = simulaion_result$MI_DDC_mcar[[i]][-1,1]
    BOOT_coeff_mcar[,i] =simulaion_result$BOOT_mcar[[i]][-1,1]
    BOOT_DDC_coeff_mcar[,i] = simulaion_result$BOOT_DDC_mcar[[i]][-1,1]
    
  }

  
  mcar_MSE = rbind(check_point(MI_coeff_mcar,beta),check_point(MI_DDC_coeff_mcar,beta), check_point(BOOT_coeff_mcar,beta),
                   check_point(BOOT_DDC_coeff_mcar,beta))
  
  mcar_fall_rate = rbind(check_interval(MI_coeff_mcar,beta),check_interval(MI_DDC_coeff_mcar,beta), check_interval(BOOT_coeff_mcar,beta),
                         check_interval(BOOT_DDC_coeff_mcar,beta))
  
 
   MI_coeff_mar = matrix(NA,p,simulation_number)
  MI_DDC_coeff_mar = matrix(NA,p,simulation_number)
  BOOT_coeff_mar = matrix(NA,p,simulation_number)
  BOOT_DDC_coeff_mar = matrix(NA,p,simulation_number)
  
  for (i in 1:simulation_number){
    MI_coeff_mar[,i] = simulaion_result$MI_mar[[i]][-1,1]
    MI_DDC_coeff_mar[,i] = simulaion_result$MI_DDC_mar[[i]][-1,1]
    BOOT_coeff_mar[,i] =simulaion_result$BOOT_mar[[i]][-1,1]
    BOOT_DDC_coeff_mar[,i] = simulaion_result$BOOT_DDC_mar[[i]][-1,1]
    
  }
  mar_MSE = rbind(check_point(MI_coeff_mar,beta),check_point(MI_DDC_coeff_mar,beta), check_point(BOOT_coeff_mar,beta),
                   check_point(BOOT_DDC_coeff_mar,beta))
  
  mar_fall_rate = rbind(check_interval(MI_coeff_mar,beta),check_interval(MI_DDC_coeff_mar,beta), check_interval(BOOT_coeff_mar,beta),
                         check_interval(BOOT_DDC_coeff_mar,beta))
  
  
  output = list(  'mcar_MSE' = mcar_MSE,   'mcar_fall_rate' = mcar_fall_rate, 'mar_MSE' = mar_MSE, 'mar_fall_rate' = mar_fall_rate)
}



performance_result1 = check_performance(simulation_result_with_cellwise,simulation_number,p ,beta)
performance_result2 = check_performance(simulation_result_without_cellwise,simulation_number,p ,beta)


