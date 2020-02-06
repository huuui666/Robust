library(ggplot2)
library(mvtnorm)

add_outlier = function(data, contamination_alpha, type=c('good leverage','vertical','bad leverage')){#generate different type of outlier regardless of dimension
  n = nrow(data)
  p = ncol(data)-1
  contamination_size = n* contamination_alpha
  contan_index = sample(n,contamination_size)
  if (type == 'vertical'){
    y_max = max(data[,1])
    y_con = rnorm(contamination_size, mean = 2*y_max, 0.25^2)
    data[contan_index,1] = y_con
    return(data)
  } else if(type=='good leverage'){
    if (p==1){
      x_max=max(data[,-1])
    }else{
      x_max = apply(data[,-1],2,max)}
    x_con = rmvnorm(contamination_size,mean=5*x_max,sigma = diag(0.25^2,p))
    data[contan_index,-1] = x_con
    return(data)
  }else if(type == 'bad leverage'){
    y_max = max(data[,1])
    y_con = rnorm(contamination_size, mean = 2*y_max, 0.25^2)
    if (p==1){
      x_max=max(data[,-1])
    }else{
      x_max = apply(data[,-1],2,max)}
    x_con = rmvnorm(contamination_size,mean=5*x_max,sigma = diag(0.25^2,p))
    contan_data = cbind(y_con,x_con)
    data[contan_index,] = contan_data
    return(data)
  }
}

new_generate_outlier_graph = function(data, data_initial){# first column is the dependent variable
  z=as.data.frame(data)
  z_clean = as.data.frame(data_initial)
  #z=data
  p = ncol(data)
  n = nrow(data)
  clean_fit = lm(data_initial[,1]~data_initial[,-1])
  clean_beta = clean_fit$coefficients
  prediction = predict(clean_fit, as.data.frame(z[,-1])) #with clean fit on contaminated data
  residuals = z[,1] - prediction
  clean_res_stdev = sigma(clean_fit) #get residual standard dev from model on clean data
  standardized_res = residuals/clean_res_stdev
  
  standard_residual = matrix(standardized_res,n,1)
  if (p>2){
    z_ex = z[,-1]
    z_ex_clean = z_clean[,-1]
    # use mean and covariance of clean data, apply on contaminated data
    ma = mahalanobis(z_ex,apply(z_ex_clean,2,mean),cov(z_ex_clean)) #for mulitple dimension
    mahalanobis_distance = matrix(ma,n,1)
  }else{
    ma = cov(z_clean)[2,2]*(z[,2]-mean(z_clean[,2]))^2
    mahalanobis_distance = matrix(ma,n,1)
  }
  zz = cbind(z,standard_residual,mahalanobis_distance)
  library(ggplot2)
  ggplot(data=zz, aes(x=mahalanobis_distance, y=standard_residual)) +
    geom_point() +
    labs(title="Outlier Detection", x="Mahalanobis Distance", y="Standardized Residual")+ geom_hline(yintercept=2.5)+ geom_hline(yintercept=-2.5)+geom_vline(xintercept = qchisq(0.975,p-1))
}


set.seed(1234)
n=100
p=1
ss = matrix(0.2,p,p)#change correlation as you want
diag(ss)=1
x =  rmvnorm(n,rep(0,p),sigma = ss)
beta=matrix(runif(p,0,5),p,1)
y= x%*%beta + matrix(rnorm(n,0,1),n,1)
data = cbind(y,x)  
data_initial = data
con_data = add_outlier(data,0.2,type='bad leverage')
new_generate_outlier_graph(con_data,data_initial)
