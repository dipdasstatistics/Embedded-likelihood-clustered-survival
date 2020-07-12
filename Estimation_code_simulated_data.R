#estimation functions and implementation on simulated data 
#load packages
library(survival)
library(mice)

#functions required to run simulated data
#Transform the data using AFT model
transformed_time <- function(time,x_1, x_2, b1, b2, cluster_effect){
  time * exp(-x_1*b1- x_2* b2 -cluster_effect)
}


#Delta_0 estimation (Nelson-Aalen estimate)
Lambda_0 <- function(b1, b2, x_1, x_2, time, cluster_effect, status){
  dat <- data.frame(new_time = transformed_time(time, x_1, x_2, b1, b2, cluster_effect), status) 
  nelsonaalen(data = dat, timevar = new_time, statusvar = status)
}


#score function for beta1 of AFT-PH model
score_b1 <- function(b1, b2, x_1, x_2, time, cluster_effect, status){
  abs(sum(na.omit(x_1*(status- Lambda_0(b1 = b1, b2 = b2,
  time= time, x_1 = x_1, x_2= x_2, cluster_effect, status = status)))))
}

#score function for beta2 of AFT-PH model
score_b2 <- function(b1, b2, x_1, x_2, time, cluster_effect, status){
  abs(sum(na.omit(x_2*(status- Lambda_0(b1 = b1, b2 = b2,
  time= time, x_1 = x_1, x_2= x_2, cluster_effect, status = status)))))
}

#score function for b1 of AFT-Proportional odds
score_PO_b1 <- function(b1, b2, cluster_effect, x1, x2, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2)
  dat_po <- data.frame(t , x1, x2, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x1[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}


#score function for b2 of AFT-Proportional odds
score_PO_b2 <- function(b1, b2,  cluster_effect, x1, x2, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2)
  dat_po <- data.frame(t , x1, x2, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x2[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}

# Estimating function for AFT-PH model
estimated_ph <- function(x_1, x_2, time, status, cluster){
  n <-as.vector(data.frame(table(cluster))[,2])
  par_ph <- matrix(NA, nrow = 2, ncol = 10)
  estimated_ph_temp <- survreg(Surv(time, status)~ x_1 + x_2 + frailty(cluster, distribution = "gaussian"), dist = "lognormal")
  par_ph[1,1] <- estimated_ph_temp$coefficients[2] 
  par_ph[2,1]<- estimated_ph_temp$coefficients[3]   
 
  for(i in 1:9){
    new_time_ph <- abs(time - x_1*par_ph[1,i] - x_2 * par_ph[2,i])
    temp_ph_surv <- survreg(Surv(new_time_ph, status)~ frailty(cluster, distribution = "gaussian"), dist = "lognormal")
    par_ph[1,i+1] <- optimise(f = score_b1, interval = c(floor(par_ph[1,1]),ceiling(par_ph[1,1])), b2 = par_ph[2,i],  
                    x_1 = x_1, x_2 = x_2,time = time, cluster_effect = rep(temp_ph_surv$frail, 
                    times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    par_ph[2,i+1] <- optimise(f = score_b2, interval = c(floor(par_ph[2,1]),ceiling(par_ph[2,1])), b1 = par_ph[1,i],
                              x_1 = x_1, x_2 = x_2,time = time, cluster_effect = rep(temp_ph_surv$frail, 
                              times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    
    
    if(abs(par_ph[1, i+1]- par_ph[1, i]) < 0.05 && abs(par_ph[2, i+1]- par_ph[2, i]) < 0.05) break
    
  }
  par_ph[,i+1]    
}




# Estimating function for AFT-PO model
estimated_po <- function(x_1,x_2, time, status, cluster){
  n <-as.vector(data.frame(table(cluster))[,2])
  par_po <- matrix(NA, nrow = 2, ncol = 10)
  estimated_po_temp <- survreg(Surv(time, status)~ x_1 + x_2 + frailty(cluster, distribution = "gaussian"), dist = "lognormal")
  par_po[1,1] <- estimated_po_temp$coefficients[2] 
  par_po[2,1]<- estimated_po_temp$coefficients[3]   
  
  for(i in 1:9){
    new_time_po <- abs(time - x_1*par_po[1,i] - x_2 * par_po[2,i])
    temp_po_surv <- survreg(Surv(new_time_po, status)~ frailty(cluster, distribution = "gaussian"), dist = "lognormal")
    par_po[1,i+1] <- optimise(f = score_PO_b1, interval = c(floor(par_po[1,1]),ceiling(par_po[1,1])), b2 = par_po[2,i],
                              x1 = x_1, x2 = x_2, time = time, cluster_effect = rep(temp_po_surv$frail,
                              times = n[which(n != 0)]), status = status, maximum = F)$minimum
    par_po[2,i+1] <- optimise(f = score_PO_b2, interval = c(floor(par_po[2,1]),ceiling(par_po[2,1])), b1 = par_po[1,i],
                              x1 = x_1, x2 = x_2,time = time, cluster_effect = rep(temp_po_surv$frail,
                             times = n[which(n != 0)]), status = status, maximum = F)$minimum
   
    
    if(abs(par_po[1, i+1]- par_po[1, i]) < 0.05 && abs(par_po[2, i+1]- par_po[2, i]) < 0.05) break
  }
  par_po[,i+1]
}

#estimation of beta1, beta2 using AFT-PH, AFT-PO and Survreg function
estimated_ph(x_1 = x1, x_2 = x2, time = time, status = censor_status, cluster = init_clus)
estimated_po(x_1 = x1, x_2 = x2, time = time, status = censor_status, cluster = init_clus)
survreg(Surv(time, censor_status)~ x1+ x2+ frailty(cluster, distribution = "gaussian"), dist = "lognormal")
