rm(list = ls())
#This functions and codes were run to derive the estimate for Lung and Bronchus cancer data

#load packages
library(survival)
library(mice)

#functions
#Transform the data using AFT model
transformed_time <- function(time,x_1, x_2, x_3, x_4, x_5, b1, b2, b3, b4, b5, cluster_effect){
  time * exp(-x_1*b1- x_2* b2 - x_3* b3 - x_4*b4- x_5* b5 -cluster_effect)
}


#Delta_0 estimation (Nelson Aalen estimate)
Lambda_0 <- function(b1, b2, b3, b4, b5, x_1, x_2, x_3, x_4, x_5, time, cluster_effect, status){
  dat <- data.frame(new_time = transformed_time(time, x_1, x_2, x_3, x_4, x_5, b1, b2, b3, b4, b5, cluster_effect), status) 
  nelsonaalen(data = dat, timevar = new_time, statusvar = status)
}


#score function for beta1 of AFT-PH model
score_b1 <- function(b1, b2, b3, b4, b5, x_1, x_2, x_3, x_4, x_5, time, cluster_effect, status){
  abs(sum(na.omit(x_1*(status- Lambda_0(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5,
                                        time= time, x_1 = x_1, x_2= x_2, x_3 = x_3, x_4 = x_4, x_5, cluster_effect, status = status)))))
}

#score function for beta2 of AFT-PH model
score_b2 <- function(b1, b2, b3, b4, b5, x_1, x_2, x_3, x_4, x_5, time, cluster_effect, status){
  abs(sum(na.omit(x_2*(status- Lambda_0(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5,
                                        time= time, x_1 = x_1, x_2= x_2, x_3 = x_3, x_4 = x_4, x_5, cluster_effect, status = status)))))
}

#score function for beta3 of AFT-PH model
score_b3 <- function(b1, b2, b3, b4, b5, x_1, x_2, x_3, x_4, x_5, time, cluster_effect, status){
  abs(sum(na.omit(x_3*(status- Lambda_0(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5,
                                        time= time, x_1 = x_1, x_2= x_2, x_3 = x_3, x_4 = x_4, x_5, cluster_effect, status = status)))))
}
#score function for beta4 of AFT-PH model
score_b4 <- function(b1, b2, b3, b4, b5, x_1, x_2, x_3, x_4, x_5, time, cluster_effect, status){
  abs(sum(na.omit(x_4*(status- Lambda_0(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, 
                                        time= time, x_1 = x_1, x_2= x_2, x_3 = x_3, x_4 = x_4, x_5, cluster_effect, status = status)))))
}


#score function for beta5 of AFT-PH model
score_b5 <- function(b1, b2, b3, b4, b5, x_1, x_2, x_3, x_4, x_5, time, cluster_effect, status){
  abs(sum(na.omit(x_5*(status- Lambda_0(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5,
                                        time= time, x_1 = x_1, x_2= x_2, x_3 = x_3, x_4 = x_4, x_5, cluster_effect, status = status)))))
}


#score function for beta1 of AFT-PO model
score_PO_b1 <- function(b1, b2, b3, b4, b5, cluster_effect, x1, x2, x3, x4, x5, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2, x_3 = x3, x_4 = x4, x_5 =x5)
  dat_po <- data.frame(t , x1, x2, x3, x4, x5, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x1[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}


#score function for beta2 of AFT-PO model
score_PO_b2 <- function(b1, b2, b3, b4, b5, cluster_effect, x1, x2, x3, x4, x5, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2, x_3 = x3, x_4 = x4, x_5 =x5)
  dat_po <- data.frame(t , x1, x2, x3, x4, x5, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x2[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}

#score function for beta3 of AFT-PO model
score_PO_b3 <- function(b1, b2, b3, b4, b5, cluster_effect, x1, x2, x3, x4, x5, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2, x_3 = x3, x_4 = x4, x_5 =x5)
  dat_po <- data.frame(t , x1, x2, x3, x4, x5, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x3[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}

#score function for beta4 of AFT-PO model
score_PO_b4 <- function(b1, b2, b3, b4, b5, cluster_effect, x1, x2, x3, x4, x5, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2, x_3 = x3, x_4 = x4, x_5 =x5)
  dat_po <- data.frame(t , x1, x2, x3, x4, x5, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x4[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}

#score function for beta5 of AFT-PO model
score_PO_b5 <- function(b1, b2, b3, b4, b5, cluster_effect, x1, x2, x3, x4, x5, time, status){
  t <- transformed_time(b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5, cluster_effect = cluster_effect, time = time, x_1 = x1, x_2 = x2, x_3 = x3, x_4 = x4, x_5 =x5)
  dat_po <- data.frame(t , x1, x2, x3, x4, x5, status)
  data_PO <- dat_po[order(dat_po$t),]
  
  F_0 <- function(data_PO){
    km <- survfit(Surv(t, status) ~ 1, data = data_PO)
    survest <- stepfun(km$time, c(1, km$surv))
    1- survest(data_PO$t)
  }
  n <- length(F_0(data_PO))
  abs(sum(data_PO$x5[1:n]* (data_PO$status[1:n] - (rep(1, n)+ data_PO$status[1:n]) * F_0(data_PO))))
}



# Estimating function for AFT-PH model
estimated_ph <- function(x_1, x_2, x_3, x_4, x_5, time, status, cluster){
  n <-as.vector(data.frame(table(cluster))[,2])
  par_ph <- matrix(NA, nrow = 5, ncol = 10)
  estimated_ph_temp <- survreg(Surv(time, status)~ x_1 + x_2 + x_3 + x_4 + x_5 + frailty(cluster, distribution = "gaussian"), dist = "lognormal")
  par_ph[1,1] <- estimated_ph_temp$coefficients[2] 
  par_ph[2,1]<- estimated_ph_temp$coefficients[3]   
  par_ph[3,1]<- estimated_ph_temp$coefficients[4]
  par_ph[4,1]<- estimated_ph_temp$coefficients[5]
  par_ph[5,1]<- estimated_ph_temp$coefficients[6]
  for(i in 1:9){
    new_time_ph <- time - x_1*par_ph[1,i] - x_2 * par_ph[2,i]- x_3 * par_ph[3,i]- x_4*par_ph[4,i] - x_5 * par_ph[5,i]
    temp_ph_surv <- survreg(Surv(new_time_ph, status)~ frailty(cluster, distribution = "gaussian"), dist = "lognormal")
    par_ph[1,i+1] <- optimise(f = score_b1, interval = c(floor(par_ph[1,1]),ceiling(par_ph[1,1])), b2 = par_ph[2,i], b3 = par_ph[3,i], b4 = par_ph[4,i], b5 = par_ph[5,i],  
                              x_1 = x_1, x_2 = x_2, x_3 = x_3, x_4 = x_4, x_5 = x_5,
                              time = time, cluster_effect = rep(temp_ph_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    par_ph[2,i+1] <- optimise(f = score_b2, interval = c(floor(par_ph[2,1]),ceiling(par_ph[2,1])), b1 = par_ph[1,i], b3 = par_ph[3,i], b4 = par_ph[4,i], b5 = par_ph[5,i],
                              x_1 = x_1, x_2 = x_2, x_3 = x_3, x_4 = x_4,  x_5 = x_5,
                              time = time, cluster_effect = rep(temp_ph_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    par_ph[3,i+1] <- optimise(f = score_b3, interval = c(floor(par_ph[3,1]),ceiling(par_ph[3,1])), b1 = par_ph[1,i], b2 = par_ph[2,i], b4 = par_ph[4,i], b5 = par_ph[5,i],
                              x_1 = x_1, x_2 = x_2, x_3 = x_3, x_4 = x_4,  x_5 = x_5,
                              time = time, cluster_effect = rep(temp_ph_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    par_ph[4,i+1] <- optimise(f = score_b4, interval = c(floor(par_ph[4,1]),ceiling(par_ph[4,1])), b1 = par_ph[1,i], b2 = par_ph[2,i], b3 = par_ph[3,i], b5 = par_ph[5,i],
                              x_1 = x_1, x_2 = x_2, x_3 = x_3, x_4 = x_4,  x_5 = x_5, 
                              time = time, cluster_effect = rep(temp_ph_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    par_ph[5,i+1] <- optimise(f = score_b5, interval = c(floor(par_ph[5,1]),ceiling(par_ph[5,1])), b1 = par_ph[1,i], b2 = par_ph[2,i], b3 = par_ph[3,i], b4 = par_ph[4,i],
                              x_1 = x_1, x_2 = x_2, x_3 = x_3, x_4 = x_4,  x_5 = x_5,
                              time = time, cluster_effect = rep(temp_ph_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    
    if(abs(par_ph[1, i+1]- par_ph[1, i]) < 0.05 && abs(par_ph[2, i+1]- par_ph[2, i]) < 0.05 && abs(par_ph[3, i+1]- par_ph[3, i]) < 0.05 && abs(par_ph[4, i+1]- par_ph[4, i]) < 0.05 
       && abs(par_ph[5, i+1]- par_ph[5, i])< 0.05) break
    
  }
  par_ph[,i+1]    
}




# Estimating function for AFT-PO model
estimated_po <- function(x_1,x_2, x_3, x_4, x_5, time, status, cluster){
  n <-as.vector(data.frame(table(cluster))[,2])
  par_po <- matrix(NA, nrow = 5, ncol = 10)
  estimated_po_temp <- survreg(Surv(time, status)~ x_1 + x_2 + x_3+ x_4+ x_5+ frailty(cluster, distribution = "gaussian"), dist = "lognormal")
  par_po[1,1] <- estimated_po_temp$coefficients[2] 
  par_po[2,1]<- estimated_po_temp$coefficients[3]   
  par_po[3,1]<- estimated_po_temp$coefficients[4]  
  par_po[4,1] <- estimated_po_temp$coefficients[5] 
  par_po[5,1]<- estimated_po_temp$coefficients[6]   
  for(i in 1:9){
    new_time_po <- time - x_1*par_po[1,i] - x_2 * par_po[2,i] - x_3 * par_po[3,i]- x_4*par_po[4,i] - x_5 * par_po[5,i]
    temp_po_surv <- survreg(Surv(new_time_po, status)~ frailty(cluster, distribution = "gaussian"), dist = "lognormal")
    par_po[1,i+1] <- optimise(f = score_PO_b1, interval = c(floor(par_po[1,1]),ceiling(par_po[1,1])), b2 = par_po[2,i], b3 = par_po[3,i], b4 = par_po[4,i], b5 = par_po[5,i],
                              x1 = x_1, x2 = x_2, x3 = x_3, x4 = x_4, x5 = x_5,
                              time = time, cluster_effect = rep(temp_po_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    par_po[2,i+1] <- optimise(f = score_PO_b2, interval = c(floor(par_po[2,1]),ceiling(par_po[2,1])), b1 = par_po[1,i], b3 = par_po[3,i], b4 = par_po[4,i], b5 = par_po[5,i],
                              x1 = x_1, x2 = x_2, x3 = x_3, x4 = x_4, x5 = x_5,
                              time = time, cluster_effect = rep(temp_po_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    par_po[3,i+1] <- optimise(f = score_PO_b3, interval = c(floor(par_po[3,1]),ceiling(par_po[3,1])), b1 = par_po[1,i], b2 = par_po[2,i], b4 = par_po[4,i], b5 = par_po[5,i],
                              x1 = x_1, x2 = x_2, x3 = x_3, x4 = x_4, x5 = x_5, 
                              time = time, cluster_effect = rep(temp_po_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    par_po[4,i+1] <- optimise(f = score_PO_b4, interval = c(floor(par_po[4,1]),ceiling(par_po[4,1])), b1 = par_po[1,i], b2 = par_po[2,i], b3 = par_po[3,i], b5 = par_po[5,i],
                              x1 = x_1, x2 = x_2, x3 = x_3, x4 = x_4, x5 = x_5,
                              time = time, cluster_effect = rep(temp_po_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    par_po[5,i+1] <- optimise(f = score_PO_b5, interval = c(floor(par_po[5,1]),ceiling(par_po[5,1])), b1 = par_po[1,i], b2 = par_po[2,i], b3 = par_po[3,i], b4 = par_po[4,i],
                              x1 = x_1, x2 = x_2, x3 = x_3, x4 = x_4, x5 = x_5, 
                              time = time, cluster_effect = rep(temp_po_surv$frail, times = n[which(n != 0)]), status = status, maximum = F)$minimum
    
    if(abs(par_po[1, i+1]- par_po[1, i]) < 0.05 && abs(par_po[2, i+1]- par_po[2, i]) < 0.05 && abs(par_po[3, i+1]- par_po[3, i]) < 0.05   &&
       abs(par_po[4, i+1]- par_po[4, i]) < 0.05 && abs(par_po[5, i+1]- par_po[5, i]) < 0.05) break
  }
  par_po[,i+1]
}

#FOllowing codes are for cleaning real data (Lung and Bronchus cancer data)
setwd("D:\\MS\\Thesis\\SEER\\") 
seer_data <- na.omit(read.table(file = "final_data/RESPIR/CA_state.txt", header = TRUE, sep = ";", colClasses = c("FIPS" = "character", "state_fips" = "character")))

seer_clean_x1 <- subset(subset(seer_data, time_1 >0), race != "unknown")
seer_clean_x3 <- subset(seer_clean_x1, marital_status != "unknown")
seer <- seer_clean_x3[order(seer_clean_x3$FIPS),]
seer_x1 <- ifelse(seer$race == "black", 0, 1) #black=0, white & other=1
seer_x2 <- ifelse(seer$sex == "female", 1, 0)#female=1, male = 0
seer_x3 <- ifelse(seer$marital_status == "single" | seer$marital_status == "unmarried",0, 1)# ever married = 1, never got married=0
seer_x4 <- seer$adult_smoking 
seer_x5 <- seer$excessive_drinking

seer_time <- seer$time_1
seer_status <- seer$status
seer_cluster <- seer$FIPS
#Data is prepared to implement

real_data_survreg <- survreg(Surv(seer_time, seer_status)~ seer_x1 + seer_x2+ seer_x3 + seer_x4 + seer_x5 + frailty(seer_cluster, distribution = "gaussian"), dist = "lognormal", data = seer)
real_data_ph <- estimated_ph(x_1 = seer_x1, x_2 = seer_x2, x_3 = seer_x3, x_4 = seer_x4, x_5 = seer_x5, time = seer_time, status = seer_status, cluster = seer_cluster)
real_data_po <- estimated_po(x_1 = seer_x1, x_2 = seer_x2, x_3 = seer_x3, x_4 = seer_x4, x_5 = seer_x5, time = seer_time, status = seer_status, cluster = seer_cluster)


temp_time_ph_final <- seer_time - real_data_ph[1] * seer_x1 - real_data_ph[2] * seer_x2- real_data_ph[3] * seer_x3- real_data_ph[4] * seer_x4- real_data_ph[5] * seer_x5
temp_time_po_final <- seer_time - real_data_po[1] * seer_x1 - real_data_po[2] * seer_x2- real_data_po[3] * seer_x3- real_data_po[4] * seer_x4- real_data_po[5] * seer_x5



#bootstrap standard errors for AFT-PH, AFT-PO, Survreg function
data_boot <- data.frame(seer_time, seer_status, seer_x1, seer_x2, seer_x3, seer_x4, seer_x5, seer_cluster)
loop_ph_b1 <- NULL
loop_ph_b2 <- NULL
loop_ph_b3 <- NULL
loop_ph_b4 <- NULL
loop_ph_b5 <- NULL
loop_po_b1 <- NULL
loop_po_b2 <- NULL
loop_po_b3 <- NULL
loop_po_b4 <- NULL
loop_po_b5 <- NULL
loop_survreg_b1 <- NULL
loop_survreg_b2 <- NULL
loop_survreg_b3 <- NULL
loop_survreg_b4 <- NULL
loop_survreg_b5 <- NULL
loop_frail_ph <- NULL
loop_frail_po <- NULL
loop_frail_survreg <- NULL

for(i in 1:10){
  dat <- data_boot[sample(x = nrow(data_boot), size = nrow(data_boot), replace = T),]
  data <- dat[order(dat$seer_cluster),]
  loop_x1 <- data$seer_x1
  loop_x2 <- data$seer_x2
  loop_x3 <- data$seer_x3
  loop_x4 <- data$seer_x4
  loop_x5 <- data$seer_x5
  loop_time <- data$seer_time
  loop_status <- data$seer_status 
  loop_cluster <- data$seer_cluster
  loop_ph <- estimated_ph(x_1 = loop_x1, x_2 = loop_x2,x_3 = loop_x3,x_4 = loop_x4, x_5 = loop_x5, time = loop_time, status = loop_status, cluster = loop_cluster)
  loop_po <- estimated_po(x_1 = loop_x1, x_2 = loop_x2, x_3 = loop_x3, x_4 = loop_x4, x_5 = loop_x5,time = loop_time, status = loop_status, cluster = loop_cluster)
  loop_survreg <- survreg(Surv(loop_time, loop_status)~ loop_x1 + loop_x2 + loop_x3 +loop_x4+  loop_x5+ frailty(loop_cluster, distribution = "gaussian"), dist = "lognormal", data = data)
  loop_time_ph <- loop_time - loop_ph[1] * loop_x1 - loop_ph[2] * loop_x2 - loop_ph[3] * loop_x3- loop_ph[4] * loop_x4 - loop_ph[5] * loop_x5
  loop_time_po <- loop_time - loop_po[1] * loop_x1 - loop_po[2] * loop_x2- loop_po[3] * loop_x3- loop_po[4] * loop_x4- loop_po[5] * loop_x5
  
  loop_ph_b1[i] <- loop_ph[1]
  loop_ph_b2[i] <- loop_ph[2]
  loop_ph_b3[i] <- loop_ph[3]
  loop_ph_b4[i] <- loop_ph[4]
  loop_ph_b5[i] <- loop_ph[5]
  loop_po_b1[i] <- loop_po[1]
  loop_po_b2[i] <- loop_po[2]
  loop_po_b3[i] <- loop_po[3]
  loop_po_b4[i] <- loop_po[4]
  loop_po_b5[i] <- loop_po[5]
  loop_survreg_b1[i] <- loop_survreg$coefficients[2]
  loop_survreg_b2[i] <- loop_survreg$coefficients[3]
  loop_survreg_b3[i] <- loop_survreg$coefficients[4]
  loop_survreg_b4[i] <- loop_survreg$coefficients[5]
  loop_survreg_b5[i] <- loop_survreg$coefficients[6]
  
  loop_frail_ph[i] <- survreg(Surv(loop_time_ph, loop_status)~frailty(loop_cluster, distribution = "gaussian"), dist = "lognormal")$history$`frailty(loop_cluster, distribution = "gaussian")`$`theta`
  
  loop_frail_po[i] <- survreg(Surv(loop_time_po, loop_status)~frailty(loop_cluster, distribution = "gaussian"), dist = "lognormal")$history$`frailty(loop_cluster, distribution = "gaussian")`$`theta`
  
  loop_frail_survreg[i] <- loop_survreg$history$`frailty(loop_cluster, distribution = "gaussian")`$`theta`
}


#results for Ph model
real_data_ph[1]# b1 from ph
real_data_ph[2]# b2 from ph
real_data_ph[3]# b3 from ph
real_data_ph[4]# b4 from ph
real_data_ph[5]# b5 from ph

#results for PO model
real_data_po[1]# b1 from po
real_data_po[2]# b2 from po
real_data_po[3]# b3 from po
real_data_po[4]# b4 from po
real_data_po[5]# b5 from po

#se for betas
c(sd(loop_ph_b1), sd(loop_ph_b2), sd(loop_ph_b3), sd(loop_ph_b4), sd(loop_ph_b5))
c(sd(loop_po_b1), sd(loop_po_b2), sd(loop_po_b3), sd(loop_po_b4), sd(loop_po_b5))
c(sd(loop_survreg_b1), sd(loop_survreg_b2), sd(loop_survreg_b3), sd(loop_survreg_b4), sd(loop_survreg_b5))



#write down frailty and scale for ph model
survreg(Surv(temp_time_ph_final, seer_status)~frailty(seer_cluster, distribution = "gaussian"), dist = "lognormal")

#write down frailty and scale for po model
survreg(Surv(temp_time_po_final, seer_status)~frailty(seer_cluster, distribution = "gaussian"), dist = "lognormal")

#results for parameter estimates, frailty, scale from survreg
real_data_survreg




#se for frailty
sd(loop_frail_ph)
sd(loop_frail_po)
sd(loop_frail_survreg)



