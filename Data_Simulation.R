set.seed(9)

#set number of clusters and members in each cluster
g<- 10
n <- 10
censor_rate <- 0.75

#run function.r to run the codes in this file

my_frailty_ph <- NULL
my_scale_ph <- NULL
est.b2.ph <- NULL
est.b1.ph <- NULL

my_frailty_po <- NULL
my_scale_po <- NULL
est.b2.po <- NULL
est.b1.po <- NULL

survreg_b1 <- NULL
survreg_b2 <- NULL
survreg_scale <- NULL
survreg_frailty <- NULL

# set initital parameter values
init_b1 <- log(2)
init_b2 <- log(3)

for(i in 1:100){
#generate cluster effect
clus_effect <- rnorm(n = g, mean = 0, sd = sqrt(1)) 
 
#generate clusters and it's effects
cluster <- rep(1:g, each = n)
init_clus<- rep(clus_effect, each = n) #cluster effects for all cluster

#generate covariate
x1 <- rbinom(n = g * n, size = 1, prob = 0.8)
x2 <- runif(n = g * n, min = 0, max = 1)


#generate random error
e <- rnorm(n = g * n, mean = 0, sd = 0.5)

#wang
b <- c(rep(rnorm(20, mean=0, sd=sqrt(1)),10)) #theta=1
e <- c(rlnorm(200, mean=0, sd=.5)) #scale=0.5
tt <- exp(1+x1-x2*0.5+x3+ b)* e
cc <- c(rexp(200, 1/quantile(tt, .945)))
status <- ifelse(cc <= tt, 0, 1)
t <- pmin(tt,cc)
data2 <- data.frame(t, status, x1,x2,x3,b)
pcc <- 1-sum(status)/200 # 20% censoring
pcc


#generate censoring
censor <- rbinom( n = g * n, size = 1, prob = censor_rate)

# generate time from aft model
time <- exp(x1*init_b1+x2*init_b2 + init_clus +e)

#ordering data by cluster
dat <- data.frame(cluster, time, censor, x1, x2, init_clus, e)
data_PH <- dat[order(dat$cluster),]


#parameters for aft model from survreg function
survreg_result <- survreg(Surv(time, censor) ~ x1 + x2 + frailty(cluster, distribution = "gaussian"), dist = "lognormal")
survreg_b1[i] <- survreg_result$coefficients[2]
survreg_b2[i] <- survreg_result$coefficients[3]
survreg_scale[i] <- survreg_result$scale
survreg_frailty[i] <- survreg_result$`history`$`frailty(cluster, distribution = "gaussian")`$`theta`


est.b1.ph[i] <- optimise(f = score_b1, interval = c(0,2), b2 = init_b2, cluster_effect = rep(survreg_result$frail, each = n), x_1 = x1, x_2 = x2, time = time, status= censor, maximum = F)$minimum
est.b1.po[i] <- optimise(f = score_PO_b1, interval = c(0, 2),b2 = init_b2, cluster_effect = rep(survreg_result$frail, each = n), x1 = x1, x2 =x2, time = time, status = censor, maximum = F)$minimum

est.b2.ph[i] <- optimise(f = score_b2, interval = c(0,2), b1 = init_b1, cluster_effect = rep(survreg_result$frail, each = n), x_1 = x1, x_2 = x2, time = time, status= censor, maximum = F)$minimum
est.b2.po[i] <- optimise(f = score_PO_b2, interval = c(0, 2),b1 = init_b1, cluster_effect = rep(survreg_result$frail, each = n), x1 = x1, x2 =x2, time = time, status = censor, maximum = F)$minimum

#frailty and scale result from ph
time_new_ph <- exp(log(time)- est.b1.ph[i]*x1 - est.b2.ph[i]*x2)
result_ph <- survreg(Surv(time_new_ph, censor)~frailty(cluster, distribution = "gaussian"), dist = "lognormal")
my_scale_ph[i] <- result_ph$scale
my_frailty_ph[i] <- result_ph$`history`$`frailty(cluster, distribution = "gaussian")`$`theta`


#frailty and scale result from po
time_new_po <- exp(log(time)- est.b1.po[i]*x1 - est.b2.po[i]*x2)
result_po <- survreg(Surv(time_new_po, censor)~frailty(cluster, distribution = "gaussian"), dist = "lognormal")
my_scale_po[i] <- result_po$scale
my_frailty_po[i] <- result_po$`history`$`frailty(cluster, distribution = "gaussian")`$`theta`

}


#Results comparison



#bias for beta for ph
round(c(mean(est.b1.ph)-init_b1, mean(est.b2.ph)-init_b2),3)#bias for PH
round(c(sd(est.b1.ph), sd(est.b2.ph)),3)#ese for ph
round(c(mean(est.b1.po)-init_b1, mean(est.b2.po)-init_b2),3)#bias for PO
round(c(sd(est.b1.po), sd(est.b2.po)),3)#ese for po
round(c(mean(survreg_b1)-init_b1, mean(survreg_b2)-init_b2),3) #bias for survreg
round(c(sd(survreg_b1), sd(survreg_b2)),3)#ese for survreg
round(mean(my_frailty_ph)-1,3) #frailty bias for PH
round(mean(my_frailty_po)-1,3)#frailty bias for PO
round(sd(my_frailty_ph),3)#ese for frailty for ph
round(sd(my_frailty_po),3)#ese for frailty for po
round(mean(survreg_frailty)- 1,3) #frailty bias for survreg
round(sd(survreg_frailty),3)
round(mean(survreg_scale),3)#scale for Survreg
round(mean(my_scale_ph),3)#scale for PH
round(mean(my_scale_po),3)#scale for PO
            




