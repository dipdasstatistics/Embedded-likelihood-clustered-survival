#set number of clusters and members in each cluster
g<- 10
n <- 10

# set initial parameter values
init_b1 <- log(2)
init_b2 <- log(3)

#generate cluster effect
clus_effect <- rnorm(n = g, mean = 0, sd = sqrt(1)) 
  
#generate cluster number and it's effects
cluster <- rep(1:g, each = n)
init_clus<- rep(clus_effect, each = n) #cluster effects for all cluster
  
#generate covariates
x1 <- rbinom(n = g * n, size = 1, prob = 0.8)
x2 <- runif(n = g * n, min = 0, max = 1)
  
#generate random error
e <- rnorm(n = g * n, mean = 0, sd = 1)

# generate time from aft model and censoring time
aft_time <- exp(x1*init_b1+x2*init_b2 + init_clus +0.5 * e)
censor_time <- runif(n = 100, min = 0, max = 25)

#survival time with censor_status and censoring rate
time <- pmin(aft_time,censor_time)
censor_status <- ifelse(censor_time <= time, 0, 1)
censor_rate <- 1-sum(censor_status)/(g*n) 
censor_rate

#simulated data frame
sim_dat <- data.frame(time, censor_status, x1, x2, init_clus)
