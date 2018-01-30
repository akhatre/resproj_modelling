#parsing data ----

library(jsonlite)
alldata <- fromJSON(txt="~/Downloads/resproj-27787-export (20).json") #worked on 3, 5, 9

ratios_sinusoid <- ratios_block <- data.frame(asd = rep(0,99))
predictions_block_df <- predictions_sinusoid_df <- locations_block_df <- locations_sinusoid_df <- data.frame(asd = rep(0,100))
users <- names(alldata$predictions)
for (user in users) {
  set <- alldata$users[[user]]$set$value
  locations <- alldata$locations[[set]]
  predictions <- alldata$predictions[[user]]$value
  ratio <- data.frame(set = rep(0,99))
  upd <- err <- array()
  for (i in 1:99) {
    upd[i] <- as.numeric(predictions[i+1]) - as.numeric(predictions[i])
    err[i] <- as.numeric(locations[i]) - as.numeric(predictions[i])
    if (err[i] != 0) {ratio[i, 1] <- upd[i]/err[i]} else {ratio[i, 1] <- 0}
  }
  
  if (as.numeric(substr(set,4,nchar(set)))%%2 == 0) {
    ratios_sinusoid[[user]] <- ratio[,1]
    predictions_sinusoid_df[[user]] <- as.numeric(predictions)
    locations_sinusoid_df[[user]] <- as.numeric(locations)
  } else {
    ratios_block[[user]] <- ratio[,1]
    predictions_block_df[[user]] <- as.numeric(predictions)
    locations_block_df[[user]] <- as.numeric(locations)
  }
}

#simple ratio analysis ----

#block
block_low <- block_high <- array()
i <- 1
for (user in names(ratios_block)) {
  if (user == "asd") {
    next
  }
  ratios <- ratios_block[[user]]
  # plot(1:99, ratios, type = "l", ylab = "Ratio", xlab = "Day", ylim = c(0,5), main = paste("B_R:", user))
  # abline(v=c(25,50,75), col="purple")
  # plot(1:100, predictions_block_df[[user]], col="blue", type = "l", ylab = "Preds/loccations", xlab = "Day", main = paste("B_PL:", user))
  # lines(1:100, locations_block_df[[user]], col="red")
  # abline(v=c(25,50,75), col="purple")
  block_low[i] <- (sum(ratios[1:25]) + sum(ratios[51:75]))/50
  block_high[i] <- (sum(ratios[26:50]) + sum(ratios[76:99]))/49
  i <- i + 1
}

sum(block_low)/length(block_low)
sum(block_high)/length(block_high)
plot(1:length(block_low), block_low, col = "blue", type = "l")
lines(1:length(block_low), block_high, col = "red")

#sinusoid

sin_low <- sin_high <- array()
i <- 1
for (user in names(ratios_sinusoid)) {
  if (user == "asd") {
    next
  }
  ratios <- ratios_sinusoid[[user]]
  # plot(1:99, ratios, type = "l", ylab = "Ratio", xlab = "Day", ylim = c(0,5), main = paste("S_R:", user))
  # abline(v=c(25,50,75), col="purple")
  # plot(1:100, predictions_sinusoid_df[[user]], col="blue", type = "l", ylab = "Preds/loccations", xlab = "Day", main = paste("S_PL:", user))
  # lines(1:100, locations_sinusoid_df[[user]], col="red")
  # abline(v=c(25,50,75), col="purple")
  sin_low[i] <- (sum(ratios[1:25]) + sum(ratios[51:75]))/50
  sin_high[i] <- (sum(ratios[26:50]) + sum(ratios[76:99]))/49
  i <- i + 1
}

sum(sin_low)/length(sin_low)
sum(sin_high)/length(sin_high)
plot(1:length(sin_low), sin_low, col = "blue", type = "l", ylim = c(0,5))
lines(1:length(sin_low), sin_high, col = "red")

#regression analysis (mixed effects model) ----
volatility_block <- c(rep(3,24),rep(30,25),rep(3,25),rep(30,26))
volatility_sinusoid <- sin(seq(-0.5*pi, 2.5*pi, length.out = 100))*(27)/2 + (((27)/2) + 3)

long <- data.frame()
for (name in names(predictions_block_df)) {
  if (name == "asd") {
    next
  }
  long <- rbind(long, data.frame(id = name, prediction = predictions_block_df[,name], update = c(diff(predictions_block_df[,name]),0), locations = locations_block_df[,name], volatility = volatility_block, set = "block"))
}

for (name in names(predictions_sinusoid_df)) {
  if (name == "asd") {
    next
  }
  long <- rbind(long, data.frame(id = name, prediction = predictions_sinusoid_df[,name], update = c(diff(predictions_sinusoid_df[,name]),0), locations = locations_sinusoid_df[,name], volatility = volatility_sinusoid, set = "sinusoid"))
}

long$error <- long$locations - long$prediction

long$set <- factor(long$set)
contrasts(long$set) <- contr.helmert(2)

library(lmerTest)
model <- lmer(update ~ error * volatility * set + (error|id), data = long, REML = FALSE)
summary(model)

model1 <- lmer(update ~ error * volatility * set + (error*volatility|id), data = long)
summary(model1)

model2 <- lmer(update ~ error + error : (volatility + set + volatility:set) + (error + error:volatility|id), REML = FALSE, data = long)
summary(model2)

#estimating learning rate from the model
long$set_as_number[long$set == "block"] <- -1
long$set_as_number[long$set == "sinusoid"] <- 1
long$eta <- 0.735  + 0.006 * long$volatility - 0.081 * long$set_as_number + 0.004 * long$set_as_number * long$volatility

ids <- unique(long$id)
etas <- data.frame(long$eta[long$id == ids[1]])
plot(1:100, etas[,1], type = "l", ylim = c(0.5,1))
for (n in 2:length(ids)) {
  etas <-  data.frame(etas, long$eta[long$id == ids[n]])
  lines(1:100, long$eta[long$id == ids[n]])
}
etas <- rowSums(etas)/length(names(etas))
plot(1:100, etas, type = "l")

#fitting bootstrap filter ----
library(truncnorm)
resample_systematic <- function(weights) {
  # input: weights is a vector of length N with (unnormalized) importance weights
  # output: a vector of length N with indices of the replicated particles
  N <- length(weights)
  weights <- weights/sum(weights)# normalize weights
  csum <- cumsum(weights)
  u1 <- runif(1,min=0,max=1/N) # draw a single uniform number
  u <- c(0, seq(1/N, (N-1)/N, length = N-1)) + u1
  idx <- vector("integer",length=length(weights))
  j <- 1
  for(i in 1:N) {
    while (u[i] > csum[j]) {
      j <- j + 1
    }
    idx[i] <- j
  }
  return(idx)
}
ll_pf_bootstrap <- function(par, ta0, y, r, num_of_particles = 1000) { #that also should take y, r, and num_of_pacticles as input variables. Though, num_of_particles is not to be optimsed at the moment.
  # set parametrers
  
  sd_y = exp(par[1]) #noise of observations, real one was 4
  sd_r = 1 #exp(par[2]) #noise of participants' responses, unknown
  if (ta0) {
    sd_ta = 0
  } else {
    sd_ta = exp(par[2]) #sd of particle update distribution
  }
  # Bootstrap filter
  # create matrices to store the particle values and weights
  n <- length(y)
  Ta <- Wa <- s <- k <- estim_state <- matrix(NA, ncol = num_of_particles, nrow = n) 
  
  # draw the particles for the initial state from the prior distribution
  
  Ta[1,] <- runif(num_of_particles, 0, 30^2)
  estim_state[1,] <- rep(0, num_of_particles)
  k[1:2,] <- rep(0, num_of_particles)
  s[1:2,] <- rep(1000, num_of_particles)
  
  Wa[1,] <- 1/num_of_particles
  
  # loop over time
  for(t in 1:(n-1)) {
    
    # sample particles according to the transition distribution
    
    Ta[t+1,] <- rtruncnorm(num_of_particles, a = 0, b = Inf, mean=Ta[t,], sd=sd_ta)
    
    k[t+1,] <- (s[t,] + Ta[t,]) / (s[t,] + Ta[t,] + 4^2)
    s[t+1,] <- (1 - k[t+1,])*(s[t,] + Ta[t,])
    estim_state[t+1,] <- estim_state[t,] + k[t+1,] * (y[t] - estim_state[t,])
    
    # compute the weights
    Wa[t+1,] <- dnorm(y[t+1], mean = estim_state[t,], sd = sqrt(s[t+1,] + Ta[t+1,] + sd_y^2)) * Wa[t,] # could have multiplied by W, but not necessary with uniform weights
    Wa[t+1,] <- Wa[t+1,]/sum(Wa[t+1,]) # normalize
    
    # draw indices of particles with systematic resampling
    idx <- resample_systematic(Wa[t+1,])
    
    # implicitly W <- 1/num_of_particles
    Ta[t+1,] <- Ta[t+1,idx]
    k[t+1,] <- k[t+1,idx]
    s[t+1,] <- s[t+1,idx]
    estim_state[t+1,] <- estim_state[t+1,idx]
    
    # reset the weights
    Wa[t+1,] <- 1/num_of_particles
  }
  
  ## implement Kalman filter to give mu (predicted means)
  estim_vol <- rowSums(sqrt(Ta)*Wa)
  eta <- 0
  pred <- 0
  n <- length(y) - 1
  
  for (t in 1:n) {
    eta[t+1] <- estim_vol[t+1]^2 / (estim_vol[t+1]^2 + sd_y)
    pred[t+1] <- pred[t] + eta[t+1] * (y[t] - pred[t])
  }
  
  mu <- pred
  # return(pred)
  
  mse <- sum((pred-y)^2)/100
  logLik <- -2*sum(dnorm(r, mean=mu, sd=sd_r, log=TRUE))	
  return(mse)
}

#serial optimisation
estimated_par <- list()
for(i in unique(long$id)) {
  data <- subset(long, id==i)
  estimated_par[[i]] <- optim(log(c(1,4,4)), fn = ll_pf_bootstrap, y = data$locations, r = data$prediction)
}

pars_btsrp_sd_ta_0 <- data.frame()
for (i in names(estimated_par)) {
  pars_btsrp_sd_ta_0 <- rbind(pars_btsrp_sd_ta_0, data.frame(noise = estimated_par[[i]]$par[1], sd_r = estimated_par[[i]]$par[2], logLik = estimated_par[[i]]$value))
}
summary(pars_btsrp_sd_ta_0$logLik)
sd(pars_btsrp_sd_ta_0$logLik)

estimated_par <- estimated_par_bootstrap_model
pars_btsrp <- data.frame()
for (i in names(estimated_par)) {
  pars_btsrp <- rbind(pars_btsrp, data.frame(sd_ta = estimated_par[[i]]$par[1], noise = estimated_par[[i]]$par[2], sd_r = estimated_par[[i]]$par[3], logLik = estimated_par[[i]]$value))
}
summary(pars_btsrp$logLik[1:75])
sd(pars_btsrp$logLik[1:75])

#parallel optimisation
library(parallel)
prl_optimisation_btsrt_var_ta <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(log(c(4,1)), method = "SANN", fn = ll_pf_bootstrap, ta0 = FALSE, y = data$locations, r = data$prediction))
  return(est_par)
}
prl_optimisation_btsrt_fixed_ta <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(log(c(4,4)), fn = ll_pf_bootstrap, ta0 = TRUE, y = data$locations, r = data$prediction))
  return(est_par)
}

no_cores <- detectCores() - 1
cluster <- makeCluster(no_cores, type = "FORK")
params_btsrt_var_ta_sann <- parLapply(cluster, unique(long$id), fun = prl_optimisation_btsrt_var_ta)
params_btsrt_fixed_ta <- parLapply(cluster, unique(long$id), fun = prl_optimisation_btsrt_fixed_ta)
stopCluster(cluster)

logLik_btsrt_fixed_ta <- logLik_btsrt_var_ta <- vector()
for (i in 1:length(params_btsrt_fixed_ta)) {
  logLik_btsrt_fixed_ta[i] <- as.numeric(params_btsrt_fixed_ta[[i]]$value)
  logLik_btsrt_var_ta[i] <- as.numeric(params_btsrt_var_ta[[i]]$value)
}
summary(2*2 + logLik_btsrt_fixed_ta)
sd(2*2 + logLik_btsrt_fixed_ta)

summary(2*3 + logLik_btsrt_var_ta)
sd(2*3 + logLik_btsrt_var_ta)

#filters (from ResProj_0.4, but now changed) ----
