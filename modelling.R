library(jsonlite)
alldata <- fromJSON(txt="~/Downloads/resproj-27787-export (14).json") #worked on 3, 5, 9

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

#block ratio analysis
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

#sinusoid ratio analysis
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

#regressions for block
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
model <- lmer(update ~ error * volatility * set + (error|id), data = long)
summary(model)

#filters (from ResProj_0.4, but now changed)

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

pf_bootstrap <- function(y, num_of_particles) {
  # Bootstrap filter
  # create matrices to store the particle values and weights
  n <- length(y)
  logTa <- Ta <- Wa <- s <- k <- estim_state <- matrix(NA, ncol = num_of_particles, nrow = n) 
  
  # draw the particles for the initial state from the prior distribution
  
  logTa[1,] <- rnorm(num_of_particles, mean=log(5), sd=1)
  # set Ta to e^logTa
  Ta[1,] <- exp(logTa[1,])
  
  estim_state[1,] <- rep(0, num_of_particles)
  k[1:2,] <- rep(0, num_of_particles)
  s[1:2,] <- rep(1000, num_of_particles)
  
  Wa[1,] <- 1/num_of_particles
  
  # loop over time
  for(t in 1:(n-1)) {
    
    # sample particles according to the transition distribution
    
    logTa[t+1,] <- rnorm(num_of_particles, mean=logTa[t,], sd = .1)
    # set Ta to e^logTa
    Ta[t+1,] <- exp(logTa[t+1,])
    
    k[t+1,] <- (s[t,] + Ta[t,]) / (s[t,] + Ta[t,] + 4^2)
    s[t+1,] <- (1 - k[t+1,])*(s[t,] + Ta[t,])
    estim_state[t+1,] <- estim_state[t,] + k[t+1,] * (y[t] - estim_state[t,])
    
    # compute the weights
    Wa[t+1,] <- dnorm(y[t+1], mean = estim_state[t,], sd = sqrt(s[t+1,] + Ta[t+1,] + 4^2)) * Wa[t,] # could have multiplied by W, but not necessary with uniform weights
    Wa[t+1,] <- Wa[t+1,]/sum(Wa[t+1,]) # normalize
    
    # draw indices of particles with systematic resampling
    idx <- resample_systematic(Wa[t+1,])
    
    # implicitly W <- 1/num_of_particles
    logTa[t+1,] <- logTa[t+1,idx]
    Ta[t+1,] <- Ta[t+1,idx]
    k[t+1,] <- k[t+1,idx]
    s[t+1,] <- s[t+1,idx]
    estim_state[t+1,] <- estim_state[t+1,idx]
    
    # reset the weights
    Wa[t+1,] <- 1/num_of_particles
  }
  
  return(rowSums(sqrt(Ta)*Wa))
}


model_tracking_volatility <- function(estim_noise, observations, num_of_particles) {
  estim_vol <- pf_bootstrap(observations, num_of_particles)
  eta <- 0
  pred <- 0
  n <- length(observations) - 1
  
  for (t in 1:n) {
    eta[t+1] <- estim_vol[t+1]^2 / (estim_vol[t+1]^2 + estim_noise)
    pred[t+1] <- pred[t] + eta[t+1] * (observations[t] - pred[t])
  }
  return(pred)
}

model_ignoring_volatility <- function(estim_noise, learning_rate, observations) {
  eta <- learning_rate
  pred <- 0
  n <- length(observations) - 1
  for (t in 1:n) {
    pred[t+1] <- pred[t] + eta * (observations[t] - pred[t])
  }
  return(pred)
}

library(stats4)
require(truncnorm)

tvm_errors <- ivm_errors <- lrs <- array()
a <- 1
for (subject in names(locations_block_df)) {
  # subject <- "YfKE878IgpYml15t6RlqJi8Hc9j1"
  observations <- locations_block_df[[subject]]*300
  subject_predictions <- predictions_block_df[[subject]]*300
  tvm_predictions <- model_tracking_volatility(4, observations, 2000)
  tvm_errors[a] <- sum((subject_predictions-tvm_predictions)^2)/100
  
  n <- 1
  learning_rates <- ivm_to_data_error <- ivm_to_subject_error <- array()
  for (i in seq(0, 1, 0.01)) {
    ivm_predictions <- model_ignoring_volatility(4, i, observations)
    learning_rates[n] <- i
    ivm_to_data_error[n] <- sum((observations - ivm_predictions)^2)/100
    ivm_to_subject_error[n] <- sum((subject_predictions - ivm_predictions)^2)/100
    n = n + 1
  }
  
  lrs[a] <- learning_rates[match(min(ivm_to_subject_error), ivm_to_subject_error)]
  ivm_errors[a] <- min(ivm_to_subject_error)
  
  a <- a + 1
}

sum(ivm_errors)/length(ivm_errors)
sum(tvm_errors)/length(tvm_errors)
sum(lrs)/length(lrs)

#ivm fitted to data
learning_rate_to_data <- learning_rates[match(min(ivm_to_data_error), ivm_to_data_error)]
ivm_predictions <- model_ignoring_volatility(4, learning_rate_to_data, observations)
ivm_fitted_to_data_error <- sum((subject_predictions - ivm_predictions)^2)/100
#ivm fitted to subject
learning_rate_to_subject <- learning_rates[match(min(ivm_to_subject_error), ivm_to_subject_error)]
ivm_predictions <- model_ignoring_volatility(4, learning_rate_to_subject, observations)
ivm_fitted_to_subject_error <- sum((subject_predictions - ivm_predictions)^2)/100

#checking
tvm_error
learning_rate_to_data
ivm_fitted_to_data_error
learning_rate_to_subject
ivm_fitted_to_subject_error



est_vol_bootstrap <- pf_bootstrap(y = as.numeric(locations)*300, num_of_particles = 2000)
plot(1:100, est_vol_bootstrap, type = "l", ylab = "est vol bs", xlab = "Day (t)")

pred_tv <- participant_tracking_volatility(estim_noise = 4, estim_vol = est_vol_bootstrap, observations = as.numeric(locations))

mse_ignore_vol_model <- 0
n <- 1
for (i in seq(0, 1, 0.01)) {
  pred_ignore_vol <- participant_ignoring_volatility(estim_noise = 4, learning_rate = i, observations = as.numeric(locations))
  mse_ignore_vol_model[n] <- sum((as.numeric(locations) - pred_ignore_vol)^2)/100
  n = n + 1
}
mse_ignore_vol_model <- min(mse_ignore_vol_model)

pred_ntv <- participant_ignoring_volatility(estim_noise = 4, learning_rate = 0.5, observations = as.numeric(locations))

mse_tv_vs_participant <- sum((as.numeric(predictions) - pred_tv)^2)/100
mse_ntv_vs_participant <- sum((as.numeric(predictions) - pred_ntv)^2)/100

mse_participant <- sum((as.numeric(locations) - as.numeric(predictions))^2)/100
mse_tv_vs_data <- sum((as.numeric(locations) - pred_tv)^2)/100
mse_ntv_vs_data <- sum((as.numeric(locations) - pred_ntv)^2)/100
