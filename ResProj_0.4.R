require(truncnorm)
set.seed(3)

# Generating the data ----

generate_data_volatiliy_linear_increase <- function(num_of_trials, noise, increase, start_state) {
  data <- data.frame(Volatility = rep(4, num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  for (t in 2:num_of_trials) {
    data$Volatility[t] <- data$Volatility[t-1] + increase
    data$State[t] <- data$State[t-1] + rnorm(1, mean = 0, sd = data$Volatility[t])
    data$Observation[t] <- data$State[t] + rnorm(1, mean = 0, sd = noise)
  }
  return(data)
}

generate_data_volatiliy_sinusoid <- function(num_of_trials, noise, max, min, start_state) {
  data <- data.frame(Volatility = rep(4, num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  x <- seq(-0.5*pi, 2.5*pi, length.out = num_of_trials)
  data$Volatility <- sin(x)*(max-min)/2 + (((max-min)/2) + min)
  for (t in 2:num_of_trials) {
    data$State[t] <- data$State[t-1] + rnorm(1, mean = 0, sd = data$Volatility[t])
    data$Observation[t] <- data$State[t] + rnorm(1, mean = 0, sd = noise)
  }
  return(data)
}

generate_data_volatiliy_block <- function(num_of_trials, noise, block1_vol, block2_vol, start_state) {
  data <- data.frame(Volatility = rep(block1_vol, num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  for (t in 2:num_of_trials) {
    
    if (t < num_of_trials/4) {
      data$Volatility[t] <- block1_vol
    } else if (t < num_of_trials/2) {
      data$Volatility[t] <- block2_vol
    } else if (t < num_of_trials/4*3) {
      data$Volatility[t] <- block1_vol
    } else {
      data$Volatility[t] <- block2_vol
    }
    
    sd = data$Volatility[t]
    data$State[t] <- data$State[t-1] + rnorm(1, mean = 0, sd = sd)
    data$Observation[t] <- data$State[t] + rnorm(1, mean = 0, sd = noise)
  }
  return(data)
}

# data <- generate_data_volatiliy_linear_increase(num_of_trials = 100, noise = 4, increase = 0.3, start_state = 0)
# data <- generate_data_volatiliy_sinusoid(num_of_trials = 100, noise = 4, max = 30, min = 3, start_state = 0)
data <- generate_data_volatiliy_block(num_of_trials = 100, noise = 4, block1_vol = 3, block2_vol = 30, start_state = 0)

plot(1:100, data$Observation+100, type = "l", ylab = "Share price", xlab = "Day (t)")

# Plotting generated data ----

par(mfrow=c(1,1))
plot(1:length(data$Observation), data$Observation, type = "l", col = "red")
lines(data$Volatility, col = "green")
legend(c(0, max(data$Observation)), lty = c(1,1), col = c("red", "green"), legend = c("Observations", "Volatility (true)"), bty = "n")

# Filters ----

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
    
    logTa[t+1,] <- rnorm(num_of_particles,mean=logTa[t,], sd = .1)
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

pf_straight <- function(y, num_of_particles) {
  # The same filter without resampling (straightforward SIS)
  # create matrices to store the particle values and weights
  n <- length(y)
  Tb <- Wb <- s <- k <- estim_state <- matrix(NA, ncol = num_of_particles, nrow = n) 
  
  # draw the particles for the initial state from the prior distribution
  Tb[1,] <- runif(num_of_particles, 0, 30^2)
  Wb[1,] <- 1/num_of_particles
  estim_state[1,] <- rep(0, num_of_particles)
  k[1,] <- rep(0, num_of_particles)
  s[1,] <- rep(1000, num_of_particles)
  
  for(t in 1:(n-1)) {
    k[t+1,] <- (s[t,] + Tb[t,]) / (s[t,] + Tb[t,] + 4^2)
    s[t+1,] <- (1 - k[t+1,])*(s[t,] + Tb[t,])
    estim_state[t+1,] <- estim_state[t,] + k[t+1,] * (y[t] - estim_state[t,])
    
    # sample particles according to the transition distribution
    Tb[t+1,] <- rtruncnorm(num_of_particles, a = 0, b = Inf, mean = Tb[t,], sd = 10)
    
    # compute the weights
    Wb[t+1,] <- dnorm(y[t+1], mean = estim_state[t,], sd = sqrt(s[t+1,] + Tb[t+1,] + 4^2)) * Wb[t,]
    Wb[t+1,] <- Wb[t+1,]/sum(Wb[t+1,])
    
    
  }
  
  return(rowSums(sqrt(Tb)*Wb))
  
}

est_vol_bootstrap <- pf_bootstrap(y = data$Observation, num_of_particles = 2000)
est_vol_straight <- pf_straight(y = data$Observation, num_of_particles = 2000)

# Plotting filters ----

plot(1:length(data$Volatility), data$Volatility, ylab="Value", xlab="Time point (t)", type = "l")
lines(est_vol_bootstrap, col = "red")
lines(est_vol_straight, col = "green")
legend(c(0, max(data$Volatility)), lty = c(1,1,1), col = c("black", "red", "green"), legend = c("Actual vol","Bootstrap filter (c=1)","SIS (c=0)"), bty="n")

# Simulating participant ----

estim_vol <- est_vol_bootstrap #or straight

participant_tracking_volatility <- function(estim_noise, estim_vol, observations) {
  eta <- 0
  pred <- 0
  n <- length(observations) - 1
  
  for (t in 1:n) {
    eta[t+1] <- estim_vol[t+1]^2 / (estim_vol[t+1]^2 + estim_noise)
    pred[t+1] <- pred[t] + eta[t+1] * (observations[t] - pred[t])
  }
  return(pred)
}

participant_ignoring_volatility <- function(estim_noise, learning_rate, observations) {
  eta <- learning_rate
  pred <- 0
  n <- length(observations) - 1
  for (t in 1:n) {
    pred[t+1] <- pred[t] + eta * (observations[t] - pred[t])
  }
  return(pred)
}

# pred <- participant_tracking_volatility(estim_noise = 4, estim_vol = est_vol_bootstrap, observations = data$Observation)
pred <- participant_ignoring_volatility(estim_noise = 4, learning_rate = 0.5, observations = data$Observation)

par(mfrow=c(1,1))
plot(data$Observation, type = "l", lty = 2)
lines(data$State)
lines(pred, col = "green")
legend(c(0, max(data$Observation)), lty = c(2,1,1), col = c("black", "black", "green"), legend = c("Observations","States","Predictions"), bty="n")


# Recovering model ----

set.seed(1)
data <- generate_data_volatiliy_block(num_of_trials = 100, noise = 4, block1_vol = 4, block2_vol = 16, start_state = 0)

set.seed(2)
est_vol_bootstrap <- pf_bootstrap(y = data$Observation, num_of_particles = 500)
simulated_participant_response <- participant_tracking_volatility(estim_noise = 4, estim_vol = est_vol_bootstrap, observations = data$Observation)

set.seed(3)
est_vol_bootstrap <- pf_bootstrap(y = data$Observation, num_of_particles = 500)
model_track_vol <- participant_tracking_volatility(estim_noise = 4, estim_vol = est_vol_bootstrap, observations = data$Observation)
mse_track_vol_model <- sum((model_track_vol - simulated_participant_response)^2)/100

mse_ignore_vol_model <- 0
n <- 1
for (i in seq(0, 1, 0.01)) {
  pred_ignore_vol <- participant_ignoring_volatility(estim_noise = 4, learning_rate = i, observations = data$Observation)
  mse_ignore_vol_model[n] <- sum((pred_ignore_vol - simulated_participant_response)^2)/100
  n = n + 1
}
mse_ignore_vol_model <- min(mse_ignore_vol_model)


