require(truncnorm)
set.seed(3)

# Generating the data ----

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

generate_data_fixed_vol <- function(num_of_trials, noise, vol, start_state) {
  data <- data.frame(Volatility = rep(vol, num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  for (t in 2:num_of_trials) {
    
    data$Volatility[t] <- vol
    sd = data$Volatility[t]
    data$State[t] <- data$State[t-1] + rnorm(1, mean = 0, sd = sd)
    data$Observation[t] <- data$State[t] + rnorm(1, mean = 0, sd = noise)
  }
  return(data)
}

# data <- generate_data_volatiliy_linear_increase(num_of_trials = 100, noise = 4, increase = 0.3, start_state = 0)
set.seed(123)

par(mfcol = c(3,2), oma = c(0,0,0,0), mar = c(4,4,1,1))
set.seed(123)
data <- generate_data_volatiliy_sinusoid(num_of_trials = 100, noise = 4, max = 30, min = 3, start_state = 0)
plot(1:100, data$Volatility, type = "l", ylab = "Volatility", xlab = "Trial",  col = "blue")
plot(1:100, (data$State-min(data$State))/(max(data$State)-min(data$State)), type = "l", ylab = "Mean fish location", xlab = "Trial",  col = "red")
plot(1:100, (data$Observation-min(data$Observation))/(max(data$Observation)-min(data$Observation)), type = "l", ylab = "Fish location", xlab = "Trial")

set.seed(123)
data <- generate_data_volatiliy_block(num_of_trials = 100, noise = 4, block1_vol = 3, block2_vol = 30, start_state = 0)
plot(1:100, data$Volatility, type = "l", ylab = "Volatility", xlab = "Trial",  col = "blue")
plot(1:100, (data$State-min(data$State))/(max(data$State)-min(data$State)), type = "l", ylab = "Mean fish location", xlab = "Trial",  col = "red")
plot(1:100, (data$Observation-min(data$Observation))/(max(data$Observation)-min(data$Observation)), type = "l", ylab = "Fish location", xlab = "Trial")



par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(4,3,1,1))
set.seed(123)
n <- 50
data <- generate_data_fixed_vol(n,5,1,0)
plot(1:n, data$Observation+100, lty = 2, type = "l", xlab = "Trial", ylab = "", yaxt = "n")
axis(side=2,labels=F) 
lines(1:n, data$State+100, ylab = "Share price", xlab = "Day (t)", col = "red")
lines(1:n, delta_subject(0.1, data$Observation+100), col = "blue")

set.seed(123)
n <- 50
data <- generate_data_fixed_vol(n,5,1,0)
plot(1:n, data$Observation+100, lty = 2, type = "l",  xlab = "Trial", ylab = "", yaxt = "n")
axis(side=2,labels=F) 
lines(1:n, data$State+100,ylab = "Share price", xlab = "Day (t)", col = "red")
lines(1:n, delta_subject(0.9, data$Observation+100), col = "blue")

set.seed(123)
n <- 50
data <- generate_data_fixed_vol(n,1,5,0)
plot(1:n, data$Observation+100, lty = 2, type = "l", xlab = "Trial", ylab = "", yaxt = "n")
axis(side=2,labels=F) 
lines(1:n, data$State+100, ylab = "Share price", xlab = "Day (t)", col = "red")
lines(1:n, delta_subject(0.9, data$Observation+100), col = "blue")

set.seed(123)
data <- generate_data_fixed_vol(n,1,5,0)
plot(1:n, data$Observation+100, lty = 2, type = "l", xlab = "Trial", ylab = "", yaxt = "n")
axis(side=2,labels=F) 
lines(1:n, data$State+100, ylab = "Share price", xlab = "Day (t)", col = "red")
lines(1:n, delta_subject(0.1, data$Observation+100), col = "blue")

delta_subject <- function(k, y) {
  
  #PARAMS
  n = length(y)
  estim_state <- 100
  
  #body
  for (t in 1:(n-1)) {
    estim_state[t+1] <- estim_state[t] + k * (y[t] - estim_state[t])
  }

  return(estim_state)
}

generate_data_behrens <- function(num_of_trials, noise, value1, value2, start_state) {
  data <- data.frame(Volatility = rep(NA, 0.5*num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  data$Volatility[1:0.5*num_of_trials] <- value1
  data$Volatility[1:0.5*num_of_trials] <- value1
  for (t in 2:num_of_trials) {
    if (runif(1, min = 0, max = 1) < 0.05) {
      data$State[t] <- data$State[t-1] + rnorm(1, mean = 0, sd = vol)
    } else {
      data$State[t] <- data$State[t-1]
    }
    data$Observation[t] <- data$State[t] + rnorm(1, mean = 0, sd = noise)
  }
  return(data)
}

generate_data_nassar <- function(num_of_trials, noise, vol, start_state) {
  data <- data.frame(Volatility = rep(NA, num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  for (t in 2:num_of_trials) {
    if (runif(1, min = 0, max = 1) < 0.04) {
      data$State[t] <- data$State[t-1] + rnorm(1, mean = 0, sd = vol)
    } else {
      data$State[t] <- data$State[t-1]
    }
    data$Observation[t] <- data$State[t] + rnorm(1, mean = 0, sd = noise)
  }
  return(data)
}

generate_data_volatiliy_2_block <- function(num_of_trials, noise, block1_vol, block2_vol, start_state) {
  data <- data.frame(Volatility = rep(block1_vol, num_of_trials), State = start_state, Observation = rnorm(1, mean = start_state, sd = noise))
  for (t in 2:num_of_trials) {
    
    if (t < num_of_trials/2) {
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

par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(4,3,1,1))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

set.seed(3)
n <- 600
data <- generate_data_nassar(n,5,5,0)
plot(1:n, data$State+100, type = "l", xlab = "Trial", ylab = "", col = "red", yaxt = "n")
axis(side=2,labels=F) 
# lines(1:n, data$State+100, ylab = "Share price", xlab = "Day (t)", col = "red")
# lines(1:n, delta_subject(0.1, data$Observation+100), col = "blue")

n <- 290
plot(1:n, c(rep(75,120), rep(20,40), rep(80,40), rep(20,30), rep(80,30), rep(20,30)), type = "l", xlab = "Trial", ylab = "", col = "red", yaxt = "n")
axis(side=2,labels=F) 
polygon(c(120,120,290,290), c(20, 80, 80, 20), col="azure2", border = NA)
polygon(c(0,0,120,120), c(20, 80, 80, 20), col="azure", border = NA)
lines(1:n, c(rep(75,120), rep(20,40), rep(80,40), rep(20,30), rep(80,30), rep(20,30)), col = "red")

set.seed(2)
n <- 100
data <- generate_data_volatiliy_2_block(n,10,7,30,0)
plot(1:n, data$State, type = "l", xlab = "Trial", ylab = "", col = "red", yaxt = "n")
axis(side=2,labels=F) 
polygon(c(50,50,100,100), c(min(data$State), max(data$State), max(data$State), min(data$State)), col="azure2", border = NA)
polygon(c(0,0,50,50), c(min(data$State), max(data$State), max(data$State), min(data$State)), col="azure", border = NA)
lines(1:n, data$State, type = "l", xlab = "Trial", ylab = "", col = "red")
