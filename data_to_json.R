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

i <- a <- 0
sets <- "{"
q <- 1
for (i in 1:10) {
  set.seed(i)
  for (a in 0:1) {

    
    if (a == 0) {
      data <- generate_data_volatiliy_block(num_of_trials = 100, noise = 4, block1_vol = 3, block2_vol = 30, start_state = 0)
    } else {
      data <- generate_data_volatiliy_sinusoid(num_of_trials = 100, noise = 4, max = 30, min = 3, start_state = 0)
    }
    
    plot(1:100, data$Observation+100, type = "l", ylab = "Share price", xlab = "Day (t)")
    
    json_string <- paste("\"set", toString(q), "\": {", sep = "")
    n <- 0
    for (i in data$Observation[1:99]) {
      json_string <- paste(json_string, "\"", toString(n), "\":\"", toString(round((i-min(data$Observation))/(max(data$Observation)-min(data$Observation)), digits = 3)), "\",", sep = "");
      n <- n + 1
    }
    json_string <- paste(json_string, "\"", "99", "\":\"", toString(round((data$Observation[100]-min(data$Observation))/(max(data$Observation)-min(data$Observation)), digits = 3)), "\"}", sep = "")
    
    sets <- paste(sets, json_string, ",", sep = "")
    
    q <- q + 1
  }
}

sets <- substring(sets, 1, nchar(sets)-1)
sets <- paste(sets, "}", sep = "")
write(sets, file = "sets.json")


