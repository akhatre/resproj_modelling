library(parallel)
library(truncnorm)
library(lmerTest)
library(DEoptim)

#parsing data ----

library(jsonlite)
alldata <- fromJSON(txt="~/Downloads/resproj-27787-export (26).json") #worked on 3, 5, 9
excluded_players <- data.frame(id = c("TGvns0Ntq1au9I8kfO4Rx5GG5Io2", "UYdMLVpENYbwc9LWRt8YozxMZgv1", "ZTrBkwUb5tfAOtHS2xVNAVHe4Rg2", "SDVPRphxWbcOP4KuM5DFXvX0ZC22"), reason = c("under 18", "second attempt", "second attempt","porbably second attempt"))

ratios_sinusoid <- ratios_block <- data.frame(asd = rep(0,99))
rts_block_df <- rts_sinusoid_df <- predictions_block_df <- predictions_sinusoid_df <- locations_block_df <- locations_sinusoid_df <- data.frame(asd = rep(0,100))
users <- names(alldata$predictions) [!(names(alldata$predictions) %in% excluded_players$id)]
for (user in users) {
  set <- alldata$users[[user]]$set$value
  locations <- alldata$locations[[set]]
  predictions <- alldata$predictions[[user]]$value
  rts <- alldata$rts[[user]]$value
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
    rts_sinusoid_df[[user]] <- as.numeric(rts)
  } else {
    ratios_block[[user]] <- ratio[,1]
    predictions_block_df[[user]] <- as.numeric(predictions)
    locations_block_df[[user]] <- as.numeric(locations)
    rts_block_df[[user]] <- as.numeric(rts)
  }
}



volatility_block <- c(rep(3,24),rep(30,25),rep(3,25),rep(30,26))
volatility_sinusoid <- sin(seq(-0.5*pi, 2.5*pi, length.out = 100))*(27)/2 + (((27)/2) + 3)

long <- data.frame()
for (name in names(predictions_block_df)) {
  if (name == "asd") {
    next
  }
  long <- rbind(long, data.frame(id = name, prediction = predictions_block_df[,name], rts = rts_block_df[,name], update = c(diff(predictions_block_df[,name]),0), locations = locations_block_df[,name], volatility = volatility_block, set = "block"))
}

for (name in names(predictions_sinusoid_df)) {
  if (name == "asd") {
    next
  }
  long <- rbind(long, data.frame(id = name, prediction = predictions_sinusoid_df[,name], rts = rts_sinusoid_df[,name], update = c(diff(predictions_sinusoid_df[,name]),0), locations = locations_sinusoid_df[,name], volatility = volatility_sinusoid, set = "sinusoid"))
}

long$error <- long$locations - long$prediction

long$set <- factor(long$set)
contrasts(long$set) <- contr.helmert(2)

#analysing accuracy and rts ----

avrg_err <- aggregate(long$error, by = list (id = long$id), FUN = mean)
avrg_rts <- aggregate(long$rts, by = list (id = long$id), FUN = mean)

# avrg_err_outliers <- avrg_err$id[abs(avrg_err$x) > (mean(avrg_err$x) + 2*sd(avrg_err$x))]
# excluded_players <- rbind(excluded_players, data.frame(id = avrg_err_outliers, reason = rep("avrg error out of 2SD", length(avrg_err_outliers))))
avrg_rts_outliers <- avrg_rts$id[abs(avrg_rts$x) > (mean(avrg_rts$x) + 2*sd(avrg_rts$x))]
excluded_players <- rbind(excluded_players, data.frame(id = avrg_rts_outliers, reason = rep("avrg rt out of 2SD", length(avrg_rts_outliers))))

mse <- aggregate(long$error, by = list (id = long$id), FUN = function(x){sum(x^2)/100})
mse_outliers <- mse$id[abs(mse$x) > (mean(mse$x) + 2*sd(mse$x))]
excluded_players <- rbind(excluded_players, data.frame(id = mse_outliers, reason = rep("mse out of 2SD", length(mse_outliers))))

long <- subset(long, !(id %in% as.character(excluded_players$id)))

ratios_block <- ratios_block[,-c(asd,excluded_players$id[5:8])]
ratios_sinusoid <- ratios_sinusoid[,-c(asd,excluded_players$id[5:8])]


y <- rowMeans(ratios_block)
x <- 1:99
lo <- loess(y ~ x)
plot(x,y, ylim = c(-0.5,2), type = "l", ylab = "Learning rate", xlab = "Trial", col = "gray")
lines(predict(lo), col = "red")
# abline(v=c(25,50,75), col="black", lty = "dashed")
lines(x,rep(1,99), lty="dashed", col = "cyan")
lines(x,rep(0,99), lty="dashed", col = "cyan")
lines(x, volatility_block[1:99]/35, col = "blue")


y <- rowMeans(ratios_sinusoid)
x <- 1:99
lo <- loess(y ~ x)
plot(x,y, ylim = c(-5,5), type = "l", ylab = "Learning rate", xlab = "Trial", col = "gray")
lines(predict(lo), col = "red")
# abline(v=c(25,50,75), col="black", lty = "dashed")
lines(x,rep(1,99), lty="dashed", col = "cyan")
lines(x,rep(0,99), lty="dashed", col = "cyan")
lines(x, volatility_sinusoid[1:99]/15, col = "blue")

y <- rowMeans(data.frame(ratios_block,ratios_sinusoid))
x <- 1:99
lo <- loess(y ~ x)
plot(x,y, type = "l", ylim = c(-5,5), ylab = "Learning rate", xlab = "Trial")
lines(predict(lo), col = "red")
abline(v=c(25,50,75), col="purple")
lines(x, volatility_sinusoid[1:99]/15)
#scaling

s <- 100
long$prediction <- long$prediction*s
long$locations <- long$locations*s
long$update <- long$update*s
long$error <- long$error*s
long$volatility <- long$volatility*s

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
  plot(1:100, rts_block_df[[user]], type="l")
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
  plot(1:100, rts_sinusoid_df[[user]], type="l")
  sin_low[i] <- (sum(ratios[1:25]) + sum(ratios[51:75]))/50
  sin_high[i] <- (sum(ratios[26:50]) + sum(ratios[76:99]))/49
  i <- i + 1
}

sum(sin_low)/length(sin_low)
sum(sin_high)/length(sin_high)
plot(1:length(sin_low), sin_low, col = "blue", type = "l", ylim = c(0,5))
lines(1:length(sin_low), sin_high, col = "red")

#regression analysis (mixed effects model) ----



model <- lmer(update ~ error * volatility * set + (error|id), data = long)
summary(model)

model1 <- lmer(update ~ error * volatility * set + (error*volatility|id), data = long)
summary(model1)

model2 <- lmer(update ~ error + error : (volatility + set + volatility:set) + (error + error:volatility|id), data = long, REML = FALSE)
summary(model2)

model_null <- lmer(update ~ error + (error|id), data = long, REML = FALSE)
summary(model_null)

anova(model2, model_null)

long_grad_only <- long[long$set == "sinusoid",]
model2_for_grad_vol <- lmer(update ~ error + error:volatility + (error + error:volatility|id), REML = FALSE, data = long_grad_only)
model_for_grad_vol <- lmer(update ~ error * volatility + (error|id), data = long_grad_only, REML = FALSE)

#estimating learning rate from the model
long$set_as_number[long$set == "block"] <- -1
long$set_as_number[long$set == "sinusoid"] <- 1
long$eta <- 0.827  + 0.004 * long$volatility - 0.056 * long$set_as_number + 0.003 * long$set_as_number * long$volatility

plot(1:100, long$eta[long$id == unique(long$id)[1]], type = "l", ylim=c(0.7,1), col = "blue", ylab = "Learning rate", xlab = "Trial")
lines(1:100, long$eta[long$id == unique(long$id)[60]], col = "red")

ids <- unique(long$id)
etas <- data.frame(long$eta[long$id == ids[1]])
plot(1:100, etas[,1], type = "l", ylim = c(0.5,1))
for (n in 2:length(ids)) {
  etas <-  data.frame(etas, long$eta[long$id == ids[n]])
  lines(1:100, long$eta[long$id == ids[n]])
}
etas <- rowSums(etas)/length(names(etas))
plot(1:100, etas, type = "l")

#models ----
resample_systematic <- function(weights) {
  # input: weights is a vector of length N with (unnormalized) importance weights
  # output: a vector of length N with indices of the replicated particles
  N <- length(weights)
  wght <- weights
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

ll_pf_bootstrap <- function(par, ta0 = FALSE, y, r, num_of_particles = 1000, random_seed=12345) { #that also should take y, r, and num_of_pacticles as input variables. Though, num_of_particles is not to be optimsed at the moment.

  # set random seed
  set.seed(random_seed)
  
  # # set parameters
  # sd_y <- exp(par[1]) #noise of observations, real one was 4
  # sd_r <- exp(par[2]) #noise of participants' responses, unknown
  # logTa_sd0 <- exp(par[3])
  # logTa_mean0 <- par[4]
  # if(ta0) {
  #   sd_ta <- 0
  # } else {
  #   sd_ta <- exp(par[5]) #sd of particle update distribution
  # }
  
  # set parameters
  sd_y <- par[1] #noise of observations, real one was 4
  sd_r <- par[2] #noise of participants' responses, unknown
  logTa_sd0 <- par[3]
  logTa_mean0 <- par[4]
  if(ta0) {
    sd_ta <- 0
  } else {
    sd_ta <- par[5] #sd of particle update distribution
  }
  # other defaults
  estim_mean0 <- 50 # CHANGED was 0.5 should set this to the middle of the screen!
  estim_var0 <- 100 # CHANGED was 1
  
  # Bootstrap filter
  # create matrices to store the particle values and weights
  n <- length(y)
  logTa <- Ta <- Wa <- s <- k <- estim_state <- matrix(NA, ncol = num_of_particles, nrow = n)
  mu_out <- rep(0.0,n)
  
  # draw the particles for the initial state from the prior distribution
  logTa[1,] <- rnorm(num_of_particles, mean=logTa_mean0, sd=logTa_sd0)
  # derived other particle values
  Ta[1,] <- exp(logTa[1,])
  estim_state[1,] <- rep(estim_mean0, num_of_particles)
  k[1:2,] <- rep(0.0, num_of_particles)
  s[1:2,] <- rep(estim_var0, num_of_particles)
  # particle weights
  Wa[1,] <- 1/num_of_particles
  # save mu for output
  mu_out[1] <- sum(Wa[1,]*estim_state[1,])
  # loop over time
  for(t in 1:(n-1)) {
    
    # sample particles according to the transition distribution
    logTa[t+1,] <- rnorm(num_of_particles, mean=logTa[t,], sd=sd_ta)
    
    # compute derived other particle values
    Ta[t+1,] <- exp(logTa[t+1,])
    
    k[t+1,] <- (s[t,] + Ta[t,]) / (s[t,] + Ta[t,] + sd_y^2)
    s[t+1,] <- (1 - k[t+1,])*(s[t,] + Ta[t,])
    estim_state[t+1,] <- estim_state[t,] + k[t+1,] * (y[t] - estim_state[t,])
    
    # compute the weights
    Wa[t+1,] <- dnorm(y[t+1], mean = estim_state[t,], sd = sqrt(s[t+1,] + Ta[t+1,] + sd_y^2)) * Wa[t,] # could have multiplied by W, but not necessary with uniform weights
    Wa[t+1,] <- Wa[t+1,]/sum(Wa[t+1,]) # normalize
    
    if (is.nan(sum(Wa[t+1,])) || sum(Wa[t+1,] == 0)) {
      return(10000)
    }
    
    # save mu for output
    mu_out[t+1] <- sum(Wa[t+1,]*estim_state[t+1,])
    
    # draw indices of particles with systematic resampling
    asd <<- Wa[t+1,]
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

  logLik <- -2*sum(dnorm(r, mean=mu_out, sd=sd_r, log=TRUE))
  pf_vol <<- rowSums(sqrt(Ta)*Wa)
  return(logLik)
}

pure_kalman_filter <- function(par, y, r, random_seed=12345) {
  set.seed(random_seed)
  n = length(y)
  
  # set parameters
  sd_y <- par[1] #noise of observations, real one was 4
  sd_r <- par[2] #noise of participants' responses, unknown
  sd_vol <- par[3]
  
  estim_state <- 50
  k <- 0
  s <- 100
  for (t in 1:(n-1)) {
    k[t+1] <- (s[t] + sd_vol^2) / (s[t] + sd_vol^2 + sd_y^2)
    s[t+1] <- (1 - k[t+1])*(s[t] + sd_vol^2)
    estim_state[t+1] <- estim_state[t] + k[t+1] * (y[t] - estim_state[t])
  }
  logLik <- -2*sum(dnorm(r, mean=estim_state, sd=sd_r, log=TRUE))
  return(logLik)
}

delta_rule <- function(par, y, r, random_seed=12345) {
  set.seed(random_seed)
  n = length(y)
  
  # set parameters
  k <- par[1] #noise of observations, real one was 4
  sd_r <- par[2]
  
  estim_state <- 50
  
  for (t in 1:(n-1)) {
    estim_state[t+1] <- estim_state[t] + k * (y[t] - estim_state[t])
  }
  
  logLik <- -2*sum(dnorm(r, mean=estim_state, sd=sd_r, log=TRUE))
  return(logLik)
}

lme_model <- function(par, y, r, eta, random_seed=12345) {
  set.seed(random_seed)
  n = length(y)
  
  # set parameters
  sd_r <- par[1] #noise of observations, real one was 4
  
  estim_state <- 50
  
  for (t in 1:(n-1)) {
    estim_state[t+1] <- estim_state[t] + eta[t] * (y[t] - estim_state[t])
  }
  
  logLik <- -2*sum(dnorm(r, mean=estim_state, sd=sd_r, log=TRUE))
  return(logLik)
}


#serial optimisation----
set.seed(122334455)
random_seeds <- sample(1:10000000,length(unique(long$id)))
estimated_par <- list()
for(i in unique(long$id)) {
  data <- subset(long, id==i)
  estimated_par[[i]] <- optim(c(log(sqrt(.4)),log(sqrt(1)),log(sqrt(.1))), fn = ll_pf_bootstrap, y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)])
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

#parallel optimisation----
set.seed(122334455)
random_seeds <- sample(1:10000000,length(unique(long$id)))

prl_optimisation_btstrp_var_ta <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(c(log(sqrt(4)), log(sqrt(1)), log(sqrt(10)), 10, log(sqrt(10))), fn = ll_pf_bootstrap, y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)]))
  return(est_par)
}
prl_optimisation_btstrp_fixed_ta <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(c(log(sqrt(4)),log(sqrt(1)),log(sqrt(10)),10), ta0 = TRUE, fn = ll_pf_bootstrap, y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)]))
  return(est_par)
}
prl_optimisation_kalman <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(c(log(sqrt(4)),log(sqrt(1)),log(sqrt(10))), fn = pure_kalman_filter, y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)]))
  return(est_par)
}
prl_optimisation_delta <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(c(0.5,log(sqrt(1))), fn = delta_rule, y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)]))
  return(est_par)
}
prl_optimisation_lme <- function(i) {
  data <- subset(long, id==i)
  est_par <- as.character(unique(data$id))
  est_par <- c(est_par, optim(c(1), fn = lme_model, y = data$locations, r = data$prediction, eta = data$eta, random_seed = random_seeds[which(unique(long$id) == i)]))
  return(est_par)
}

no_cores <- detectCores()
cluster <- makeCluster(no_cores, type = "FORK")
optim_res_btstrp_var_ta <- parLapply(cluster, unique(long$id), fun = prl_optimisation_btstrp_var_ta)
optim_res_btstrp_fixed_ta <- parLapply(cluster, unique(long$id), fun = prl_optimisation_btstrp_fixed_ta)
optim_res_kalman <- parLapply(cluster, unique(long$id), fun = prl_optimisation_kalman)
optim_res_delta <- parLapply(cluster, unique(long$id), fun = prl_optimisation_delta)
optim_res_heuristic <- parLapply(cluster, unique(long$id), fun = prl_optimisation_heuristic)
stopCluster(cluster)

params_heuristic <- params_delta <- params_btstrp_fixed_ta <- params_btstrp_var_ta <- params_kalman <- data.frame(Deviance = rep(0, length(optim_res_btstrp_var_ta)))
for (i in 1:length(optim_res_btstrp_var_ta)) {
  params_btstrp_fixed_ta$Deviance[i] <- optim_res_btstrp_fixed_ta[[i]]$value
  params_btstrp_fixed_ta$AIC[i] <- optim_res_btstrp_fixed_ta[[i]]$value + 2*4
  params_btstrp_fixed_ta$BIC[i] <- optim_res_btstrp_fixed_ta[[i]]$value + log(100)*4
  params_btstrp_fixed_ta$sd_y[i] <- round(exp(optim_res_btstrp_fixed_ta[[i]]$par[1]), digits = 3)
  params_btstrp_fixed_ta$sd_r[i] <- round(exp(optim_res_btstrp_fixed_ta[[i]]$par[2]), digits = 3)
  params_btstrp_fixed_ta$logTa_sd[i] <- round(exp(optim_res_btstrp_fixed_ta[[i]]$par[3]), digits = 2)
  params_btstrp_fixed_ta$logTa_mean[i] <- round(optim_res_btstrp_fixed_ta[[i]]$par[4], digits = 2)
  params_btstrp_fixed_ta$id[i] <- optim_res_btstrp_fixed_ta[[i]][[1]]
  
  params_btstrp_var_ta$Deviance[i] <- optim_res_btstrp_var_ta[[i]]$value
  params_btstrp_var_ta$AIC[i] <- optim_res_btstrp_var_ta[[i]]$value + 2*5
  params_btstrp_var_ta$BIC[i] <- optim_res_btstrp_var_ta[[i]]$value + log(100)*5
  params_btstrp_var_ta$sd_y[i] <- round(exp(optim_res_btstrp_var_ta[[i]]$par[1]), digits = 3)
  params_btstrp_var_ta$sd_r[i] <- round(exp(optim_res_btstrp_var_ta[[i]]$par[2]), digits = 3)
  params_btstrp_var_ta$logTa_sd[i] <- round(exp(optim_res_btstrp_var_ta[[i]]$par[3]), digits = 2)
  params_btstrp_var_ta$logTa_mean[i] <- round(optim_res_btstrp_var_ta[[i]]$par[4], digits = 2)
  params_btstrp_var_ta$sd_ta[i] <- round(exp(optim_res_btstrp_var_ta[[i]]$par[5]), digits = 2)
  params_btstrp_var_ta$id[i] <- optim_res_btstrp_var_ta[[i]][[1]]
  
  params_kalman$Deviance[i] <- optim_res_kalman[[i]]$value
  params_kalman$AIC[i] <- optim_res_kalman[[i]]$value + 2*3
  params_kalman$BIC[i] <- optim_res_kalman[[i]]$value + log(100)*3
  params_kalman$sd_y[i] <- round(exp(optim_res_kalman[[i]]$par[1]), digits = 3)
  params_kalman$sd_r[i] <- round(exp(optim_res_kalman[[i]]$par[2]), digits = 3)
  params_kalman$sd_vol[i] <- round(exp(optim_res_kalman[[i]]$par[3]), digits = 3)
  params_kalman$id[[i]] <- optim_res_kalman[[i]][[1]]
  
  params_delta$Deviance[i] <- optim_res_delta[[i]]$value
  params_delta$AIC[i] <- optim_res_delta[[i]]$value + 2*2
  params_delta$BIC[i] <- optim_res_delta[[i]]$value + log(100)*2
  params_delta$k[i] <- optim_res_delta[[i]]$par[1]
  params_delta$sd_r[i] <- exp(optim_res_delta[[i]]$par[2])
  params_delta$id[i] <- optim_res_delta[[i]][[1]]
  
  params_heuristic$Deviance[i] <- optim_res_heuristic[[i]]$value
  params_heuristic$AIC[i] <- optim_res_heuristic[[i]]$value + 2*2
  params_heuristic$BIC[i] <- optim_res_heuristic[[i]]$value + log(100)*2
  params_heuristic$a[i] <- optim_res_heuristic[[i]]$par[1]
  params_heuristic$sd_r[i] <- exp(optim_res_heuristic[[i]]$par[2])
  params_heuristic$id[i] <- optim_res_heuristic[[i]][[1]]
}

#evolutionary optimisation ----
deoptim_res_btstrp_var_ta <- deoptim_res_btstrp_fixed_ta <- list()
ids <- unique(long$id)
for (i in ids[11:96]) {
  data <- subset(long, id == i)
  deoptim_res_btstrp_var_ta[[i]] <- DEoptim(fn = ll_pf_bootstrap, lower = c(0, 0, 0, log(0.001), 0), upper = c(10, 10, 10, 5, 30), y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
  deoptim_res_btstrp_fixed_ta[[i]] <- DEoptim(fn = ll_pf_bootstrap, lower = c(0, 0, 0, log(0.001)), upper = c(10, 10, 10, 5), y = data$locations, r = data$prediction, ta0 = TRUE, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

deoptim_res_delta <- list()
for (i in ids) {
  data <- subset(long, id == i)
  deoptim_res_delta[[i]] <- DEoptim(fn = delta_rule, lower = c(0, 0), upper = c(10, 10), y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

deoptim_res_lme <- list()
for (i in ids) {
  data <- subset(long, id == i)
  deoptim_res_lme[[i]] <- DEoptim(fn = lme_model, lower = c(0), upper = c(10), y = data$locations, r = data$prediction, eta = data$eta, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

deoptim_res_kalman <- list()
for (i in ids) {
  data <- subset(long, id == i)
  deoptim_res_kalman[[i]] <- DEoptim(fn = pure_kalman_filter, lower = c(0,0,0), upper = c(10,10,40), y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

params_btstrp_fixed_ta <- params_btstrp_var_ta <- data.frame(id = ids)
for(i in ids) {
  params_btstrp_var_ta$Deviance[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestval
  params_btstrp_var_ta$AIC[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestval + 2*5
  params_btstrp_var_ta$BIC[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestval +log(100)*5
  params_btstrp_var_ta$sd_y[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestmem[[1]]
  params_btstrp_var_ta$sd_r[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestmem[[2]]
  params_btstrp_var_ta$log_sd_ta[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestmem[[3]]
  params_btstrp_var_ta$log_mean_ta[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestmem[[4]]
  params_btstrp_var_ta$sd_ta[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_var_ta[[i]]$optim$bestmem[[5]]
  
  params_btstrp_fixed_ta$Deviance[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestval
  params_btstrp_fixed_ta$AIC[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestval + 2*4
  params_btstrp_fixed_ta$BIC[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestval +log(100)*4
  params_btstrp_fixed_ta$sd_y[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestmem[[1]]
  params_btstrp_fixed_ta$sd_r[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestmem[[2]]
  params_btstrp_fixed_ta$log_sd_ta[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestmem[[3]]
  params_btstrp_fixed_ta$log_mean_ta[which(params_btstrp_var_ta$id == i)] <- deoptim_res_btstrp_fixed_ta[[i]]$optim$bestmem[[4]]
}

params_delta <- params_lme <- data.frame(id = ids)
for(i in ids) {
  params_delta$Deviance[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestval
  params_delta$AIC[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestval + 2*2
  params_delta$BIC[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestval +log(100)*2
  params_delta$k[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestmem[[1]]
  params_delta$sd_r[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestmem[[2]]
  
  params_lme$Deviance[which(params_lme$id == i)] <- deoptim_res_lme[[i]]$optim$bestval
  params_lme$AIC[which(params_lme$id == i)] <- deoptim_res_lme[[i]]$optim$bestval + 2
  params_lme$BIC[which(params_lme$id == i)] <- deoptim_res_lme[[i]]$optim$bestval + log(100)
  params_lme$sd_r[which(params_lme$id == i)] <- deoptim_res_lme[[i]]$optim$bestmem[[1]]
}

params_kalman <- data.frame(id = ids)
for(i in ids) {
  params_kalman$Deviance[which(params_kalman$id == i)] <- deoptim_res_kalman[[i]]$optim$bestval
  params_kalman$AIC[which(params_kalman$id == i)] <- deoptim_res_kalman[[i]]$optim$bestval + 2*3
  params_kalman$BIC[which(params_kalman$id == i)] <- deoptim_res_kalman[[i]]$optim$bestval + log(100)*3
  params_kalman$sd_y[which(params_kalman$id == i)] <- deoptim_res_kalman[[i]]$optim$bestmem[[1]]
  params_kalman$sd_r[which(params_kalman$id == i)] <- deoptim_res_kalman[[i]]$optim$bestmem[[2]]
  params_kalman$sd_vol[which(params_kalman$id == i)] <- deoptim_res_kalman[[i]]$optim$bestmem[[3]]
}

#comparing models ----
models <- data.frame(model = c("params_btstrp_var_ta", "params_btstrp_fixed_ta", "params_kalman", "params_delta", "params_lme"))
for (model in models$model) {
  model_params <- get(model)
  other_models <- as.character(models$model[models$model!=model])
  row_index <- which(models$model == model)
  
  models$dAIC_mean[row_index] <- mean(params_delta$AIC - model_params$AIC)
  models$dAIC_sd[row_index] <- sd(params_delta$AIC - model_params$AIC)
  models$nAIC[row_index] <- table(model_params$AIC < get(other_models[1])$AIC & model_params$AIC < get(other_models[2])$AIC & model_params$AIC < get(other_models[3])$AIC & model_params$AIC < get(other_models[4])$AIC)["TRUE"]
 
  models$dBIC_mean[row_index] <- mean(params_delta$BIC - model_params$BIC)
  models$dBIC_sd[row_index] <- sd(params_delta$BIC - model_params$BIC)
  models$nBIC[row_index] <- table(model_params$BIC < get(other_models[1])$BIC & model_params$BIC < get(other_models[2])$BIC & model_params$BIC < get(other_models[3])$BIC & model_params$BIC < get(other_models[4])$BIC)["TRUE"]
  
  models$dDeviance_mean[row_index] <- mean(params_delta$Deviance - model_params$Deviance)
  models$dDeviance_sd[row_index] <- sd(params_delta$Deviance - model_params$Deviance)
  models$nDeviance[row_index] <- table(model_params$Deviance < get(other_models[1])$Deviance & model_params$Deviance < get(other_models[2])$Deviance & model_params$Deviance < get(other_models[3])$Deviance & model_params$Deviance < get(other_models[4])$Deviance)["TRUE"]
}


block_ids <- params_btstrp_var_ta$id %in% unique(long$id[long$set == "block"])
param <- "Deviance"

points <- points1 <- points2 <- points2 <- rep("red",96)
points[block_ids] <- "blue"
get_params <- function(param, sort = FALSE, model) {
  if (sort) {
    result <- sort(get(model)[[param]])
  } else {
    result <- get(model)[[param]]
  }
  return(result)
}

plot_params <- function(param_name, sort = FALSE, model_names, ylim = FALSE, xlim = FALSE) {
  set.seed(3)
  colors <- c("red","blue","green","grey","purple","orange","magenta","brown","pink","cyan")
  i <- 1
  par(cex = 1)
  for (model_name in model_names) {
    model <- get(model_name)
    n <- length(model$id)
    block_ids <- model$id %in% unique(long$id[long$set == "block"])
    points <- rep("0",n)
    points[-block_ids] <- colors[i]
    points[block_ids] <- colors[i+1]
    param <- get_params(param_name, sort = sort, model_name)
    if (ylim==FALSE && xlim==FALSE){ylim <- c(min(param),max(param)); xlim <- c(1,n)}
    if (which(model_name == model_names) == 1) {plot(1:n, param, col = points, ylim=ylim, xlim=xlim, type = "l")}
    else {lines(1:n, param, col = points)}
    i <- i + 2
  }
  par(cex = 0.5)
  legend("bottomright", legend = rep(model_names, each = 2), col = colors, lty = 1)
}

models_names <- c("params_btstrp_var_ta", "params_btstrp_fixed_ta", "params_kalman", "params_delta")
plot_params("AIC", model_names = models_names, sort = TRUE)

?plot(1:96, get_params, col = points)

points1[-block_ids] <- "green"
points1[block_ids] <- "magenta"
points(1:96, sort(params_btstrp_fixed_ta[[param]]), col = points1)

points2[-block_ids] <- "cyan"
points2[block_ids] <- "brown"
points(1:96, sort(params_kalman[[param]]), col = points2)

points3[-block_ids] <- "grey"
points3[block_ids] <- "orange"
points(1:96, sort(params_delta[[param]]), col = points3)

#debugging
test_id <- "05lfqZKIxbPaV1kS9Go8wXzeqVI2"
prl_optimisation_btstrp_var_ta(test_id)
prl_optimisation_btstrp_fixed_ta(test_id)
prl_optimisation_kalman(test_id)


plot(1:100, pure_kalman_filter(c(6.186476, 1.408480, 7.444722), y = long$locations[long$id == test_id], r = long$prediction[long$id == test_id]), type = "l")

ll_pf_for_all <- function(par) {
  ids <- unique(long$id)
  deviance <- 0
  for (i in 1:length(ids)) {
    data <- subset(long, id==ids[i])
    deviance[i] <- ll_pf_bootstrap(par, y = data$locations, r = data$prediction)
  }
  return(sum(deviance))
}




# watching learning rate ----
subject_pf_bootstrap <- function(par, ta0 = FALSE, y, num_of_particles = 1000, random_seed=12345) { #that also should take y, r, and num_of_pacticles as input variables. Though, num_of_particles is not to be optimsed at the moment.
  
  # set random seed
  set.seed(random_seed)
  
  # set parameters
  sd_y <- par[1] #noise of observations, real one was 4
  sd_r <- par[2] #noise of participants' responses, unknown
  logTa_sd0 <- par[3]
  logTa_mean0 <- par[4]
  if(ta0) {
    sd_ta <- 0
  } else {
    sd_ta <- par[5] #sd of particle update distribution
  }
  # other defaults
  estim_mean0 <- 50 # CHANGED was 0.5 should set this to the middle of the screen!
  estim_var0 <- 100 # CHANGED was 1
  
  # Bootstrap filter
  # create matrices to store the particle values and weights
  n <- length(y)
  logTa <- Ta <- Wa <- s <- k <- estim_state <- matrix(NA, ncol = num_of_particles, nrow = n)
  mu_out <- k_out <- ta_out <- rep(0.0,n)
  
  # draw the particles for the initial state from the prior distribution
  logTa[1,] <- rnorm(num_of_particles, mean=logTa_mean0, sd=logTa_sd0)
  # derived other particle values
  Ta[1,] <- exp(logTa[1,])
  estim_state[1,] <- rep(estim_mean0, num_of_particles)
  k[1:2,] <- rep(0.0, num_of_particles)
  s[1:2,] <- rep(estim_var0, num_of_particles)
  # particle weights
  Wa[1,] <- 1/num_of_particles
  # save mu for output
  mu_out[1] <- sum(Wa[1,]*estim_state[1,])
  # loop over time
  for(t in 1:(n-1)) {
    
    # sample particles according to the transition distribution
    logTa[t+1,] <- rnorm(num_of_particles, mean=logTa[t,], sd=sd_ta)
    
    # compute derived other particle values
    Ta[t+1,] <- exp(logTa[t+1,])
    
    k[t+1,] <- (s[t,] + Ta[t,]) / (s[t,] + Ta[t,] + sd_y^2)
    s[t+1,] <- (1 - k[t+1,])*(s[t,] + Ta[t,])
    estim_state[t+1,] <- estim_state[t,] + k[t+1,] * (y[t] - estim_state[t,])
    
    # compute the weights
    Wa[t+1,] <- dnorm(y[t+1], mean = estim_state[t,], sd = sqrt(s[t+1,] + Ta[t+1,] + sd_y^2)) * Wa[t,] # could have multiplied by W, but not necessary with uniform weights
    Wa[t+1,] <- Wa[t+1,]/sum(Wa[t+1,]) # normalize
    
    if (is.nan(sum(Wa[t+1,])) || sum(Wa[t+1,] == 0)) {
      return(10000)
    }
    
    # save mu for output
    mu_out[t+1] <- sum(Wa[t+1,]*estim_state[t+1,])
    k_out[t+1] <- sum(Wa[t+1,]*k[t+1,])
    ta_out[t+1] <- sum(Wa[t+1,]*sqrt(Ta[t+1,]))
    
    # draw indices of particles with systematic resampling
    asd <<- Wa[t+1,]
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
  
  pf_vol <<- rowSums(sqrt(Ta)*Wa)
  results <- data.frame(mu = mu_out, k = k_out, pfvol = pf_vol, ta = ta_out)
  return(results)
}
subject_kalman_filter <- function(par, y, random_seed=12345) {
  set.seed(random_seed)
  n = length(y)
  
  # set parameters
  sd_y <- par[1] #noise of observations, real one was 4
  sd_r <- par[2] #noise of participants' responses, unknown
  sd_vol <- par[3]
  
  estim_state <- 50
  k <- 0
  s <- 100
  for (t in 1:(n-1)) {
    k[t+1] <- (s[t] + sd_vol^2) / (s[t] + sd_vol^2 + sd_y^2)
    s[t+1] <- (1 - k[t+1])*(s[t] + sd_vol^2)
    estim_state[t+1] <- estim_state[t] + k[t+1] * (y[t] - estim_state[t])
  }
  results <- data.frame(k = k, mu = estim_state)
  return(results)
}
subject_heuristic_model <- function(par, y, random_seed=12345) {
  set.seed(random_seed)
  n = length(y)
  
  # set parameters
  a <- round(par[1])
  
  estim_state <- 50
  k <- c(1,1)
  
  e <- y[1] - estim_state[1]
  estim_state[2] <- estim_state[1] + k[2] * e[1]
  
  sigmoid_k <- 1/(1+exp(10:-10))
  
  for (t in 2:(n-1)) {
    e[t] <- y[t] - estim_state[t]
    
    if (abs(e[t]) < abs(e[t-1])) {
      if(a >= 1 & a < 21){a <- a + 1}
      k[t+1] <- sigmoid_k[a]
    } else {
      if(a > 1 & a <= 21){a <- a - 1}
      k[t+1] <- sigmoid_k[a]
    }
    
    estim_state[t+1] <- estim_state[t] + k[t+1] * e[t]
  }
  
  results <- data.frame(k = k, mu = estim_state)
  return(results)
}

colors = colors()[sample(1:length(colors()), length(ids), replace=FALSE)]

for (i in ids) {
  data <- subset(long, id == i)
  par <- as.numeric(params_btstrp_var_ta[which(params_btstrp_var_ta$id == i), c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE)])
  results <- subject_pf_bootstrap(par, y = data$locations, random_seed = random_seeds[which(unique(long$id) == i)])
  if (i == ids[1]) {
    # plot(1:100, results$k, type = "l", col = colors[which(ids == i)])
    pf_var_ta_etas <<- data.frame(results$k)
  } else {
    # lines(1:100, results$k, col = colors[which(ids == i)])
    pf_var_ta_etas <<- data.frame(pf_var_ta_etas, data.frame(results$k))
  }
}
plot(1:100, rowMeans(pf_var_ta_etas), type = "l", col = "red", ylim=c(0,1))
lines(1:100, rowMeans(pf_var_ta_etas)+rowSds(as.matrix(pf_var_ta_etas)), type = "l")
lines(1:100, rowMeans(pf_var_ta_etas)-rowSds(as.matrix(pf_var_ta_etas)), type = "l")

for (i in ids) {
  data <- subset(long, id == i)
  par <- as.numeric(params_btstrp_fixed_ta[which(params_btstrp_fixed_ta$id == i), c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE)])
  results <- subject_pf_bootstrap(par, y = data$locations, ta0 = TRUE)
  if (i == ids[1]) {
    plot(1:100, results$k, type = "l", col = colors[which(ids == i)])
  } else {
    lines(1:100, results$k, col = colors[which(ids == i)])
  }
}
for (i in ids) {
  data <- subset(long, id == i)
  par <- as.numeric(params_kalman[which(params_kalman$id == i), c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE)])
  results <- subject_kalman_filter(par, y = data$locations)
  if (i == ids[1]) {
    plot(1:100, results$k, type = "l", col = colors[which(ids == i)])
  } else {
    lines(1:100, results$k, col = colors[which(ids == i)])
  }
}
for (i in ids[1:5]) {
  data <- subset(long, id == i)
  par <- as.numeric(params_heuristic[which(params_heuristic$id == i), c(FALSE,FALSE,FALSE,TRUE,TRUE,FALSE)])
  results <- subject_heuristic_model(par, y = data$locations)
  if (i == ids[1]) {
    plot(1:100, results$k, type = "l", col = colors[which(ids == i)], ylim = c(0,1))
  } else {
    lines(1:100, results$k, col = colors[which(ids == i)])
  }
}

DEoptim(fn = heuristic_model, lower = c(2, log(0.5)), upper = c(20, log(20)), y = data$locations, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == test_id)])
