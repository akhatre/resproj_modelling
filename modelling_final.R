library(parallel)
library(truncnorm)
library(DEoptim)
library(jsonlite)
library(lmerTest)

#parsing the data ----
alldata <- fromJSON(txt="~/ENTER_YOUR_PATH/raw_data.json")
excluded_players <- data.frame(id = c("TGvns0Ntq1au9I8kfO4Rx5GG5Io2", "UYdMLVpENYbwc9LWRt8YozxMZgv1", "ZTrBkwUb5tfAOtHS2xVNAVHe4Rg2", "SDVPRphxWbcOP4KuM5DFXvX0ZC22"), reason = c("under 18", "second attempt", "second attempt","porbably second attempt"))
users <- names(alldata$predictions)[!(names(alldata$predictions) %in% excluded_players$id)]
volatility_block <- c(rep(3,24),rep(30,25),rep(3,25),rep(30,26))
volatility_sinusoid <- sin(seq(-0.5*pi, 2.5*pi, length.out = 100))*(27)/2 + (((27)/2) + 3)
long <- data.frame()

for (user in users) {
  set <- alldata$users[[user]]$set$value
  locations <- as.numeric(alldata$locations[[set]])
  predictions <- as.numeric(alldata$predictions[[user]]$value)
  rts <- as.numeric(alldata$rts[[user]]$value)
  if (as.numeric(substr(set,4,nchar(set)))%%2 == 0) {
    condition <- "sinusoid"
    volatility <- volatility_sinusoid
  } else {
    condition <- "block"
    volatility <- volatility_block
  }
  
  long <- rbind(long, data.frame(id = user, condition = condition, volatility = volatility, location = locations, prediction = predictions, update = c(diff(predictions),0), rts = rts))
}

long$error <- long$location - long$prediction
long$condition <- factor(long$condition)
contrasts(long$condition) <- contr.helmert(2)
long$condition_as_number[long$condition == "block"] <- -1
long$condition_as_number[long$condition == "sinusoid"] <- 1

avrg_rts <- aggregate(long$rts, by = list (id = long$id), FUN = mean)
avrg_rts_outliers <- avrg_rts$id[abs(avrg_rts$x) > (mean(avrg_rts$x) + 2*sd(avrg_rts$x))]
excluded_players <- rbind(excluded_players, data.frame(id = avrg_rts_outliers, reason = rep("avrg rt out of 2SD", length(avrg_rts_outliers))))

mse <- aggregate(long$error, by = list (id = long$id), FUN = function(x){sum(x^2)/100})
mse_outliers <- mse$id[abs(mse$x) > (mean(mse$x) + 2*sd(mse$x))]
excluded_players <- rbind(excluded_players, data.frame(id = mse_outliers, reason = rep("mse out of 2SD", length(mse_outliers))))

long <- subset(long, !(id %in% as.character(excluded_players$id)))
avrg_rts <- avrg_rts[avrg_rts$id %in% unique(long$id),]
mse <- mse[mse$id %in% unique(long$id),]

#scaling the data ----
scale0to100 <- function(){
  if (range(long$location)[2] == 1) {
    long$prediction <<- long$prediction*100
    long$location <<- long$location*100
    long$update <<- long$update*100
    long$error <<- long$error*100
  }
}
scale0to1 <- function(){
  if (range(long$location)[2] == 100) {
    long$prediction <<- long$prediction*0.01
    long$location <<- long$location*0.01
    long$update <<- long$update*0.01
    long$error <<- long$error*0.01
  }
}

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
pf <- function(par, y, r, random_seed, num_of_particles = 1000) {
  
  # PARAMS
  set.seed(random_seed)
  sd_y <- par[1]
  sd_ta <- par[2]
  
  logTa_mean0 <- 0
  estim_mean0 <- 50
  estim_var0 <- 100
  
  # screate matrices to store the particle values and weights
  n <- length(y)
  logTa <- Ta <- Wa <- s <- k <- estim_state <- matrix(NA, ncol = num_of_particles, nrow = n)
  mu_out <- rep(0.0,n)
  
  # draw the particles for the initial state from the prior distribution
  logTa[1,] <- rnorm(num_of_particles, mean=logTa_mean0, sd=sd_ta)
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
  
  logLik <- -2*sum(dnorm(r, mean=mu_out, sd=sd_y, log=TRUE))
  return(logLik)
}
delta <- function(par, y, r, random_seed) {
  
  #PARAMS
  set.seed(random_seed)
  n = length(y)

  k <- par[1]
  sd_r <- par[2]
  
  estim_state <- 50
  
  #body
  for (t in 1:(n-1)) {
    estim_state[t+1] <- estim_state[t] + k * (y[t] - estim_state[t])
  }
  
  logLik <- -2*sum(dnorm(r, mean=estim_state, sd=sd_r, log=TRUE))
  return(logLik)
}
kf <- function(par, random_seed, y, r) {
  
  set.seed(random_seed)
  n = length(y)
  sd_y <- par[1]
  sd_vol <- par[2]
  
  estim_state <- 50
  k <- 0
  s <- 100
  for (t in 1:(n-1)) {
    k[t+1] <- (s[t] + sd_vol^2) / (s[t] + sd_vol^2 + sd_y^2)
    s[t+1] <- (1 - k[t+1])*(s[t] + sd_vol^2)
    estim_state[t+1] <- estim_state[t] + k[t+1] * (y[t] - estim_state[t])
  }
  logLik <- -2*sum(dnorm(r, mean=estim_state, sd=sd_y, log=TRUE))
  return(logLik)
}
lme_delta <- function(par, i, random_seed) {
  
  #DATA
  y <- subset(long, id == i)$location
  r <- subset(long, id == i)$prediction
  v <- subset(long, id == i)$volatility
  c <- subset(long, id == i)$condition_as_number[1]
  
  rnd_error <- random_effects[i,2]
  rnd_erorr_volatility <- random_effects[i,3]
  
  fxd_error <- fixed_effects[[2]]
  fxd_error_volatility <- fixed_effects[[3]]
  fxd_error_condition <- fixed_effects[[4]]
  fxd_error_volatility_condition <- fixed_effects[[5]]
  
  #PARAMS
  set.seed(random_seed)
  n = length(y)
  
  sd_r <- par[1]
  
  estim_state <- 50
  k <- 1
  
  #body
  for (t in 1:(n-1)) {
    k[t+1] <- fxd_error + fxd_error_volatility*v[t] + fxd_error_condition*c + fxd_error_volatility_condition*v[t]*c + rnd_error + rnd_erorr_volatility*v[t]
    estim_state[t+1] <- estim_state[t] + k[t+1] * (y[t] - estim_state[t])
  }
  
  logLik <- -2*sum(dnorm(r, mean=estim_state, sd=sd_r, log=TRUE))
  return(logLik)
}

#lme models----
scale0to1()

lme_model_alt <- lmer(update ~ error + error : (volatility + condition + volatility:condition) + (error + error:volatility|id), data = long, REML = FALSE)
summary(lme_model_alt)

lme_model_null <- lmer(update ~ error + (error|id), data = long, REML = FALSE)
summary(lme_model_null)

random_effects <- ranef(lme_model_alt)[[1]]
fixed_effects <- fixef(lme_model_alt)

anova(lme_model_null, lme_model_alt)

#fitting models ----
scale0to100()
set.seed(1122334455)
ids <- unique(long$id)
random_seeds <- sample(1:10000000, length(ids))
num_of_cores <- detectCores()
cluster <- makeCluster(num_of_cores, type = "FORK")

deoptim_res_pf <- list()
for (i in ids) {
  print(paste(as.character(which(unique(long$id) == i)), "/96", sep = ""))
  data <- subset(long, id == i)
  deoptim_res_pf[[i]] <- DEoptim(fn = pf, lower = c(0, 0), upper = c(50, 4), y = data$location, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

deoptim_res_delta <- list()
for (i in ids) {
  print(paste(as.character(which(unique(long$id) == i)), "/96", sep = ""))
  data <- subset(long, id == i)
  deoptim_res_delta[[i]] <- DEoptim(fn = delta, lower = c(0, 0), upper = c(2, 50), y = data$location, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

deoptim_res_kf <- list()
for (i in ids) {
  print(paste(as.character(which(unique(long$id) == i)), "/96", sep = ""))
  data <- subset(long, id == i)
  deoptim_res_kf[[i]] <- DEoptim(fn = kf, lower = c(0, 0), upper = c(50, 50), y = data$location, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

deoptim_res_lme_delta <- list()
for (i in ids) {
  print(paste(as.character(which(unique(long$id) == i)), "/96", sep = ""))
  deoptim_res_lme_delta[[i]] <- DEoptim(fn = lme_delta, lower = c(0), upper = c(50), i = i, random_seed = random_seeds[which(unique(long$id) == i)], control = DEoptim.control(cluster = cluster))
}

params_pf <- params_delta <- params_kf <- params_lme_delta <- data.frame(id = ids)
for(i in ids) {
  params_kf$Deviance[which(params_kf$id == i)] <- deoptim_res_kf[[i]]$optim$bestval
  params_kf$AIC[which(params_kf$id == i)] <- deoptim_res_kf[[i]]$optim$bestval + 2*2
  params_kf$BIC[which(params_kf$id == i)] <- deoptim_res_kf[[i]]$optim$bestval +log(100)*2
  params_kf$sd_y[which(params_kf$id == i)] <- deoptim_res_kf[[i]]$optim$bestmem[[1]]
  params_kf$sd_vol[which(params_kf$id == i)] <- deoptim_res_kf[[i]]$optim$bestmem[[2]]
  
  params_pf$Deviance[which(params_pf$id == i)] <- deoptim_res_pf[[i]]$optim$bestval
  params_pf$AIC[which(params_pf$id == i)] <- deoptim_res_pf[[i]]$optim$bestval + 2*2
  params_pf$BIC[which(params_pf$id == i)] <- deoptim_res_pf[[i]]$optim$bestval +log(100)*2
  params_pf$sd_y[which(params_pf$id == i)] <- deoptim_res_pf[[i]]$optim$bestmem[[1]]
  params_pf$sd_ta[which(params_pf$id == i)] <- deoptim_res_pf[[i]]$optim$bestmem[[2]]
  
  params_delta$Deviance[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestval
  params_delta$AIC[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestval + 2*2
  params_delta$BIC[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestval +log(100)*2
  params_delta$k[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestmem[[1]]
  params_delta$sd_r[which(params_delta$id == i)] <- deoptim_res_delta[[i]]$optim$bestmem[[2]]
  
  params_lme_delta$Deviance[which(params_lme_delta$id == i)] <- deoptim_res_lme_delta[[i]]$optim$bestval
  params_lme_delta$AIC[which(params_lme_delta$id == i)] <- deoptim_res_lme_delta[[i]]$optim$bestval + 2*2
  params_lme_delta$BIC[which(params_lme_delta$id == i)] <- deoptim_res_lme_delta[[i]]$optim$bestval +log(100)*2
  params_lme_delta$sd_r[which(params_lme_delta$id == i)] <- deoptim_res_lme_delta[[i]]$optim$bestmem[[1]]
}

models <- data.frame(model = c("params_pf", "params_kf", "params_delta", "params_lme_delta"))
for (model in models$model) {
  model_params <- get(model)
  other_models <- as.character(models$model[models$model!=model])
  row_index <- which(models$model == model)
  
  models$dAIC_mean[row_index] <- mean(params_delta$AIC - model_params$AIC)
  models$dAIC_sd[row_index] <- sd(params_delta$AIC - model_params$AIC)
  models$nAIC[row_index] <- table(model_params$AIC < get(other_models[1])$AIC & model_params$AIC < get(other_models[2])$AIC & model_params$AIC < get(other_models[3])$AIC)["TRUE"]  # & model_params$AIC < get(other_models[4])$AIC
  
  models$dBIC_mean[row_index] <- mean(params_delta$BIC - model_params$BIC)
  models$dBIC_sd[row_index] <- sd(params_delta$BIC - model_params$BIC)
  models$nBIC[row_index] <- table(model_params$BIC < get(other_models[1])$BIC & model_params$BIC < get(other_models[2])$BIC & model_params$BIC < get(other_models[3])$BIC)["TRUE"] # & model_params$BIC < get(other_models[4])$BIC
  
  models$dDeviance_mean[row_index] <- mean(params_delta$Deviance - model_params$Deviance)
  models$dDeviance_sd[row_index] <- sd(params_delta$Deviance - model_params$Deviance)
  models$nDeviance[row_index] <- table(model_params$Deviance < get(other_models[1])$Deviance & model_params$Deviance < get(other_models[2])$Deviance & model_params$Deviance < get(other_models[3])$Deviance)["TRUE"] # & model_params$Deviance < get(other_models[4])$Deviance
}

#rounding
for (i in 2:length(names(models))) {
  models[,i] <- round(models[,i], digits = 2)
}

boxplot(params_pf$BIC, params_kf$BIC, params_delta$BIC, params_lme_delta$BIC, names = c("PF", "KF", "DR", "LME"), ylab = "BIC", xlab = "Model")

# verification of results ----
verification_clean_up <- function(ind) {
  if (exists("fitted_deviance") & exists("actual_deviance")) {
    rm(fitted_deviance,actual_deviance, envir = .GlobalEnv)
  }
  if (ind == 1) {
    scale0to100()
    failures <<- 0
  }
}
#pf
for (i in ids) {
  verification_clean_up(ind = which(i == ids))
  row <- which(params_pf$id == i)
  par <- as.numeric(params_pf[row, c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)])
  fitted_deviance <- as.numeric(params_pf[row, 2])
  data <- subset(long, id == i)
  actual_deviance <- pf(par = par, y = data$location, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)])
  
  if(round(fitted_deviance, 3) != round(actual_deviance, 3)) {
    print("FAILURE")
    print(fitted_deviance)
    print(actual_deviance)
    failures <- failures + 1
  } else {
    print(paste("FD:", round(fitted_deviance, 3), "| AD:", round(actual_deviance, 3), "| PARS:", round(par[1], 3), round(par[2], 3), "| y[1]:", data$location[1], "| r[1]:", data$prediction[1]))
  }
  if(i == ids[96]){print(paste("Failures count:", failures))}
}
#delta
for (i in ids) {
  verification_clean_up(ind = which(i == ids))
  row <- which(params_lme_delta$id == i)
  sd_r <- params_lme_delta[row, 5]
  fitted_deviance <- as.numeric(params_lme_delta[row, 2])
  actual_deviance <- lme_delta(par = sd_r, i = i, random_seed = random_seeds[which(unique(long$id) == i)])

  if(round(fitted_deviance, 3) != round(actual_deviance, 3)) {
    print("FAILURE")
    print(fitted_deviance)
    print(actual_deviance)
    failures <- failures + 1
  } else {
    print(paste("FD:", round(fitted_deviance, 3), "| AD:", round(actual_deviance, 3), "| PARS:", round(sd_r[1], 3), "| i:", i))
  }
  if(i == ids[96]){print(paste("Failures count:", failures))}
}
#delta
for (i in ids) {
  verification_clean_up(ind = which(i == ids))
  row <- which(params_delta$id == i)
  par <- as.numeric(params_delta[row, c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)])
  data <- subset(long, id == i)
  fitted_deviance <- as.numeric(params_delta[row, 2])
  actual_deviance <- delta(par = par, y = data$location, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)])

  if(round(fitted_deviance, 3) != round(actual_deviance, 3)) {
    print("FAILURE")
    print(fitted_deviance)
    print(actual_deviance)
    failures <- failures + 1
  } else {
    print(paste("FD:", round(fitted_deviance, 3), "| AD:", round(actual_deviance, 3), "| PARS:", round(par[1], 3), round(par[2], 3), "| y[1]:", data$location[1], "| r[1]:", data$prediction[1]))
  }
  if(i == ids[96]){print(paste("Failures count:", failures))}
}
#kf
for (i in ids) {
  verification_clean_up(ind = which(i == ids))
  row <- which(params_kf$id == i)
  par <- as.numeric(params_kf[row, c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)])
  data <- subset(long, id == i)
  fitted_deviance <- as.numeric(params_kf[row, 2])
  actual_deviance <- kf(par = par, y = data$location, r = data$prediction, random_seed = random_seeds[which(unique(long$id) == i)])

  if(round(fitted_deviance, 3) != round(actual_deviance, 3)) {
    print("FAILURE")
    print(fitted_deviance)
    print(actual_deviance)
    failures <- failures + 1
  } else {
    print(paste("FD:", round(fitted_deviance, 3), "| AD:", round(actual_deviance, 3), "| PARS:", round(par[1], 3), round(par[2], 3), "| y[1]:", data$location[1], "| r[1]:", data$prediction[1]))
  }
  if(i == ids[96]){print(paste("Failures count:", failures))}
}

