### Headers...
#
## Libraries
library(fitdistrplus)
library(patchwork)

### ### ### Set true distributions of parameters

confidence <- seq(from = 40, to= 99, by = 5)

set.seed(123)
### Normal: abundance 
truenorm <- rnorm(1000, mean = .30, sd = .05)

### Log-Normal: Length 

### Beta: Survival Prob

### Logist-Normal: Breeding Prob


### ### ### Draw samples from normal..
#### Draw 1 sample
drawconfidence = sample(confidence, 1)
samp1 <- sample(truenorm, drawconfidence) #need to randomly generate number of samples
samp1mean <- mean(samp1)
#min value = random draw from bottom quantile 
minquanSamp1 <- quantile(samp1, .25) # find quantile value
samp1[which(samp1 < minquanSamp1)] # subset 
lowestvalSamp1 <- sample(samp1[which(samp1 < minquanSamp1)], 1)

# max value = random draw from top quanitle 
maxquanSamp1 <- quantile(samp1, .75) # find quantile value
samp1[which(samp1 > maxquanSamp1)] # subset 
highestvalSamp1 <- sample(samp1[which(samp1 > maxquanSamp1)], 1)

#### automate
expertvaluesNorm <- empty.dump()
for(i in 1:7){ ## set number of experts here
  drawconfidence = sample(confidence, 1)
  samp1 <- sample(truenorm, drawconfidence) #need to randomly generate number of samples
  samp1mean <- mean(samp1)
  #min value = random draw from bottom quantile 
  minquanSamp1 <- quantile(samp1, .25) # find quantile value
  samp1[which(samp1 < minquanSamp1)] # subset 
  lowestvalSamp1 <- sample(samp1[which(samp1 < minquanSamp1)], 1)
  
  # max value = random draw from top quanitle 
  maxquanSamp1 <- quantile(samp1, .75) # find quantile value
  samp1[which(samp1 > maxquanSamp1)] # subset 
  highestvalSamp1 <- sample(samp1[which(samp1 > maxquanSamp1)], 1)
  ## lower, mean, highest, confidence
  sample <- c(i, lowestvalSamp1, samp1mean, highestvalSamp1, drawconfidence)
  if(i == 1){
    expertvaluesNorm <- sample
  }
  expertvaluesNorm <- rbind(expertvaluesNorm, sample)
  colnames(expertvaluesNorm) <- c("expert", "lowest", "mean", "highest", "confidence")
}

### fitting to distribution (note, could just combine the two loops..)
# note, x is determined by true distribution 
x <- seq(min(truenorm), max(truenorm), by = .01)
FittedvaluesNorm <- empty.dump()
FittedvaluesNorm_LOP <- empty.dump()
for(i in 1:7){ # for each expert... 
  values <- expertvaluesNorm[i, 2:4]
  pullProb <- expertvaluesNorm[i, 5]
  AdjustProb <- (100 - pullProb)/2
  lowerprob <- AdjustProb / 100
  higherprob <- (pullProb + AdjustProb) / 100
  est <- qmedist(values, "norm", probs = c(lowerprob, higherprob))
  params <- est$estimate
  
  ## For LOP
  prob = dnorm(x, params[1], params[2])
  expert <- rep(i, length(prob))
  expertdist <- cbind(x, prob, expert)
  
  ## For Vin
  hold_LOP <- cbind(params[1], params[2], i)
  
  if(i == 1){
    FittedvaluesNorm <- expertdist
    FittedvaluesNorm_LOP <- hold_LOP
  } else {
  ## For LOP
  FittedvaluesNorm <- rbind(FittedvaluesNorm, expertdist)
  colnames(FittedvaluesNorm) <- c("value", "prob", "expert")
  ## For VIN
  FittedvaluesNorm_LOP <- rbind(FittedvaluesNorm_LOP, hold_LOP)
  colnames(FittedvaluesNorm_LOP) <- c("mean", "sd", "expert")
  }
}

LOP_pdf <- FittedvaluesNorm %>%
  as.data.frame() %>% 
  group_by(value) %>% 
  summarize(prob = mean(prob))

vin_pdf <- data.frame(value = x, 
                      prob = dnorm(x, 
                                   mean(FittedvaluesNorm_LOP[,1]), 
                                   c(mean(FittedvaluesNorm_LOP[,2]))))

TRUE_pdf <- data.frame(value = x, 
                      prob = dnorm(x, 
                                   mean = .30, 
                                   sd = .05))

FittedvaluesNorm <- as.data.frame(FittedvaluesNorm)

both_pdf_norm <- ggplot(data = FittedvaluesNorm,
                   aes(x = value, y = prob)) +
  geom_line(aes(group = expert), color = "black", size = 0.5) +
  geom_line(data = LOP_pdf, aes(color = "Linear Opinion Pool"), size = 1) +
  geom_line(data = vin_pdf, aes(color = "Vincent average"), size = 1) +
  geom_line(data = TRUE_pdf, aes(color = "True Distribution"), size = 1) +
  labs(x = "bat wing span (mm)", 
       y = "cumulative probability")+
  scale_color_manual(values = c("#377EB8", "green", "#ff7f00"))+
  theme_classic()+
  theme(legend.position = c(0.2,0.9), 
        legend.title = element_blank())
both_pdf_norm


### ### Fit to the wrong distribution... Log-Normal


FittedvaluesNorm_log <- empty.dump()
FittedvaluesNorm_LOP_log <- empty.dump()
for(i in 1:7){ # for each expert... 
  values <- expertvaluesNorm[i, 2:4]
  pullProb <- expertvaluesNorm[i, 5]
  AdjustProb <- (100 - pullProb)/2
  lowerprob <- AdjustProb / 100
  higherprob <- (pullProb + AdjustProb) / 100
  est <- qmedist(values, "lnorm", probs = c(lowerprob, higherprob))
  params <- est$estimate
  
  ## For LOP
  prob = dlnorm(x, params[1], params[2])
  expert <- rep(i, length(prob))
  expertdist <- cbind(x, prob, expert)
  
  ## For Vin
  hold_LOP <- cbind(params[1], params[2], i)
  
  if(i == 1){
    FittedvaluesNorm_log <- expertdist
    FittedvaluesNorm_LOP_log <- hold_LOP
  } else {
    ## For LOP
    FittedvaluesNorm_log <- rbind(FittedvaluesNorm_log, expertdist)
    colnames(FittedvaluesNorm_log) <- c("value", "prob", "expert")
    ## For VIN
    FittedvaluesNorm_LOP_log <- rbind(FittedvaluesNorm_LOP_log, hold_LOP)
    colnames(FittedvaluesNorm_LOP_log) <- c("mean", "sd", "expert")
  }
}

LOP_pdf <- FittedvaluesNorm_log %>%
  as.data.frame() %>% 
  group_by(value) %>% 
  summarize(prob = mean(prob))

vin_pdf <- data.frame(value = x, 
                      prob = dlnorm(x, 
                                   mean(FittedvaluesNorm_LOP_log[,1]), 
                                   c(mean(FittedvaluesNorm_LOP_log[,2]))))

FittedvaluesNorm_log <- as.data.frame(FittedvaluesNorm_log)

both_pdf_lognorm <- ggplot(data = FittedvaluesNorm_log,
                   aes(x = value, y = prob)) +
  geom_line(aes(group = expert), color = "black", size = 0.5) +
  geom_line(data = LOP_pdf, aes(color = "Linear Opinion Pool"), size = 1) +
  geom_line(data = vin_pdf, aes(color = "Vincent average"), size = 1) +
  geom_line(data = TRUE_pdf, aes(color = "True Distribution"), size = 1) +
  labs(x = "bat wing span (mm)", 
       y = "cumulative probability")+
  scale_color_manual(values = c("#377EB8", "green", "#ff7f00"))+
  theme_classic()+
  theme(legend.position = "none")
both_pdf_lognorm



### ### Fit to the wrong distribution... Beta


FittedvaluesNorm_beta <- empty.dump()
FittedvaluesNorm_LOP_beta <- empty.dump()
for(i in 1:7){ # for each expert... 
  values <- expertvaluesNorm[i, 2:4]
  pullProb <- expertvaluesNorm[i, 5]
  AdjustProb <- (100 - pullProb)/2
  lowerprob <- AdjustProb / 100
  higherprob <- (pullProb + AdjustProb) / 100
  est <- qmedist(values, "beta", probs = c(lowerprob, higherprob))
  params <- est$estimate
  
  ## For LOP
  prob = dbeta(x, params[1], params[2])
  expert <- rep(i, length(prob))
  expertdist <- cbind(x, prob, expert)
  
  ## For Vin
  hold_LOP <- cbind(params[1], params[2], i)
  
  if(i == 1){
    FittedvaluesNorm_beta <- expertdist
    FittedvaluesNorm_LOP_beta <- hold_LOP
  } else {
    ## For LOP
    FittedvaluesNorm_beta <- rbind(FittedvaluesNorm_beta, expertdist)
    colnames(FittedvaluesNorm_beta) <- c("value", "prob", "expert")
    ## For VIN
    FittedvaluesNorm_LOP_beta <- rbind(FittedvaluesNorm_LOP_beta, hold_LOP)
    colnames(FittedvaluesNorm_LOP_beta) <- c("mean", "sd", "expert")
  }
}

LOP_pdf <- FittedvaluesNorm_beta %>%
  as.data.frame() %>% 
  group_by(value) %>% 
  summarize(prob = mean(prob))

vin_pdf <- data.frame(value = x, 
                      prob = dbeta(x, 
                                    mean(FittedvaluesNorm_LOP_beta[,1]), 
                                    c(mean(FittedvaluesNorm_LOP_beta[,2]))))

FittedvaluesNorm_beta <- as.data.frame(FittedvaluesNorm_beta)

both_pdf_beta <- ggplot(data = FittedvaluesNorm_beta,
                   aes(x = value, y = prob)) +
  geom_line(aes(group = expert), color = "black", size = 0.5) +
  geom_line(data = LOP_pdf, aes(color = "Linear Opinion Pool"), size = 1) +
  geom_line(data = vin_pdf, aes(color = "Vincent average"), size = 1) +
  geom_line(data = TRUE_pdf, aes(color = "True Distribution"), size = 1) +
  labs(x = "bat wing span (mm)", 
       y = "cumulative probability")+
  scale_color_manual(values = c("#377EB8", "green", "#ff7f00"))+
  theme_classic()+
  theme(legend.position = "none")
both_pdf_beta


both_pdf_norm + both_pdf_lognorm + both_pdf_beta



### ### Fit to the wrong distribution... Logit-normal.... NOT IN THE STATS PACKAGE...!!!


FittedvaluesNorm_logit <- empty.dump()
FittedvaluesNorm_LOP_logit <- empty.dump()
for(i in 1:7){ # for each expert... 
  values <- expertvaluesNorm[i, 2:4]
  pullProb <- expertvaluesNorm[i, 5]
  AdjustProb <- (100 - pullProb)/2
  lowerprob <- AdjustProb / 100
  higherprob <- (pullProb + AdjustProb) / 100
  est <- qmedist(values, "logitnorm", probs = c(lowerprob, higherprob))
  params <- est$estimate
  
  ## For LOP
  prob = dlogitnorm(x, params[1], params[2])
  expert <- rep(i, length(prob))
  expertdist <- cbind(x, prob, expert)
  
  ## For Vin
  hold_LOP <- cbind(params[1], params[2], i)
  
  if(i == 1){
    FittedvaluesNorm_logit <- expertdist
    FittedvaluesNorm_LOP_logit <- hold_LOP
  } else {
    ## For LOP
    FittedvaluesNorm_logit <- rbind(FittedvaluesNorm_logit, expertdist)
    colnames(FittedvaluesNorm_logit) <- c("value", "prob", "expert")
    ## For VIN
    FittedvaluesNorm_LOP_logit <- rbind(FittedvaluesNorm_LOP_logit, hold_LOP)
    colnames(FittedvaluesNorm_LOP_logit) <- c("mean", "sd", "expert")
  }
}

LOP_pdf <- FittedvaluesNorm_logit %>%
  as.data.frame() %>% 
  group_by(value) %>% 
  summarize(prob = mean(prob))

vin_pdf <- data.frame(value = x, 
                      prob = dlogitnorm(x, 
                                   mean(FittedvaluesNorm_LOP_log[,1]), 
                                   c(mean(FittedvaluesNorm_LOP_log[,2]))))

FittedvaluesNorm_log <- as.data.frame(FittedvaluesNorm_log)

both_pdf <- ggplot(data = FittedvaluesNorm_log,
                   aes(x = value, y = prob)) +
  geom_line(aes(group = expert), color = "black", size = 0.5) +
  geom_line(data = LOP_pdf, aes(color = "Linear Opinion Pool"), size = 1) +
  geom_line(data = vin_pdf, aes(color = "Vincent average"), size = 1) +
  geom_line(data = TRUE_pdf, aes(color = "True Distribution"), size = 1) +
  labs(x = "bat wing span (mm)", 
       y = "cumulative probability")+
  scale_color_manual(values = c("#377EB8", "green", "#ff7f00"))+
  theme_classic()+
  theme(legend.position = c(0.2,0.9), 
        legend.title = element_blank())
both_pdf








### ### ### fit distributions (ex from textbook)

# Ex abundance
## lower, median, and upper values
w <- c(33, 35, 39)
#80% confidence so specifiy lower 10 and upper 90 percentiles
est <- qmedist(w, "norm", probs = c(0.10, 0.90))
params <- est$estimate
params
test <- rnorm(1000, params[1], params[2])
quantile(test, c(0.1, .5, .9))

# Ex surival
## lower, median, and upper values
w <- c(.25, .55, .9)
#90% confidence so specifiy lower 5 and upper 95 percentiles
est <- qmedist(w, "beta", probs = c(0.05, 0.95))
params <- est$estimate
params
test <- rbeta(1000, params[1], params[2])
quantile(test, c(0.05, .5, .95))
