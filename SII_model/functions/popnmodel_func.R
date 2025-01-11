
#' Ref - https://clas.ucdenver.edu/marcelo-perraillon/sites/default/files/attached-files/lecture_12_inf_model.pdf

## Suseptible
#' fawn.S.N = Fawn Suseptible count (scaler value >= 0),
#' 
#' juv.f.S.N = Juvenile Female Suseptible count (scaler value >= 0),
#' 
#' juv.m.S.N = Juvenile Male Suseptible count (scaler value >= 0),
#' 
#' ad.f.S.N = Adult Female Suseptible count (scaler value >= 0),
#' 
#' ad.m.S.N = Adult Male Suseptible count (scaler value >= 0),

## Exposed
#' fawn.E.N = Fawn exposed count (scaler value >= 0),
#' 
#' juv.f.E.N = Juvenile Female exposed count (scaler value >= 0),
#' 
#' juv.m.E.N = Juvenile Male exposed count (scaler value >= 0),
#' 
#' ad.f.E.N = Adult Female exposed count (scaler value >= 0),
#' 
#' ad.m.E.N = Adult Male exposed count (scaler value >= 0),


## Infected
#' 
#' juv.f.I.N = Juvenile Female infected count (scaler value >= 0),
#' 
#' juv.m.I.N = Juvenile Male infected count (scaler value >= 0),
#' 
#' ad.f.I.N = Adult Female infected count (scaler value >= 0),
#' 
#' ad.m.I.N = Adult Male infected count (scaler value >= 0),


### ### survival --- 
# question: do we want different survival rates for infected age/sex classes or flat survival rate x mortality of being infected?*
#' fawn.an.sur = annual fawn survival (scaler value between 0 and 1),  
#' 
#' juv.f.an.sur = annual juvenile female survival (scaler value between 0 and 1),  
#' 
#' juv.m.an.sur = annual juvenile male survival (scaler value between 0 and 1),  
#' 
#' ad.an.f.sur = annual adult female survival (scaler value between 0 and 1),  
#' 
#' ad.an.m.sur = annual adult male survival (scaler value between 0 and 1),

### ### reproduction --- 
## female only 
#' 
#' juv.repro = juvenile reproduction (scaler value >= 0),  
#' 
#' ad.repro = adult reproduction (scaler value >= 0),


### ### Harvest --- 
#' hunt.mort.fawn = percentage of fawns hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.juv.f = percentage of juvenile females hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.juv.m = percentage of juvenile males hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.ad.f = percentage of adult females hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.ad.m = percentage of adult males hunted (scaler value between 0 and 1),


### ### disease state transition  --- 
## e = beta * no. infected * no susceptible + E
## i = probability of transitioning 

#' e.rate = rate of all age classes moving from susceptible to exposed (scaler value between 0 and 1),  
#' 
#' i.rate = rate of all age classes moving from susceptible to exposed (scaler value between 0 and 1), 
#' 

### ### Other parameters  --- 
# Environmental Load: E_{t-1} + no. infectious individuals * shedding rate - decay/loss of prions + additional prions added through other pathways (carcass disposal)
# Hunter population
#'  
#' env = % of the population that is infected by the environment per month (scaler value between 0 and 1)
#' 
#' shedrate = the rate at which individuals infected with CWD release prions over 1 year
#' 
#' expdecayconstant = rate at which prions are lost from the environment over 1 year
#' 
#' ?beta.f = female transmission coefficient (scaler value greater than 0),  
#' 
#' ?beta.m = male transmission coefficient (scaler value greater than 0), 
#' 
#' n0 = starting population value [we are assuming that no individuals are diseased when begin model?] 
#' 
#estimate beta distribution parameters from a mean and variance. 
est_beta_params <- function(mu, var) {
  # check that the var <mu(1-mu), if greater, use the largest possible
  if(var > mu*(1-mu)){var = mu*(1-mu)-.00001}
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

est_gamma_params <- function(mu, var) {
  # Check that variance is positive and not too large for given mean
  if (var <= 0) {
    stop("Variance must be positive.")
  }
  
  # Estimate alpha (shape) and theta (scale) parameters
  alpha <- mu^2 / var
  theta <- var / mu
  
  return(params = list(alpha = alpha, theta = theta))
}


cwd_stoch_model <- function(params) {

  
  # Extract params
  # write the list objects to the local environment
  for (v in 1:length(params)){
    assign(names(params)[v], params[[v]]) 
  }
  
  ### Create list to hold output for each year
  popout <- matrix(0, nrow = nyears, ncol = 14)
  dimnames(popout) <- list(as.character(seq(1:nyears)), c('fawn.S.N', 'juv.f.S.N', "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N", "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N"))
  
  # for hunted
  hunted <- popout
  # for surviellance sampling 
  samp <- popout
  
  # for Environmental reservoir  
  envres <- matrix(0, nrow = nyears, ncol = 1)
  
  # For calibrating disease dynamics 
  disease <- matrix(0, nrow = nyears, ncol = 2)
  
  # Need to check popn vital rates (survival, repro, carrying capacity) 
  checkvitals <- matrix(0, nrow = nyears, ncol = 8)

  
  # CWD matrix 
  ### Need to set up a way to keep track of diseased individuals 
  
  stablestate <- returnstartingpopn(params) # this is ordered fawns, juv f, adults f, juv m, adult m
  
  popout[1,] <- c(stablestate[1:2], stablestate[4], stablestate[3], stablestate[5], rep(0, 9))
  hunted[1,] <- c(popout[1,1] * hunt.mort.fawn, popout[1,2] * hunt.mort.juv.f, popout[1,3] * hunt.mort.juv.m, popout[1,4] * hunt.mort.ad.f, popout[1,5] * hunt.mort.ad.m, rep(0, 9))
    
  ### Assuming diseased individuals are just added on top of pre-existing population 
  if(switchdiseasedy_on == 1){
    arrival <- c(arrival_input, rep(0, nyears - length(arrival_input))) # arrival vec should be same length as number of years, this is just a fail safe
    if(arrival[1] == 0){ # if no new individuals arrive year 1, do nothing
    } else{
      split <- as.vector(table(sample(1:4, size = arrival[1], replace = T))) #assume equally likely btw males and females and adults and yearlings
      if(length(split) <= 3){ 
        if(split[2] < 1){ split[2] <- 0}
        if(split[3] < 1){ split[3] <- 0}
        if(split[4] < 1){ split[4] <- 0}
        }# otherwise have issues
      popout[1,c(11:14)] <- split
    }
  
  }


  ##estimate vitals##
  ### Survival
  fawn.sur.beta <- est_beta_params(fawn.an.sur, fawn.an.sur.var)
  juv.f.sur.beta <- est_beta_params(juv.f.an.sur, juv.f.an.sur.var)
  juv.m.sur.beta <- est_beta_params(juv.m.an.sur, juv.m.an.sur.var)
  ad.f.sur.beta <- est_beta_params(ad.an.f.sur, ad.an.f.sur.var)
  ad.m.sur.beta <- est_beta_params(ad.an.m.sur, ad.an.m.sur.var)
  ### Hunting 
  hunt.fawn.beta <- est_beta_params(hunt.mort.fawn, hunt.var)
  hunt.juv.f.beta <- est_beta_params(hunt.mort.juv.f, hunt.var)
  hunt.juv.m.beta <- est_beta_params(hunt.mort.juv.m, hunt.var)
  hunt.f.beta <- est_beta_params(hunt.mort.ad.f, hunt.var)
  hunt.m.beta <- est_beta_params(hunt.mort.ad.m, hunt.var)
  
  ### Reproduction 
  # juv.r.beta <- est_beta_params(juv.repro/2, repro.var)
  # ad.r.beta <- est_beta_params(ad.repro/2, repro.var)
  
  juv.r.gamma <- est_gamma_params(juv.repro/2, repro.var)
  ad.r.gamma <- est_gamma_params(ad.repro/2, repro.var)
  
  
  for(t in 2:nyears){
    if(switch_stoch_survival_rates == 1){
      ## if there is too much inter annual variation, move this outside year
      fawn.sur <- rbeta(1, fawn.sur.beta$alpha, fawn.sur.beta$beta)
      juv.f.sur <- rbeta(1, juv.f.sur.beta$alpha, juv.f.sur.beta$beta)
      juv.m.sur <- rbeta(1, juv.m.sur.beta$alpha, juv.m.sur.beta$beta)
      ad.f.sur <- rbeta(1, ad.f.sur.beta$alpha, ad.f.sur.beta$beta)
      ad.m.sur <- rbeta(1, ad.m.sur.beta$alpha, ad.m.sur.beta$beta)
      } else{
        ### Survival
        fawn.sur <- fawn.an.sur
        juv.f.sur <- juv.f.an.sur
        juv.m.sur <- juv.m.an.sur
        ad.f.sur <- ad.an.f.sur
        ad.m.sur <- ad.an.m.sur 
      }
    
    if(switch_stoch_hunt_rates == 1){ 
      # stochastic hunting survival rates
      hunt.fawn <- rbeta(1, hunt.fawn.beta$alpha, hunt.fawn.beta$beta, ncp = 0)
      hunt.juv.f <- rbeta(1, hunt.juv.f.beta$alpha, hunt.juv.f.beta$beta, ncp = 0)
      hunt.juv.m <- rbeta(1, hunt.juv.m.beta$alpha, hunt.juv.m.beta$beta, ncp = 0)
      hunt.ad.f <- rbeta(1, hunt.f.beta$alpha, hunt.f.beta$beta, ncp = 0)
      hunt.ad.m <- rbeta(1, hunt.m.beta$alpha, hunt.m.beta$beta, ncp = 0)
    }else{
      ### Hunting 
      hunt.fawn <- hunt.mort.fawn 
      hunt.juv.f <- hunt.mort.juv.f 
      hunt.juv.m <- hunt.mort.juv.m 
      hunt.ad.f <- hunt.mort.ad.f
      hunt.ad.m <- hunt.mort.ad.m
    }
      
      # yearly stochastic reproductive rates
      
      if(switch_stoch_repro_rates == 1){       
      
      # Generate random samples from the Gamma distribution
      juv.repro.draw <- rgamma(1, juv.r.gamma$alpha, scale = juv.r.gamma$theta)*2 * (1 - (sum(popout[t-1,])/K))
      ad.repro.draw <- rgamma(1, ad.r.gamma$alpha, scale = ad.r.gamma$theta)*2 * (1 - (sum(popout[t-1,])/K))
    } else {
      
      ### Reproduction 
      juv.repro.draw <- juv.repro
      ad.repro.draw <- ad.repro
    }

    
    ### log vitals
    checkvitals[t,] <- c(fawn.sur, juv.f.sur, juv.m.sur, ad.f.sur, ad.m.sur, 
                            juv.repro.draw, ad.repro.draw, (sum(popout[t-1,])/K) ) 
    
    if(t == 2){ # no movement on the first time step? Don't think this is right... 
      e.rate <- 0 
    }
    

    ### Susceptible Category 
    # Fawns
    popout[t,1] <- (popout[t-1,2] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (1-e.rate) ) + 
      (popout[t-1,4] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (1-e.rate)) +
      (popout[t-1,7] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (1-e.rate) ) + 
      (popout[t-1,9] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (1-e.rate)) +
      (popout[t-1,11] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (1-e.rate) ) + 
      (popout[t-1,13] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (1-e.rate)) 
    
    # Juveniles 
    popout[t,2] <- .5 * popout[t-1,1] * (1 - hunt.fawn) * fawn.sur * (1 - e.rate)
    
    popout[t,3] <- .5 * popout[t-1,1] * (1 - hunt.fawn) * fawn.sur * (1 - e.rate)
    
    # Adults 
    popout[t,4] <- (popout[t-1,4] * (1- hunt.ad.f) * ad.f.sur * (1-e.rate)) + 
      (popout[t-1,2] * (1 - hunt.juv.f) * juv.f.sur * (1-e.rate))
    
    popout[t,5] <- (popout[t-1,5] * (1- hunt.ad.m) * ad.m.sur * (1-e.rate)) + 
      (popout[t-1, 3] * (1 - hunt.juv.m) * juv.m.sur * (1-e.rate))
    
    # ### Exposed Category 
    # # Fawns
    popout[t,6] <- (popout[t-1,2] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (e.rate) ) + 
      (popout[t-1,4] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (e.rate)) +
      (popout[t-1,7] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (e.rate) ) + 
      (popout[t-1,9] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (e.rate)) +
      (popout[t-1,11] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (e.rate) ) + 
      (popout[t-1,13] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (e.rate)) 
    
    # # Juveniles 
    popout[t,7] <- .5 * popout[t-1,1] * (1 - hunt.fawn) * fawn.sur * e.rate + 
      (.5 * popout[t-1,6] * (1 - hunt.fawn) * fawn.sur * (1 - i.rate))
    
    popout[t,8] <- .5 * popout[t-1,1] * (1 - hunt.fawn) * fawn.sur * e.rate + 
      (.5 * popout[t-1,6] * (1 - hunt.fawn) * fawn.sur * (1 - i.rate))
    
    # # Adults 
    popout[t,9] <- (popout[t-1,4] * (1- hunt.ad.f) * ad.f.sur * e.rate) +
      ( popout[t-1,9] * (1 - hunt.ad.f) * ad.f.sur * (1 - i.rate))

    popout[t,10] <- (popout[t-1,5] * (1- hunt.ad.m) * ad.m.sur * e.rate) +
      ( popout[t-1,10] * (1 - hunt.ad.m) * ad.m.sur * (1 - i.rate))
    

    ### Infected Category

    # Juveniles
    popout[t,11] <- (.5 * popout[t-1,6] * (1 - hunt.fawn) * fawn.sur * i.rate)
    
    popout[t,12] <- (.5 * popout[t-1,6] * (1 - hunt.fawn) * fawn.sur * i.rate)
    
    # Adults
    popout[t,13] <- (popout[t-1,9] * (1- hunt.ad.f) * ad.f.sur * i.rate) +
      (popout[t-1,7] * (1 - hunt.juv.f) * juv.f.sur * i.rate) +
      (popout[t-1,13] * (1- hunt.ad.f) * ad.f.sur) +
      (popout[t-1,11] * (1 - hunt.juv.f) * juv.f.sur)
    
    popout[t,14] <- (popout[t-1,10] * (1- hunt.ad.m) * ad.m.sur * i.rate) +
      (popout[t-1,8] * (1 - hunt.juv.m) * juv.m.sur * i.rate) +
      (popout[t-1,14] * (1- hunt.ad.m) * ad.m.sur) +
      (popout[t-1,12] * (1 - hunt.juv.m) * juv.m.sur)
    
    
    ### ### Disease dynamics
    if(switchdiseasedy_on == 1){ # this is turned on in params
      
    ### from individuals arriving from out of state
    split <- as.vector(table(sample(1:4, size = arrival[t], replace = T)))
    ### Need to make this conditional or get error when nothing new to add
    if(sum(split) != 0){
      if(length(split) <= 3){ 
        if(split[2] < 1 | is.na(split[2])){ split[2] <- 0}
        if(split[3] < 1 | is.na(split[3])){ split[3] <- 0}
        if(split[4] < 1 | is.na(split[4])){ split[4] <- 0}
      }# otherwise have issues
      popout[t,c(11:14)] <- popout[t,c(11:14)] + split
    }
    
    ### Update the environmental load
    no.infect.indiv <- sum(popout[t,c(11:14)])
    envres[t,] <- (no.infect.indiv * shedrate_I) + (no.infect.indiv * shedrate_H) +
      (no.infect.indiv * shedrate_ND) +(exp(expdecayconstant) * envres[t-1,])

    
    if(model1 == 1){
      Fe = ke * log(1 + (lambda * ((Pe*envres[t,])/(ke*(1 + sigma * envres[t,])))))
      # !! What do we set parameters as?? 
    }else { #model 2 
      ### unclear what this should be - just looks like a single parameter beta?
    }
    
    ### Update e
    # transmission rate (or the rate of contact) beta and the probability of infection given that contact occurred
    e.rate <- ((sum(popout[t,-c(1:5)]))/sum(popout[t,])) *(sum(popout[t,c(1:5)])) * beta + Fe
    if(is.na(e.rate)){
      e.rate <- 0 
    }
    
    ## delete later - compare relationship between Fe and e
    disease[t,1] <- e.rate
    disease[t,2] <- Fe
    
    } else {
      e.rate <- 0 # for when disease dynamics aren't turned on 
    }


    ## Need count of number of dead/removed individuals
    hunted[t,] <- c(popout[t,1] * hunt.fawn, popout[t,2] * hunt.juv.f, popout[t,3] * hunt.juv.m, popout[t,4] * hunt.ad.f, popout[t,5] * hunt.ad.m, 
                  popout[t,6] * hunt.fawn, popout[t,7] * hunt.juv.f, popout[t,8] * hunt.juv.m, popout[t,9] * hunt.ad.f, popout[t,10] * hunt.ad.m, 
                                                popout[t,11] * hunt.juv.f, popout[t,12] * hunt.juv.m, popout[t,13] * hunt.ad.f, popout[t,14] * hunt.ad.m)
    
    
    ### Surveillance starts here.. 
    # bio check 10% or less of total harvest - randomly select 10% of available harvest
    # no. individuals available for sampling
    sampfromhunt <- sum(hunted[t,])*.1
    
    ### Check if number of individuals available for sampling is more or less than target
    # target number of deer to check = nosampled (set in params)
    if(sampfromhunt < nosampled){
      availableSamples <- sampfromhunt
    } else{
      availableSamples <- nosampled
    }

    # randomly draw samples
    samp[t, ] <- randomsampling(matrix = hunted[t,], target_sum = availableSamples, emptysample = samp[t, ])

    # see how many of the samples have infected individuals
    # NOTE - assuming can only detect it in infected individuals (not including exposed)
    no.found.infected <-   sum(samp[t, grep("\\.I\\.", colnames(samp), ignore.case = TRUE)])
    #if this is more than one, trigger actions
    
    
    ### ### Farm years
    ## Calculate TRUE prevalence 
    prev <- .2
    
    x1 = prev = 0 # prevalence in the wild 
    x2 = action1 = 1 # double fencing (indicator)
    x3 = action2 = 1 # movement or export of captive individuals (action: restriction)

    mu = alpha + beta1*x1 + beta2*x2 + beta3*x3
    
    p <- plogis(mu) # need logistic tranformation function (do by hand)
    #logit(p) ~  rnorm(mu = a + bx1 + bx2… , sigma)
    farmsleft = 10 #need to make dynamic
    out <- rbinom(size = 1, n = farmsleft, prob  = p)
    ## p value returned doesn't make sense... 
    
    
    
  } ## Close timeloop   
  
  output <- list(counts = popout, removed = hunted, environemnt = envres, diseaseexplore = disease,
                 checkvitals = checkvitals)
  
  return(output)
}





