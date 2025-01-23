
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
  dimnames(popout) <- list(as.character(seq(1:nyears)), c('fawn.S.N', 'juv.f.S.N', 
  "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N", 
  "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N"))
  
  # for hunted
  hunted <- popout
  
  # for surveillance sampling 
  samp <- popout
  
  # for Environmental reservoir  
  envres <- matrix(0, nrow = nyears, ncol = 1)
  
  # For counting the number of farms that are closed 
  farmyears_noclosed <- matrix(0, nrow = nyears, ncol = 1)
  
  # For calibrating disease dynamics 
  disease <- matrix(0, nrow = nyears, ncol = 2)
  
  # For calibrating disease dynamics 
  disease <- matrix(0, nrow = nyears, ncol = 2)# Need to check popn vital rates (survival, repro, carrying capacity) 
  checkvitals <- matrix(0, nrow = nyears, ncol = 8)
  
  #jc revised 152-163
  disease.v <- matrix(0, nrow = nyears, ncol = 1)
  
  # starting population
  n0 <- runif(1, 120000, 160000)
  
  stablestate <- returnstartingpopn(params, n0) # this is ordered fawns, juv f, adults f, juv m, adult m
  
  popout[1,] <- c(stablestate[1:2], stablestate[4], stablestate[3], stablestate[5], rep(0, 9))
  
  hunted[1,] <- c(popout[1,1] * hunt.mort.fawn * (fawn.an.sur ^ (7/12)), popout[1,2] * hunt.mort.juv.f * (juv.f.an.sur ^ (7/12)), 
                  popout[1,3] * hunt.mort.juv.m * (juv.m.an.sur ^ (7/12)), popout[1,4] * hunt.mort.ad.f * (ad.an.f.sur ^ (7/12)), 
                  popout[1,5] * hunt.mort.ad.m * (ad.an.m.sur ^ (7/12)), rep(0, 9))
  
  
  ### Assuming diseased individuals are just added on top of pre-existing population 
  if(switchdiseasedy_on == 1){
    arrival <- c(burnin, arrival_input, rep(0, nyears - c(length(arrival_input) + length(burnin)) )) # arrival vec should be same length as number of years, this is just a fail safe
    # Because we have a burnin period, t1 should never have infected individuals.. 
    # The following code shouldn't ever have to run.. 
     if(arrival[1] == 0){ # if no new individuals arrive year 1, do nothing
    } else{
      split <- as.vector(table(sample(1:2, size = arrival[1], replace = T))) #assume equally likely btw males and females and adults and yearlings
      # if(length(split) <= 3){ 
      #   if(split[2] < 1){ split[2] <- 0}
      #   if(split[3] < 1){ split[3] <- 0}
      #   if(split[4] < 1){ split[4] <- 0}
      #   }# otherwise have issues
      popout[1,c(13:14)] <- split
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
  hunt.fawn.beta <- est_beta_params(hunt.mort.fawn, hunt.fawn.var) #jc addded fawn specific variance term
  hunt.juv.f.beta <- est_beta_params(hunt.mort.juv.f, hunt.var)
  hunt.juv.m.beta <- est_beta_params(hunt.mort.juv.m, hunt.var)
  hunt.f.beta <- est_beta_params(hunt.mort.ad.f, hunt.var)
  hunt.m.beta <- est_beta_params(hunt.mort.ad.m, hunt.var)
  
  ### Reproduction 
  juv.r.gamma <- est_gamma_params(juv.repro, repro.var) #jc removed the /2 on the mean repro param
  ad.r.gamma <- est_gamma_params(ad.repro, repro.var) #jc removed the /2 on the mean repro param
  
  #disease rates JC added 
  beta.d <- runif(1, beta.d.min, beta.d.max)
  beta.e <- runif(1, beta.e.min, beta.e.max)
  
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
      
  # yearly stochastic reproductive rates #jc changed to include 3 param model for cc
    if(switch_stoch_repro_rates == 1){       
      juv.repro.draw <- rgamma(1, juv.r.gamma$alpha, scale = juv.r.gamma$theta)*(1- (alpha.cc*(sum(popout[t-1,])/K)^beta.cc))
      ad.repro.draw <- rgamma(1, ad.r.gamma$alpha, scale = ad.r.gamma$theta)*(1- (alpha.cc*(sum(popout[t-1,])/K)^beta.cc))
      } else {
      juv.repro.draw <- juv.repro
      ad.repro.draw <- ad.repro
    }

    ### log vitals
    checkvitals[t,] <- c(fawn.sur, juv.f.sur, juv.m.sur, ad.f.sur, ad.m.sur, 
                            juv.repro.draw, ad.repro.draw, (sum(popout[t-1,])/K) ) 
    
    ### ### Disease dynamics
    # if(switchdiseasedy_on == 1){ # this is turned on in params
      
      # ### from individuals arriving from out of state
      # split <- as.vector(table(sample(1:4, size = arrival[t], replace = T)))
      # ### Need to make this conditional or get error when nothing new to add
      # if(sum(split) != 0){
      #   if(length(split) <= 3){ 
      #     if(split[2] < 1 | is.na(split[2])){ split[2] <- 0}
      #     if(split[3] < 1 | is.na(split[3])){ split[3] <- 0}
      #     if(split[4] < 1 | is.na(split[4])){ split[4] <- 0}
      #   }# otherwise have issues
      #   popout[t,c(11:14)] <- popout[t,c(11:14)] + split
      # }
      
      ### Update the environmental load 
      #JC updated some details but Annabelle needs to align this eqn with overleaf
      no.infect.indiv <- sum(popout[t-1,c(11:14)])
      
      envres[t,] <- (no.infect.indiv * shedrate_I) + (no.infect.indiv*0.2 * shedrate_H) +
        (no.infect.indiv*0.2 * shedrate_ND) - (exp(-expdecayconstant) * envres[t-1,])
      
      # transmission rate (or the rate of contact) beta and the probability of infection given that contact occurred
      e.rate <- ifelse(((((sum(popout[t-1,-c(1:10)]))/sum(popout[t-1,])) * beta.d) + (envres[t,] * beta.e))<1, 
                       ((((sum(popout[t-1,-c(1:10)]))/sum(popout[t-1,])) * beta.d) + (envres[t,] * beta.e)),
                       0.99)
  
    ### Susceptible Category 
    # Fawns
      popout[t,1] <- (popout[t-1,2] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (1-e.rate) ) + 
        (popout[t-1,4] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (1-e.rate)) +
        (popout[t-1,7] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (1-e.rate) ) + 
        (popout[t-1,9] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (1-e.rate)) +
        (popout[t-1,11] * (1- hunt.juv.f) * juv.repro.draw * juv.f.sur * (1-e.rate) ) + 
        (popout[t-1,13] * (1- hunt.ad.f) * ad.repro.draw * ad.f.sur * (1-e.rate)) 

    # 
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
      (popout[t-1,2] * (1 - hunt.juv.f) * juv.f.sur *  e.rate) +
      ( popout[t-1,9] * (1 - hunt.ad.f) * ad.f.sur * (1 - i.rate))

    popout[t,10] <- (popout[t-1,5] * (1- hunt.ad.m) * ad.m.sur * e.rate) +
      (popout[t-1,3] * (1 - hunt.juv.m) * juv.m.sur *  e.rate) +
      ( popout[t-1,10] * (1 - hunt.ad.m) * ad.m.sur * (1 - i.rate))
    

    ### Infected Category

    # Juveniles
    popout[t,11] <- (.5 * popout[t-1,6] * (1 - hunt.fawn) * fawn.sur)
    
    popout[t,12] <- (.5 * popout[t-1,6] * (1 - hunt.fawn) * fawn.sur)
    
    # Adults
    popout[t,13] <- (popout[t-1,9] * (1- hunt.ad.f) * ad.f.sur * i.rate) +
      (popout[t-1,7] * (1 - hunt.juv.f) * juv.f.sur * i.rate) +
      (popout[t-1,13] * (1- hunt.ad.f) * (ad.f.sur - (ad.f.sur*additive.dz.mort))) +
      (popout[t-1,11] * (1 - hunt.juv.f) * (juv.f.sur - (juv.f.sur*additive.dz.mort)))
    
    popout[t,14] <- (popout[t-1,10] * (1- hunt.ad.m) * ad.m.sur * i.rate) +
      (popout[t-1,8] * (1 - hunt.juv.m) * juv.m.sur * i.rate) +
      (popout[t-1,14] * (1- hunt.ad.m) * (ad.m.sur - (ad.m.sur*additive.dz.mort))) +
      (popout[t-1,12] * (1 - hunt.juv.m) * (juv.m.sur - (juv.m.sur*additive.dz.mort)))
    
    
    
    ### ### Disease dynamics
    if(switchdiseasedy_on == 1){ # this is turned on in params
      
      ### from individuals arriving from out of state
      split <- as.vector(table(sample(1:2, size = arrival[t], replace = T)))
      ### Need to make this conditional or get error when nothing new to add
      if(sum(split) != 0){
        # if(length(split) <= 3){ 
        #   if(split[2] < 1 | is.na(split[2])){ split[2] <- 0}
        #   if(split[3] < 1 | is.na(split[3])){ split[3] <- 0}
        #   if(split[4] < 1 | is.na(split[4])){ split[4] <- 0}
        # }# otherwise have issues
        popout[t,c(13:14)] <- popout[t,c(13:14)] + split
      }
    }



    ## Need count of number of dead/removed individuals
    hunted[t,] <- c(popout[t,1] * hunt.fawn * (fawn.an.sur), popout[t,2] * hunt.juv.f * (juv.f.an.sur), 
                    popout[t,3] * hunt.juv.m * (juv.m.an.sur), popout[t,4] * hunt.ad.f * (ad.an.f.sur), 
                    popout[t,5] * hunt.ad.m * (ad.an.m.sur), popout[t,6] * hunt.fawn * (fawn.an.sur), 
                    popout[t,7] * hunt.juv.f * (juv.f.an.sur), popout[t,8] * hunt.juv.m * (juv.m.an.sur), 
                    popout[t,9] * hunt.ad.f * (ad.an.f.sur), popout[t,10] * hunt.ad.m * (ad.an.m.sur), 
                    popout[t,11] * hunt.juv.f * (juv.f.an.sur), popout[t,12] * hunt.juv.m * (juv.m.an.sur), 
                    popout[t,13] * hunt.ad.f * (ad.an.f.sur), popout[t,14] * hunt.ad.m * (ad.an.m.sur))
    
#     ### Surveillance starts here.. 
#     # bio check 10% or less of total harvest - randomly select 10% of available harvest
#     # no. individuals available for sampling
#     sampfromhunt <- sum(hunted[t,])*.1
#     
#     ### Check if number of individuals available for sampling is more or less than target
#     # target number of deer to check = nosampled (set in params)
#     if(sampfromhunt < nosampled){
#       availableSamples <- sampfromhunt
#     } else{
#       availableSamples <- nosampled
#     }
# 
#     # randomly draw samples
#     samp[t, ] <- randomsampling(matrix = hunted[t,], target_sum = availableSamples, emptysample = samp[t, ])
#     ### need to save this to back calc detection ****
# 
#     
#     # see how many of the samples have infected individuals
#     # NOTE - assuming can only detect it in infected individuals (not including exposed)
#     no.found.infected <-   sum(samp[t, grep("\\.I\\.", colnames(samp), ignore.case = TRUE)])
#     #if this is more than one, trigger actions
#     
#     ### Actions
#     if(no.found.infected != 0){
#       if(Action_young_bucks == 1){
#         ### Target yearling bucks to reduce density demography [8 - 10%]
#         # Overwriting the season's harvest after bump harvest numbers
#         hunted[t, c(3,8,12)] <- round(hunted[t,c(3,8,12)]*1.1)
#       }
# 
# ### FIX ME!!!            
#       if(Action_lib_harvest == 1){ # produces warnings bc different lengths but functionally fine
#         ### Liberalize harvest season 
#         # totalpopn <- sum(popout[t,])
#         # totalhunted <- sum(hunted[t,])
#         # huntingtarget <- (15000/143758) * n0 # make 4x bigger 
#         # prop.popn.to.take <- ((totalhunted + huntingtarget)/totalpopn )
#         # 
#         # # want to pull a number btw 10 - 20% 
#         # huntertake <- sample(10:20, 1)*.01
#         # if(prop.popn.to.take >  huntertake){ # check to make sure not harvesting more of popn than random draw of hunter ability
#         #   if(totalhunted/totalpopn > huntertake){
#         #     huntingtarget <- 2 # effectively not increasing harvest 
#         #   }else{
#         #     #If harvesting more than draw, modify take to hit hunter draw
#         #     modify.take <- ((huntertake * totalpopn) - totalhunted)*(1/huntingtarget)
#         #     huntingtarget <- huntingtarget *modify.take
#         #   }
#         # }
#         
#       }
#     
#       if(Action_sharpshooting == 1){
#         ### Sharp shooting  
#         # Pull 300 per year until population drops to 5 deer/square mile 
#         # In the simulation that would be ?? 
#         # According to spreadsheet, target is approximately 15/square mile
#         # Assuming we have that with calibration
#         # would stop pulling when reduce population by 2/3s 
#         
#         ### Set up assumes that have more females than males.. 
#         stopvalue <- sum(popout[25,]) *1/3
#         #stopvalue = starting popn * 1/3 (taking burn in into account so start at t = 25)
#         
#         currentvalue <- sum(popout[t,])
#         
#         if(currentvalue > stopvalue){
#           checksize = currentvalue - stopvalue
#           if(checksize > 300){
#             setsize = 300
#           } else{
#             setsize = checksize
#           }
#           ## Need to pull 300 individuals but this should be random?
#           ## Can repurpose code used to pull hunting
#           
#           #Assuming sharpshooters go for an even split between ages and sexes (not true..) 
#           split <- as.vector(table(sample(1:4, size = setsize, replace = T)))
#           
#           ## For each of the infected classes, see how many can get
#           sharpshootout <- c(0,0,0,0) # reset every time period 
#           for(i in 1:4){ # for each class 
#             probadj <- popout[t,i + 10] / split[i] # adjust probability based on # infected available and how much sharpshooting effort you have
#             newprob <- (.75 * probadj)
#             if(newprob > 1){
#               newprob <- 1
#             }
#             sharpshootout[i] <- sum(rbinom(size = 1, n = popout[t,i + 10], prob  = newprob))
#             ## If not removing infected, then removing susceptible
#             othertake <- split - sharpshootout
#           }
#           
#           ### remove individuals removed by sharpshooting from gen pop
#           popout[t,11:14] <- popout[t,11:14] - sharpshootout
#           popout[t,2:5] <- popout[t,2:5] - othertake
#           ### Not adding to no harvested but could pull this out as separate data if want
#           
#         }
#       } 
#       
#       
#     
#     }
#     
#     
#     
#     
#     ### ### Farm years
#     ## Calculate TRUE prevalence 
#     prev <- sum(popout[t,c(6:14)])/sum(popout[t,])
#     
#     x1 = prev # prevalence in the wild 
#     x2 = action1 # double fencing (indicator)
#     x3 = action2 # movement or export of captive individuals (action: restriction)
# 
#     mu = alpha + beta1*x1 + beta2*x2 + beta3*x3
#     
#     p <- plogis(mu) # need logistic tranformation function (do by hand)
#     #logit(p) ~  rnorm(mu = a + bx1 + bx2… , sigma)
#     if(t == 2){
#       farmsleft = 10 #need to make dynamic
#     }
#     out <- rbinom(size = 1, n = farmsleft, prob  = p)
#     ## p value returned doesn't make sense... 
#     
#     ## Save count of farms shutdown
#     farmyears_noclosed[t] <- sum(out)
#     
#     ### Update number of farms for next calc
#     farmsleft <- farmsleft - sum(out)
#     if(farmsleft < 0){ # catch for when there are no more farms
#       farmsleft <- 0
#     }
#     
    
    
  } ## Close timeloop   
  
  output <- list(counts = popout, removed = hunted, environemnt = envres, farmyears_noclosed = farmyears_noclosed,
                 diseaseexplore = disease, detectionsamp = samp, checkvitals = checkvitals)

  return(output)

  }





