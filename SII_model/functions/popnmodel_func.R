
# Ref - https://clas.ucdenver.edu/marcelo-perraillon/sites/default/files/attached-files/lecture_12_inf_model.pdf

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

source("../VTCWD/SII_model/functions/CreateStableAgeStructure_func.R")


cwd_stoch_model <- function(params) {
  
  
  # Extract params
  # write the list objects to the local environment
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  
  ### Create list to hold output for each year
  popout <- matrix(0, nrow = nyears, ncol = 14)
  dimnames(popout) <- list(as.character(seq(1:nyears)), c('fawn.S.N', 'juv.f.S.N', "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N", "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N"))
  
  # for hunted
  hunted <- popout
  
  # for Environmental reservoir  
  envres <- matrix(0, nrow = nyears, ncol = 1)
  
  # CWD matrix 
  ### Need to set up a way to keep track of diseased individuals 
  
  stablestate <- returnstartingpopn() # this is ordered fawns, juv f, adults f, juv m, adult m
  
  popout[1,] <- c(stablestate[1:2], stablestate[4], stablestate[3], stablestate[5], rep(0, 9))
  hunted[1,] <- c(popout[1,1] * hunt.mort.fawn, popout[1,2] * hunt.mort.juv.f, popout[1,3] * hunt.mort.juv.m, popout[1,4] * hunt.mort.ad.f, popout[1,5] * hunt.mort.ad.m, rep(0, 9))
    


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
  juv.r.beta <- est_beta_params(juv.repro/2, repro.var)
  ad.r.beta <- est_beta_params(ad.repro/2, repro.var)
  ### !!!!! Is this right?? Do we want to divide by 2??? 
  
  for(t in 2:nyears){
    ## if there is too much inter annual variation, move this outside year
    fawn.sur <- rbeta(1, fawn.sur.beta$alpha, fawn.sur.beta$beta)
    juv.f.sur <- rbeta(1, juv.f.sur.beta$alpha, juv.f.sur.beta$beta)
    juv.m.sur <- rbeta(1, juv.m.sur.beta$alpha, juv.m.sur.beta$beta)
    ad.f.sur <- rbeta(1, ad.f.sur.beta$alpha, ad.f.sur.beta$beta)
    ad.m.sur <- rbeta(1, ad.m.sur.beta$alpha, ad.m.sur.beta$beta)
    
    # stochastic hunting survival rates
    hunt.fawn <- rbeta(1, hunt.fawn.beta$alpha, hunt.fawn.beta$beta, ncp = 0)
    hunt.juv.f <- rbeta(1, hunt.juv.f.beta$alpha, hunt.juv.f.beta$beta, ncp = 0)
    hunt.juv.m <- rbeta(1, hunt.juv.m.beta$alpha, hunt.juv.m.beta$beta, ncp = 0)
    hunt.ad.f <- rbeta(1, hunt.f.beta$alpha, hunt.f.beta$beta, ncp = 0)
    hunt.ad.m <- rbeta(1, hunt.m.beta$alpha, hunt.m.beta$beta, ncp = 0)
    
    # yearly stochastic reproductive rates
    
    
    juv.repro.draw <- rbeta(1, juv.r.beta$alpha, juv.r.beta$beta, ncp = 0)
    ad.repro.draw <- rbeta(1, ad.r.beta$alpha, ad.r.beta$beta, ncp = 0)

    ### Disease dynamics
    # transmission rate (or the rate of contact) beta and the probability of infection given that contact occurred

    e.rate <- (sum(popout[t,-c(1:5)])/sum(popout[t,])) *R0
    
    envres[t,] <- envres[t-1,] + (popout[t,c(11:14)]*shedrate) - (envres[t-1,] * expdecayconstant)
    
    ## envres needs to feed into indirect transmission
    # additive in beta? 
    # contact rate and contact given exposure 
    
    
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
    
    # ### Update the environmental load
    # no.infect.indiv <- sum(juv.f.I.N, juv.m.I.N, ad.f.I.N, ad.m.I.N)
    # env <- env + (no.infect.indiv * shedrate) - (expdecayconstant * env)
    # 
    # ### Update e
    # # This needs work!!! [after fitted to susceptible data..]
    # no.susceptible.indiv <- sum(juv.f.S.N, juv.m.S.N, ad.f.S.N, ad.m.S.N)
    # # beta?? 
    # # e = beta * no.infect.indiv * no.susceptible.indiv + env


    ## Need count of number of dead/removed individuals
    hunted[t,] <- c(popout[t,1] * hunt.fawn, popout[t,2] * hunt.juv.f, popout[t,3] * hunt.juv.m, popout[t,4] * hunt.ad.f, popout[t,5] * hunt.ad.m, 
                  popout[t,6] * hunt.fawn, popout[t,7] * hunt.juv.f, popout[t,8] * hunt.juv.m, popout[t,9] * hunt.ad.f, popout[t,10] * hunt.ad.m, 
                                                popout[t,11] * hunt.juv.f, popout[t,12] * hunt.juv.m, popout[t,13] * hunt.ad.f, popout[t,14] * hunt.ad.m)
    
  } ## Close timeloop   
  
  output <- list(counts = popout, removed = hunted)
  
  return(output)
}





