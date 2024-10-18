

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

#' e.fawn.rate = rate of fawns moving from susceptible to exposed (scaler value between 0 and 1),  
#' 
#' e.ja.rate = rate of juvenile and adults moving from susceptible to exposed (scaler value between 0 and 1),  
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



cwd_stoch_model <- function(params) {
  
  
  # Extract params
  # write the list objects to the local environment
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  
  # check if static parameters exist.    
  if(exists("fawn.an.sur")==FALSE){
    message("fawn survival is missing, using default value")
    fawn.an.sur <- 0.6
  }
  
  if(exists("juv.f.an.sur")==FALSE){
    message("juvenile survival is missing, using default value")
    juv.f.an.sur <- 0.8
  } 
  
  
  ### Create list to hold output for each year
  popout <- matrix(0, nrow = nyears, ncol = 14)
  dimnames(popout) <- list(as.character(seq(1:nyears)), c('fawn.S.N', 'juv.f.S.N', "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N", "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N"))
  
  # for hunted
  hunted <- popout
  
  # CWD matrix 
  ### Need to set up a way to keep track of diseased individuals 
  
  
  ### Create starting population - copied from Paul Cross model. Is there a better way to do this???
  # Create the Leslie Matrix to start the population at stable age dist
  M <- matrix(rep(0, 5*5), nrow = 5) #fawn, juv f, adult f, juv m, adult m
  
  # replace the -1 off-diagonal with the survival rates
  M[row(M) == (col(M) + 1)] <- c(juv.f.an.sur * (1 - hunt.mort.juv.f),
                                 (ad.an.f.sur *  (1 - hunt.mort.ad.f)), 0,
                                 (ad.an.m.sur * (1 - hunt.mort.ad.m)))
  M[4,1] <-                     (juv.m.an.sur * (1 - hunt.mort.juv.m)) # don't have male and female fawns so this is wonky
  # ****** Is this right??? 
  
  
  # if you want the top age category to continue to survive
  M[3, 3] <- ad.an.f.sur * (1 - hunt.mort.ad.f)
  M[5, 5] <- ad.an.m.sur * (1 - hunt.mort.ad.m)
  
  # insert the fecundity vector prebirth census - fawns don't have sex so only doing this once
  M[1, 1:3] <- c(0, juv.repro, ad.repro) * fawn.an.sur * (1 - hunt.mort.fawn)
  
  ##I don't like how this function is working right now (it's adding too many individuals into the juvenile category) so I am manually entering numbers below
  # Initializing with the stable age distribution.
  suseptible <- round(popbio::stable.stage(M)[1:5] * n0) # fawns, juv f, adults f, juv m, adult m
  
  popout[1,] <- c(42900, 15000, 15000, 61500, 13600, rep(0, 9))
  
  ##estimate vitals##
  fawn.sur.beta <- est_beta_params(fawn.an.sur, fawn.an.sur.var)
  juv.f.sur.beta <- est_beta_params(juv.f.an.sur, juv.f.an.sur.var)
  juv.m.sur.beta <- est_beta_params(juv.m.an.sur, juv.m.an.sur.var)
  ad.f.sur.beta <- est_beta_params(ad.an.f.sur, ad.an.f.sur.var)
  ad.m.sur.beta <- est_beta_params(ad.an.m.sur, ad.an.m.sur.var)
  
  for(t in 2:nyears){
    
    fawn.sur <- rbeta(1, fawn.sur.beta$alpha, fawn.sur.beta$beta)
    juv.f.sur <- rbeta(1, juv.f.sur.beta$alpha, juv.f.sur.beta$beta)
    juv.m.sur <- rbeta(1, juv.m.sur.beta$alpha, juv.m.sur.beta$beta)
    ad.f.sur <- rbeta(1, ad.f.sur.beta$alpha, ad.f.sur.beta$beta)
    ad.m.sur <- rbeta(1, ad.m.sur.beta$alpha, ad.m.sur.beta$beta)
    
    # #Pull out classes of individuals out of matrix and set as sex/age/disease state name
    # for (i in seq_along(colnames(popout))) {
    #   col_names <- colnames(popout)
    #   values <- popout[t-1,]
    #   assign(col_names[i], values[i])
    # }
    
    ### Susceptible Category 
    ##please look these over. I changed most of the relationships here##
    # Fawns
    popout[t,1] <- (popout[t-1,2] * (1- hunt.mort.juv.f) * juv.repro * juv.f.sur * (1-e.fawn.rate) ) + 
      (popout[t-1,4] * (1- hunt.mort.ad.f) * ad.repro * ad.f.sur * (1-e.fawn.rate)) +
      (popout[t-1,7] * (1- hunt.mort.juv.f) * juv.repro * juv.f.sur * (1-e.fawn.rate) ) + 
      (popout[t-1,9] * (1- hunt.mort.ad.f) * ad.repro * ad.f.sur * (1-e.fawn.rate)) +
      (popout[t-1,11] * (1- hunt.mort.juv.f) * juv.repro * juv.f.sur * (1-e.fawn.rate) ) + 
      (popout[t-1,13] * (1- hunt.mort.ad.f) * ad.repro * ad.f.sur * (1-e.fawn.rate)) 
    
    # Juveniles 
    popout[t,2] <- .5 * popout[t-1,1] * (1 - hunt.mort.fawn) * fawn.sur * (1 - e.fawn.rate)
    
    popout[t,3] <- .5 * popout[t-1,1] * (1 - hunt.mort.fawn) * fawn.sur * (1 - e.fawn.rate)
    
    # Adults 
    popout[t,4] <- (popout[t-1,4] * (1- hunt.mort.ad.f) * ad.f.sur * (1-e.ja.rate)) + 
      (popout[t-1,2] * (1 - hunt.mort.juv.f) * juv.f.sur * (1-e.ja.rate))
    
    popout[t,5] <- (popout[t-1,5] * (1- hunt.mort.ad.m) * ad.m.sur * (1-e.ja.rate)) + 
      (popout[t-1, 3] * (1 - hunt.mort.juv.m) * juv.m.sur * (1-e.ja.rate))
    
    # ### Exposed Category 
    # # Fawns
    # fawn.E.N[t, 6] <- fawn.S.N[t,2] * (1- hunt.mort.fawn) * fawn.an.sur * e.fawn.rate
    #   # (ad.f.S.N * (1- hunt.mort.juv.f) * juv.repro * juv.f.an.sur * e.fawn.rate ) +
    #   # (juv.f.E.N * (1- hunt.mort.juv.f) * juv.repro * juv.f.an.sur * e.fawn.rate ) + 
    #   # ( ad.f.E.N * (1- hunt.mort.juv.f) * juv.repro * juv.f.an.sur * e.fawn.rate ) +
    #   # (juv.f.I.N * (1- hunt.mort.juv.f) * juv.repro * juv.f.an.sur * e.fawn.rate ) + 
    #   # ( ad.f.I.N * (1- hunt.mort.juv.f) * juv.repro * juv.f.an.sur * e.fawn.rate )
    #   # 
    # # Juveniles 
    # juv.f.E.N[t,7] <- (.5 * fawn.S.N[t-1, 1] * (1 - hunt.mort.fawn) * fawn.an.sur * e.fawn.rate) + 
    #   (.5 * fawn.E.N[t-1, 6] * (1 - hunt.mort.fawn) * fawn.an.sur * (1 - i.fawn.rate)) +
    #   
    # 
    # juv.m.E.N <- (.5 * fawn.S.N * (1 - hunt.mort.fawn) * fawn.an.sur * e.fawn.rate) + 
    #   (.5 * fawn.E.N * (1 - hunt.mort.fawn) * fawn.an.sur * (1 - i.fawn.rate))
    # 
    # # Adults 
    # ad.f.E.N <- (ad.f.S.N * (1- hunt.mort.ad.f) * ad.an.f.sur * e.ja.rate) + 
    #   (ad.f.E.N * (1 - hunt.mort.ad.f) * ad.an.f.sur * (1 - i.ja.rate))
    # 
    # ad.m.E.N <- (ad.m.S.N * (1- hunt.mort.ad.m) * ad.an.m.sur * e.ja.rate) + 
    #   (ad.m.E.N * (1 - hunt.mort.ad.m) * ad.an.m.sur * (1 - i.ja.rate))
    # 
    # 
    # ### Infected Category 
    # 
    # # Juveniles 
    # juv.f.I.N <- (.5 * fawn.S.N * (1 - hunt.mort.fawn) * fawn.an.sur * i.ja.rate) 
    # 
    # juv.m.I.N <- (.5 * fawn.S.N * (1 - hunt.mort.fawn) * fawn.an.sur * i.ja.rate)
    # 
    # # Adults 
    # ad.f.I.N <- (ad.f.E.N * (1- hunt.mort.ad.f) * ad.an.f.sur * i.ja.rate) + 
    #   (juv.f.E.N * (1 - hunt.mort.juv.f) * juv.f.an.sur * i.ja.rate) + 
    #   (ad.f.I.N * (1- hunt.mort.ad.f) * ad.an.f.sur) + 
    #   (juv.f.I.N * (1 - hunt.mort.juv.f) * juv.f.an.sur) 
    # 
    # ad.m.I.N <- (ad.m.E.N * (1- hunt.mort.ad.m) * ad.an.m.sur * i.ja.rate) + 
    #   (juv.m.E.N * (1 - hunt.mort.juv.m) * juv.m.an.sur * i.ja.rate) + 
    #   (ad.m.I.N * (1- hunt.mort.ad.m) * ad.an.m.sur) + 
    #   (juv.m.I.N * (1 - hunt.mort.juv.m) * juv.m.an.sur) 
    # 
    # 
    # ### Update the environmental load
    # no.infect.indiv <- sum(juv.f.I.N, juv.m.I.N, ad.f.I.N, ad.m.I.N)
    # env <- env + (no.infect.indiv * shedrate) - (expdecayconstant * env)
    # 
    # ### Update e
    # # This needs work!!! [after fitted to susceptible data..]
    # no.susceptible.indiv <- sum(juv.f.S.N, juv.m.S.N, ad.f.S.N, ad.m.S.N)
    # # beta?? 
    # # e = beta * no.infect.indiv * no.susceptible.indiv + env
    # 
    # # group age/sex/disease stage population for t + 1 into one output 
    # counts <- list(
    #   fawn.S.N = fawn.S.N, juv.f.S.N = juv.f.S.N, juv.m.S.N = juv.m.S.N, ad.f.S.N = ad.f.S.N, ad.m.S.N = ad.m.S.N,
    #   fawn.E.N = fawn.E.N, juv.f.E.N = juv.f.E.N, juv.m.E.N = juv.m.E.N, ad.f.E.N = ad.f.E.N, ad.m.E.N = ad.m.E.N,
    #   juv.f.I.N = juv.f.I.N, juv.m.I.N = juv.m.I.N, ad.f.I.N = ad.f.I.N, ad.m.I.N = ad.m.I.N)
    # 
    # popout[t,] <- unlist(counts)
    # 
    # ## Need count of number of dead/removed individuals 
    # 
    # hunted_indiv <- list(
    #   fawn.S.N = fawn.S.N * hunt.mort.fawn, juv.f.S.N = juv.f.S.N * hunt.mort.juv.f, juv.m.S.N = juv.m.S.N * hunt.mort.juv.m, ad.f.S.N = ad.f.S.N * hunt.mort.ad.f, ad.m.S.N = ad.m.S.N * hunt.mort.ad.m,
    #   fawn.E.N = fawn.E.N * hunt.mort.fawn, juv.f.E.N = juv.f.E.N * hunt.mort.juv.f, juv.m.E.N = juv.m.E.N * hunt.mort.juv.m, ad.f.E.N = ad.f.E.N * hunt.mort.ad.f, ad.m.E.N = ad.m.E.N * hunt.mort.ad.m,
    #   juv.f.I.N = juv.f.I.N * hunt.mort.juv.f, juv.m.I.N = juv.m.I.N * hunt.mort.juv.m, ad.f.I.N = ad.f.I.N * hunt.mort.ad.f, ad.m.I.N = ad.m.I.N * hunt.mort.ad.m)
    # 
    # hunted[t,] <- unlist(hunted_indiv)
    
  } ## Close timeloop   
  
  output <- list(counts = popout)
  
  return(output)
}





