

## Additional pieces... 

#set up arrival piece as vector.. 
# currently set up so that have arrival for as long as model is running but 
#   if need to project out, might need to change this.. 

# arrival_input <- c(100,100,1,1,1,1,5, 1, 8, 14) # 10 years.. 
# 
# # number of individuals sampled is determined by strategy 
# # statquosamp <- 5
# # prv2.5_90samp <- 125
# # nosampled <- statquosamp
# # prv1_99samp <- 460
# # prv1_90samp <- 230
# 
# huntingactions <- 1

#' CWD stochastic model function
#'
#' Stochastic monthly age and sex structured model with constant
#'  environmental transmission and dynamic direct transmission. The function 
#'  conducts one run of the model. 
#'
#' @param params A list with the following parameters included: 
#' 
#' fawn.an.sur = annual fawn survival (scaler value between 0 and 1),  
#' 
#' juv.an.sur = annual juvenile survival (scaler value between 0 and 1),  
#' 
#' ad.an.f.sur = annual adult female survival (scaler value between 0 and 1),  
#' 
#' ad.an.m.sur = annual adult male survival (scaler value between 0 and 1),   
#' 
#' fawn.repro = fawn reproduction (scaler value >= 0),  
#' 
#' juv.repro = juvenile reproduction (scaler value >= 0),  
#' 
#' ad.repro = adult reproduction (scaler value >= 0),  
#' 
#' hunt.mort.fawn = percentage of fawns hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.juv.f = percentage of juvenile females hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.juv.m = percentage of juvenile males hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.ad.f = percentage of adult females hunted (scaler value between 0 and 1),  
#' 
#' hunt.mort.ad.m = percentage of adult males hunted (scaler value between 0 and 1),  
#' 
#' ini.fawn.prev = percentage of fawns infected at the start (scaler value between 0 and 1),  
#' 
#' ini.juv.prev = percentage of juveniles infected at the start (scaler value between 0 and 1),  
#' 
#' ini.ad.f.prev = percentage of adult females infected at the start (scaler value between 0 and 1),  
#' 
#' ini.ad.m.prev = percentage of adult males infected at the start (scaler value between 0 and 1),   
#' 
#' n.age.cats = number of age categories to monitor. (scaler value greater than 3).
#' The final age category includes all those of that age or greater. 
#' 
#' p = rate of movement in the infectious categories (scaler values between 0 and 1). 
#' See model documentation vignette for how this relates to disease induced mortality.   
#' 
#' env.foi = % of the population that is infected by the environment per month (scaler value between 0 and 1)
#' 
#' beta.f = female transmission coefficient (scaler value greater than 0),  
#' 
#' beta.m = male transmission coefficient (scaler value greater than 0),  
#' 
#' theta = effect of population size on transmission (1 = frequency dependence, 0 = density dependent).  
#' 
#' n0 = initial population size (scaler value greater than 0)
#' 
#' n.years = number of years to run the model (scaler value greater than 2),  
#' 
#' rel.risk = relative risk of infected individuals being hunted. A value of 1 
#' indicates no hunter preference for infected individuals  
#' 
#' repro.var = variance in the annual reproduction rates from year to year
#' 
#' fawn.sur.var = variance in the annual fawn survival rate from year to year
#' 
#' sur.var = variance in the juvenile and adult survival rates from year to year. 
#' 
#' hunt.var = variance in the proportion of individuals hunter per category per year. 
#'
#'
#' @return A list with 2 dataframes is returned as output: 
#' 
#' 1. counts of the # of individuals in the susceptible and infectious 
#' categories by over time. 
#' 
#'  Columns include: 
#' 
#'  age (in years)
#' 
#'  month of simulation,
#' 
#'  population = number of individuals
#' 
#'  category: St.f = susceptible females, St.m = susceptible males, Ixt.f = 
#'  infectious females in the x category (1-10), Ixt.m = infectious males in the 
#'  x infectious category (1-10) 
#' 
#'  sex = female or male
#' 
#' disease = yes or no for susceptible or infectious
#' 
#'   
#' 2. deaths--how individuals died over time (hunting, natural or disease).
#' 
#'  Columns include: 
#'  
#'  age in years, 
#'  
#'  month of the simulation, 
#'  
#'  population = # of individuals, 
#'  
#'  category: Ht.f = hunted females, Ht.m = hunted males, Dt.f = natural 
#'  mortality females, Dt.m = natural mortality males, CWDt.f = disease mortality 
#'  females, CWDt.m = disease mortality males.  
#'    
#'  year = year of the simulation
#'  
#'  sex
#'  
#' 3. m.R0 = basic disease reproductive number for an initially infected adult 
#' male for direct transmission only.
#'   
#' 4. f.R0 = basic disease reproductive number for an initially infected adult 
#' female for direct transmission only.
#' 
#' @importFrom popbio stable.stage
#' @importFrom stats rbeta rbinom rnbinom rgamma
#' @importFrom dplyr rename mutate
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @examples 
#' 
#' params <- list(fawn.an.sur = 0.7, juv.an.sur = 0.9, ad.an.f.sur = 0.95, 
#' ad.an.m.sur = 0.8, fawn.repro = 0, juv.repro = 0.4, ad.repro = .9, 
#' hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.2,
#' hunt.mort.ad.f = 0.15, hunt.mort.ad.m = 0.35, ini.fawn.prev = 0.01,
#' ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
#' n.age.cats = 12,  p = 0.27, env.foi = 0,  beta.f = 0.08,  beta.m = 0.08,
#' theta = 1, n0 = 1000, n.years = 10, rel.risk = 1.0, 
#' repro.var = 0.005, fawn.sur.var = 0.005, sur.var = 0.005, hunt.var = 0.005)
#' 
#' out <- cwd_stoch_model(params)
#' 
#' plot_tots(out$counts)
#' 
#' @export


# Stochastic monthly age and sex structured model that incorporates random draws
# from distibutions of natural survival, reproduction, and hunt mortality.
# Currently does not include a distribution on transmission rate.
cwd_stoch_model <- function(params) {
  
  # write the list objects to the local environment
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  
  
  # check if parameters exist.    
  if(exists("fawn.an.sur")==FALSE){
    message("fawn survival is missing, using default value")
    fawn.an.sur <- 0.6
  }
  
  if(exists("juv.an.sur")==FALSE){
    message("juvenile survival is missing, using default value")
    juv.an.sur <- 0.8
  }
  
  if(exists("ad.an.f.sur")==FALSE){
    message("adult female survival is missing, using default value")
    ad.an.f.sur <- 0.95
  }
  
  if(exists("ad.an.m.sur")==FALSE){
    message("adult male survival is missing, using default value")
    ad.an.m.sur <- 0.9
  }
  
  if(exists("fawn.repro")==FALSE){
    message("fawn repro is missing, using default value")
    fawn.repro <- 0
  }
  
  if(exists("juv.repro")==FALSE){
    message("juvenile repro is missing, using default value")
    juv.repro <- 0.6
  }
  
  if(exists("ad.repro")==FALSE){
    message("adult repro is missing, using default value")
    ad.repro <- 1
  }
  
  if(exists("hunt.mort.fawn")==FALSE){
    message("fawn hunting mortality is missing, using default value")
    hunt.mort.fawn <- 0.01
  }
  
  if(exists("hunt.mort.juv.f")==FALSE){
    message("juv. female  hunting mortality is missing, using default value")
    hunt.mort.juv.f <- 0.1
  }
  if(exists("hunt.mort.juv.m")==FALSE){
    message("juv. male hunting mortality is missing, using default value")
    hunt.mort.juv.m <- 0.1
  }
  
  if(exists("hunt.mort.ad.f")==FALSE){
    message("adult female hunting mortality is missing, using default value")
    hunt.mort.ad.f <- 0.1
  }
  if(exists("hunt.mort.ad.m")==FALSE){
    message("adult male hunting mortality is missing, using default value")
    hunt.mort.ad.m <- 0.2
  }
  if(exists("ini.fawn.prev")==FALSE){
    message("initial fawn prevalence is missing, using default value")
    ini.fawn.prev <- 0.01
  }
  if(exists("ini.juv.prev")==FALSE){
    message("initial juvenile prevalence is missing, using default value")
    ini.juv.prev <- 0.03
  }
  if(exists("ini.ad.f.prev")==FALSE){
    message("initial adult female prevalence is missing, using default value")
    ini.ad.f.prev <- 0.04
  }
  if(exists("ini.ad.m.prev")==FALSE){
    message("initial adult male prevalence is missing, using default value")
    ini.ad.m.prev <- 0.04
  }
  if(exists("n.age.cats")==FALSE){
    message("# of age categories is missing, using default value")
    n.age.cats <- 12
  }
  if(exists("n.age.cats.m")==FALSE){
    message("# of male age categories is missing, using default value")
    n.age.cats.m <- 10
  }
  
  if(exists("n.age.cats.f")==FALSE){
    message("# of female age categories is missing, using default value")
    n.age.cats.f <- 15
  }
  
  if(exists("p")==FALSE){
    message("disease mortality index p is missing, using default value")
    p <- 0.43
  }
  if(exists("env.foi")==FALSE){
    message("indirect transmission env.foi is missing, using default value")
    env.foi <- 0
  }
  
  if(exists("beta.f")==FALSE){
    message("female transmission beta.f is missing, using default value")
    beta.f <- 0.08
  }
  
  if(exists("beta.m")==FALSE){
    message("male transmission beta.m is missing, using default value")
    beta.m <- 0.08
  }
  
  if(exists("theta")==FALSE){
    message("theta is missing, using default value")
    theta <- 1
  }
  
  if(exists("n0")==FALSE){
    message("initial population size n0 is missing, using default value")
    n0 <- 10000
  }
  
  if(exists("n.years")==FALSE){
    message("n.years is missing, using default value")
    n.years <- 10
  }
  
  if(exists("rel.risk")==FALSE){
    message("rel.risk is missing, using default value")
    rel.risk <- 1
  }
  
  if(exists("repro.var")==FALSE){
    message("repro.var is missing, using default value")
    repro.var <- 0.005
  }
  if(exists("fawn.sur.var")==FALSE){
    message("fawn.sur.var is missing, using default value")
    fawn.sur.var <- 0.005
  }
  if(exists("sur.var")==FALSE){
    message("sur.var is missing, using default value")
    sur.var <- 0.005
  }
  if(exists("hunt.var")==FALSE){
    message("hunt.var is missing, using default value")
    hunt.var <- 0.005
  }
  
  ### Adding additional variance parameters here.. ###
  ## survival
  if(exists("juv.sur.var")==FALSE){
    message("juv.sur.var is missing, using default value")
    juv.sur.var <- 0.005
  }
  if(exists("ad.f.sur.var")==FALSE){
    message("ad.f.sur.var is missing, using default value")
    ad.f.sur.var <- 0.005
  }
  if(exists("ad.m.sur.var")==FALSE){
    message("ad.m.sur.var is missing, using default value")
    ad.m.sur.var <- 0.005
  }
  ## reproduction
  if(exists("juv.repro.var")==FALSE){
    message("juv.repro.var is missing, using default value")
    juv.repro.var <- 0.005
  }
  if(exists("ad.repro.var")==FALSE){
    message("ad.repro.var is missing, using default value")
    ad.repro.var <- 0.005
  }
  #note - model doesn't have fawn reproduction. 
  
  if(exists("WSI")==FALSE){
    message("WSI is missing, using default value")
    WSI <- 0 # defaults to "off"
  }
  
  ### New action params
  
  if(exists("nosampled")==FALSE){
    message("nosampled is missing, using default value")
    nosampled <- 5 ## Status Quo value
  }
  
  if(exists("arrival_input")==FALSE){
    message("arrival_input is missing, using default value")
    arrival_input <- c(0,0,0,0,0,0,0, 0, 0, 0)
  }  
  ### Hunting action switches
  if(exists("Action_young_bucks")==FALSE){
    message("Action_young_bucks is missing, using default value")
    Action_young_bucks <- 0 #hunting actions are turned off
  }
  if(exists("Action_lib_harvest")==FALSE){
    message("Action_lib_harvest is missing, using default value")
    Action_lib_harvest <- 0 #hunting actions are turned off
  }

  if(exists("Action_sharpshooting")==FALSE){
    message("Action_sharpshooting is missing, using default value")
    Action_sharpshooting <- 0 #hunting actions are turned off
  }
  
  
  ###### check parameter values ###
  if(fawn.an.sur < 0) warning("fawn survival must be positive")
  if(fawn.an.sur > 1) warning("fawn survival must be < 1")
  if(juv.an.sur < 0) warning("juvenile survival must be positive")
  if(juv.an.sur > 1) warning("juvenile survival must be < 1")
  if(ad.an.f.sur < 0) warning("adult female survival must be positive")
  if(ad.an.f.sur > 1) warning("adult female survival must be < 1")
  
  if(fawn.repro < 0) warning("fawn.repro must be >= 0")
  if(juv.repro < 0) warning("juv.repro must be >= 0 ")
  if(ad.repro  < 0) warning("ad.repro must be >= 0 ")
  
  if(hunt.mort.fawn <= 0) warning("hunt.mort.fawn must be =0")
  if(hunt.mort.fawn >= 1) warning("hunt.mort.fawn must be < 1")
  if(hunt.mort.juv.f <= 0) warning("hunt.mort.juv.f must be >0")
  if(hunt.mort.juv.f >= 1) warning("hunt.mort.juv.f must be < 1")
  if(hunt.mort.juv.m <= 0) warning("hunt.mort.juv.m must be >0")
  if(hunt.mort.juv.m >= 1) warning("hunt.mort.juv.m must be < 1")
  if(hunt.mort.ad.f <= 0) warning("hunt.mort.ad.f must be >0")
  if(hunt.mort.ad.f >= 1) warning("hunt.mort.ad.f must be < 1")
  if(hunt.mort.ad.m <= 0) warning("hunt.mort.ad.m must be >0")
  if(hunt.mort.ad.m >= 1) warning("hunt.mort.ad.m must be < 1")
  
  if(ini.fawn.prev < 0) warning("ini.fawn.prev must >=0")
  if(ini.fawn.prev > 1) warning("ini.fawn.prev must be <= 1")
  if(ini.juv.prev < 0) warning("ini.juv.prev must >=0")
  if(ini.juv.prev > 1) warning("ini.juv.prev must be <= 1")
  if(ini.ad.f.prev < 0) warning("ini.ad.f.prev must >=0")
  if(ini.ad.f.prev > 1) warning("ini.ad.f.prev must be <= 1")
  if(ini.ad.m.prev < 0) warning("ini.ad.m.prev must >=0")
  if(ini.ad.m.prev > 1) warning("ini.ad.m.prev must be <= 1")
  
  if(n.age.cats < 3) warning("n.age.cats must be 3 or more")
  if(p < 0) warning("p must be between 0 and 1")
  if(p > 1) warning("p must be between 0 and 1")
  if(env.foi < 0) warning("env.foi must be between 0 and 1")
  if(env.foi > 1) warning("env.foi must be between 0 and 1")
  if(beta.f < 0) warning("beta.f cannot be negative")
  if(beta.m < 0) warning("beta.m cannot be negative")
  if(n0 <= 0) warning("n0 must be positive")
  if(n.years <= 0) warning("n.years must be positive")
  if(rel.risk <= 0) warning("n.years must be positive")
  
  if(repro.var <= 0) warning("repro.var must be positive")
  if(fawn.sur.var <= 0) warning("fawn.sur.var must be positive")
  if(sur.var <= 0) warning("sur.var must be positive")
  if(hunt.var <= 0) warning("hunt.var must be positive")
  
  ######### CREATE INITIAL CONDITIONS##########
  
  
  months <- seq(1, n.years * 12)  # monthly timestep
  hunt.mo <- rep(0, n.years * 12)  # months in where the hunt occurs
  hunt.mo[months%%12 == 7] <- 1  # hunt.mo==1 on Nov
  
  ### Set up WSI
  WSI.mo <- rep(0, n.years * 12)
  WSI.mo[months%%12 == 8 | months%%12 == 9 | months%%12 == 10 | months%%12 == 11] <- 1  # hunt.mo==1 on Dec:March
  
  # Estimate shape and scale parameters for the Beta distribution given the user
  # input of mean and variance.  natural survival
  fawn.s.b <- est_beta_params(fawn.an.sur, fawn.sur.var)
  juv.s.b <- est_beta_params(juv.an.sur, juv.sur.var)
  ad.f.s.b <- est_beta_params(ad.an.f.sur, ad.f.sur.var)
  ad.m.s.b <- est_beta_params(ad.an.m.sur, ad.m.sur.var)
  
  # reproduction
  juv.r.b <- est_beta_params(juv.repro/2, juv.repro.var)
  ad.r.b <- est_beta_params(ad.repro/2, ad.repro.var)
  
  # hunting
  hunt.fawn.b <- est_beta_params(hunt.mort.fawn, hunt.var)
  hunt.juv.f.b <- est_beta_params(hunt.mort.juv.f, hunt.var)
  hunt.juv.m.b <- est_beta_params(hunt.mort.juv.m, hunt.var)
  hunt.f.b <- est_beta_params(hunt.mort.ad.f, hunt.var)
  hunt.m.b <- est_beta_params(hunt.mort.ad.m, hunt.var)
  
  
  # Create the Leslie Matrix to start the population at stable age dist
  # initial female prevalence
  ini.f.prev <- c(ini.fawn.prev, ini.juv.prev,
                  rep(ini.ad.f.prev, (n.age.cats.f - 2)))
  # initial male prevalence
  ini.m.prev <- c(ini.fawn.prev, ini.juv.prev,
                  rep(ini.ad.m.prev, (n.age.cats.m - 2)))
  
  # Create the Leslie Matrix to start the population at stable age dist
  M <- matrix(rep(0, (n.age.cats.f + n.age.cats.m) * (n.age.cats.f + n.age.cats.m) ), nrow = (n.age.cats.f + n.age.cats.m))
  
  # replace the -1 off-diagonal with the survival rates
  M[row(M) == (col(M) + 1)] <- c(juv.an.sur * (1 - hunt.mort.juv.f),
                                 rep(ad.an.f.sur *  (1 - hunt.mort.ad.f),
                                     n.age.cats.f - 2), 0,
                                 c(juv.an.sur * (1 - hunt.mort.juv.m),
                                   rep(ad.an.m.sur * (1 - hunt.mort.ad.m),
                                       n.age.cats.m - 2)))
  
  # if you want the top age category to continue to survive
  M[n.age.cats.f, n.age.cats.f] <- ad.an.f.sur * (1 - hunt.mort.ad.f)
  M[(n.age.cats.f + n.age.cats.m) , (n.age.cats.f + n.age.cats.m)] <- ad.an.m.sur * (1 - hunt.mort.ad.m)
  
  # insert the fecundity vector prebirth census
  M[1, 1:n.age.cats.f] <- c(0, juv.repro, rep(ad.repro, n.age.cats.f - 2)) * 0.5 *
    fawn.an.sur * (1 - hunt.mort.fawn)
  M[n.age.cats.f + 1, 1:n.age.cats.m] <- M[1, 1:n.age.cats.m]
  
  # pre-allocate the output matrices
  tmp.f <- matrix(0, nrow = n.age.cats.f, ncol = n.years * 12)
  tmp.m <- matrix(0, nrow = n.age.cats.m, ncol = n.years * 12)
  St.f <- tmp.f  # susceptible female vector
  St.m <- tmp.m  # suceptible male vector
  # infectious categories
  It.m <- array(rep(tmp.m), dim = c(n.age.cats.m, n.years * 12, 10))  # females
  It.f <- array(rep(tmp.f), dim = c(n.age.cats.f, n.years * 12, 10))  # males
  
  # tracking the # hunted
  Ht.f <- tmp.f
  Ht.m <- tmp.m
  # natural deaths
  Dt.f <- tmp.f
  Dt.m <- tmp.m
  # disease deaths
  CWDt.f <- tmp.f
  CWDt.m <- tmp.m
  
  # tracking samples
  Samp.f <- tmp.f
  Samp.m <- tmp.m
  I.Samp.f <- tmp.f
  I.Samp.m <- tmp.m
  
  # Sharpshooters..
  Sharp.f <- tmp.f
  Sharp.m <- tmp.m
  
  arrival <- c()
  for(i in arrival_input){
    x <- c(i, rep(0, 11))
    arrival <- c(arrival, x)
  }

  # group into a vector
  # initial female prevalence
  ini.f.prev <- c(ini.fawn.prev, ini.juv.prev,
                  rep(ini.ad.f.prev, (n.age.cats.f - 2)))
  # initial male prevalence
  ini.m.prev <- c(ini.fawn.prev, ini.juv.prev,
                  rep(ini.ad.m.prev, (n.age.cats.m - 2)))
  
  ## These are based on arrival vector..
  if(arrival[1] == 0){ #even split between males and females 
    set.ini.f.prev <- rep(0, n.age.cats.f)
    set.ini.m.prev <- rep(0, n.age.cats.m)
    
    # randomly allocating infecteds across ages and categories.
    It.f[, 1, 1:10] <- rbinom(n.age.cats.f * 10, round(popbio::stable.stage(M)[1:n.age.cats.f] *
                                                         n0/10), set.ini.f.prev)
    It.m[, 1, 1:10] <- rbinom(n.age.cats.m * 10, round(popbio::stable.stage(M)
                                                       [(n.age.cats.f + 1):(n.age.cats.f + n.age.cats.m)] *
                                                         n0/10), set.ini.m.prev)
    
    # Intializing with the stable age distribution.
    St.f[, 1] <- round(popbio::stable.stage(M)[1:n.age.cats.f] * n0 * (1 - set.ini.f.prev))
    St.m[, 1] <- round(popbio::stable.stage(M)[(n.age.cats.f + 1):(n.age.cats.f + n.age.cats.m)] *
                         n0 * (1 - set.ini.m.prev))
    
  } else{ #even split between males and females 
    split <- as.vector(table(sample(1:2, size = arrival[1], replace = T)))
    if(length(split) == 1){ split[2] <- 0} # otherwise have issues
    
    new.inf.f <- as.vector(rmultinom(1, size = split[1], prob = ini.f.prev))
    new.inf.m <- as.vector(rmultinom(1, size = split[2], prob = ini.m.prev))
    
    ### randomly allocating infecteds across ages and categories.
    random_allocation.f <- matrix(0, nrow = length(new.inf.f), ncol = 10)
    for(j in 1:nrow(as.matrix(new.inf.f))){
      col_index <- sample(1:10, 1)  # randomly choose column index
      random_allocation.f[j, col_index] <- random_allocation.f[j, col_index] + new.inf.f[j]  # add the value to the chosen column
    }
    
    random_allocation.m <- matrix(0, nrow = length(new.inf.m), ncol = 10)
    for(j in 1:nrow(as.matrix(new.inf.m))){
      col_index <- sample(1:10, 1)  # randomly choose column index
      random_allocation.m[j, col_index] <- random_allocation.m[j, col_index] + new.inf.m[j]  # add the value to the chosen column
    }
    
    It.f[, 1, 1:10] <- random_allocation.f
    It.m[, 1, 1:10] <- random_allocation.m
    
    ### None effected popn
    
    set.ini.f.prev <- new.inf.f/n0
    set.ini.m.prev <- new.inf.m/n0
    
    # Intializing with the stable age distribution.
    St.f[, 1] <- round(popbio::stable.stage(M)[1:n.age.cats.f] * n0 * (1 - set.ini.f.prev))
    St.m[, 1] <- round(popbio::stable.stage(M)[(n.age.cats.f + 1):(n.age.cats.f + n.age.cats.m)] *
                         n0 * (1 - set.ini.m.prev))
    
  } 
  
  
  if(sum(St.f[,1]) <= 0) {
    warning("These parameters result in a stable age structure with no surviving 
            females.")
  } 
  
  # calculate R0s for adult females and males
  # in the denominator find the average minimum survival for the 3 mortality types
  f.R0 <-  (beta.f * n0) / (n0 ^ theta) * 
    mean(apply(cbind(rnbinom(1000, 1, (1 - ad.an.f.sur^(1/12))), 
                     rnbinom(1000, 1, (1 - (1 - hunt.mort.ad.f)^(1/12))),
                     rgamma(1000, 10, p)), 1, FUN = min, na.rm = T))
  
  m.R0 <-  (beta.m * n0)  / (n0 ^ theta) *
    mean(apply(cbind(rnbinom(1000, 1, (1 - ad.an.m.sur^(1/12))), 
                     rnbinom(1000, 1, (1 - (1 - hunt.mort.ad.m)^(1/12))),
                     rgamma(1000, 10, p)), 1, FUN = min, na.rm = T))
  
  rm(M)
  
  
  ####### POPULATION MODEL############
  for (t in 2:(n.years * 12)) {
    
    # Annual Parameter Draws monthly stochastic survival rates
    fawn.sur.draw <- rbeta(1, fawn.s.b$alpha, fawn.s.b$beta, ncp = 0)^(1/12)
    juv.sur.draw <- rbeta(1, juv.s.b$alpha, juv.s.b$beta, ncp = 0)^(1/12)
    ad.f.sur.draw <- rbeta(1, ad.f.s.b$alpha, ad.f.s.b$beta, ncp = 0)^(1/12)
    ad.m.sur.draw <- rbeta(1, ad.m.s.b$alpha, ad.m.s.b$beta, ncp = 0)^(1/12)
    
    # monthly stochastic reproductive rates
    juv.preg.draw <- rbeta(1, juv.r.b$alpha, juv.r.b$beta, ncp = 0)
    ad.preg.draw <- rbeta(1, ad.r.b$alpha, ad.r.b$beta, ncp = 0)
    
    
    ## See if WSI is happening
    if (WSI == 1){ ## Check to see if WSI is turned on in params
      if (t%%12 == 2) { # In december, draw from binomial to see if winter will be sever
        WSI_draw <- rbinom(1, size = 1, prob = 0.125)
      }
    } else{
      WSI_draw = 0
    }
    
    ## If it is, redo survival parameters.. 
    if(WSI_draw == 1){
      if (WSI.mo[t] == 1) {
        
        fawn.sur.draw <- rbeta(1, fawn.s.b$alpha, fawn.s.b$beta, ncp = 0)^(1/12)*.83
        juv.sur.draw <- rbeta(1, juv.s.b$alpha, juv.s.b$beta, ncp = 0)^(1/12)*.87
        ad.f.sur.draw <- rbeta(1, ad.f.s.b$alpha, ad.f.s.b$beta, ncp = 0)^(1/12)*.96
        ad.m.sur.draw <- rbeta(1, ad.m.s.b$alpha, ad.m.s.b$beta, ncp = 0)^(1/12)*.96
        
      } else{
        
        fawn.sur.draw <- rbeta(1, fawn.s.b$alpha, fawn.s.b$beta, ncp = 0)^(1/12)
        juv.sur.draw <- rbeta(1, juv.s.b$alpha, juv.s.b$beta, ncp = 0)^(1/12)
        ad.f.sur.draw <- rbeta(1, ad.f.s.b$alpha, ad.f.s.b$beta, ncp = 0)^(1/12)
        ad.m.sur.draw <- rbeta(1, ad.m.s.b$alpha, ad.m.s.b$beta, ncp = 0)^(1/12)
      }
    }
    
    
    # group into a vector
    Sur.f <- c(fawn.sur.draw, juv.sur.draw, rep(ad.f.sur.draw, n.age.cats.f - 2))
    Sur.m <- c(fawn.sur.draw, juv.sur.draw, rep(ad.m.sur.draw, n.age.cats.m - 2))
    
    # stochastic hunting survival rates
    hunt.fawn.draw <- rbeta(1, hunt.fawn.b$alpha, hunt.fawn.b$beta, ncp = 0)
    hunt.juv.f.draw <- rbeta(1, hunt.juv.f.b$alpha, hunt.juv.f.b$beta, ncp = 0)
    hunt.juv.m.draw <- rbeta(1, hunt.juv.m.b$alpha, hunt.juv.m.b$beta, ncp = 0)
    hunt.f.draw <- rbeta(n.age.cats.f - 2, hunt.f.b$alpha, hunt.f.b$beta, ncp = 0)
    hunt.m.draw <- rbeta(n.age.cats.m - 2, hunt.m.b$alpha, hunt.m.b$beta, ncp = 0)
    
    # on birthdays add in recruits and age everyone by one year also on birthdays do
    # the random parameter draws births happen in June, model starts in May
    if (t%%12 == 2) {
      
      # the last age category remains in place and doesn't die 
      St.f[2:(n.age.cats.f - 1), t] <- St.f[1:(n.age.cats.f - 2), t - 1]
      St.f[n.age.cats.f, t] <- round((St.f[n.age.cats.f, t - 1] +
                                        St.f[(n.age.cats.f - 1), t - 1]) )
      St.m[2:(n.age.cats.m - 1), t] <- St.m[1:(n.age.cats.m - 2), t - 1]
      St.m[n.age.cats.m, t] <- round((St.m[n.age.cats.m, t - 1] +
                                        St.m[(n.age.cats.m - 1), t - 1]))
      It.f[2:(n.age.cats.f - 1), t, ] <- It.f[1:(n.age.cats.f - 2), t - 1, ]
      It.f[n.age.cats.f, t, ] <- round((It.f[n.age.cats.f, t - 1, ] +
                                          It.f[(n.age.cats.f - 1), t - 1, ]))
      It.m[2:(n.age.cats.m - 1), t, ] <- It.m[1:(n.age.cats.m - 2), t - 1, ]
      It.m[n.age.cats.m, t, ] <- round((It.m[n.age.cats.m, t - 1, ] +
                                          It.m[(n.age.cats.m - 1), t - 1, ]))
      
      # reproduction
      I_juv <- sum(It.f[2, t - 1, ])
      I_adults <- sum(It.f[3:n.age.cats.f, t - 1, ])
      
      fawns_born <- rbinom(1, (St.f[2, t - 1] + I_juv), juv.preg.draw) *2  +
        rbinom(1, (sum(St.f[3:12, t - 1]) + I_adults), ad.preg.draw) *2 
      
      St.f[1, t] <- rbinom(1, fawns_born, 0.5)
      St.m[1, t] <- fawns_born - St.f[1, t]
    }
    
    # if not June populate the current month with last months values
    if (t%%12 != 2) {
      # updating the next month
      St.f[, t] <- St.f[, t - 1]
      St.m[, t] <- St.m[, t - 1]
      It.f[, t, ] <- It.f[, t - 1, ]
      It.m[, t, ] <- It.m[, t - 1, ]
    }
    
    ## Natural Mort then hunt then disease mort Then transmission Natural Mortality
    ## susceptibles
    nat.s.f <- rbinom(n.age.cats.f, St.f[, t], (1 - Sur.f))
    nat.s.m <- rbinom(n.age.cats.m, St.m[, t], (1 - Sur.m))
    
    St.f[, t] <- St.f[, t] - nat.s.f
    St.m[, t] <- St.m[, t] - nat.s.m
    # infecteds
    nat.i.f <- matrix(rbinom(length(It.f[, t, ]), size = It.f[, t, ],
                             prob = (1 - Sur.f)), nrow = n.age.cats.f)
    nat.i.m <- matrix(rbinom(length(It.m[, t, ]), size = It.m[, t, ],
                             prob = (1 - Sur.m)), nrow = n.age.cats.m)
    
    It.f[, t, ] <- It.f[, t, ] - nat.i.f
    It.m[, t, ] <- It.m[, t, ] - nat.i.m
    
    Dt.f[, t] <- nat.s.f + rowSums(nat.i.f)
    Dt.m[, t] <- nat.s.m + rowSums(nat.i.m)
    
    # Hunt mortality
    if (hunt.mo[t] == 1) {
      Iall.f <- rowSums(It.f[, t, ])  # total # infected females
      Iall.m <- rowSums(It.m[, t, ])  # total # infected males
      Nt.f <- St.f[, t] + Iall.f  # total population of females
      Nt.m <- St.m[, t] + Iall.m  # total population of males
      
      # binomial draw on the total hunted
      Ht.f[, t] <- rbinom(n.age.cats.f, Nt.f, c(hunt.fawn.draw, hunt.juv.f.draw,
                                                hunt.f.draw))
      
      Ht.m[, t] <- rbinom(n.age.cats.m, Nt.m, c(hunt.fawn.draw, hunt.juv.m.draw,
                                                hunt.m.draw))
      
      
      ###If harvesting too many fawns, redraw harvest numbers
      if(!is.na(Ht.m[1, t]) & !is.na(Ht.f[1, t])){ # This was bugging out
        if (Ht.f[1, t] > 5000 | Ht.m[1, t] > 5000){ #If harvesting too many fawns, redraw harvest numbers
          Ht.f[, t] <- rbinom(n.age.cats.f, Nt.f, c(hunt.fawn.draw, hunt.juv.f.draw,
                                                    hunt.f.draw))
          Ht.m[, t] <- rbinom(n.age.cats.m, Nt.m, c(hunt.fawn.draw, hunt.juv.m.draw,
                                                    hunt.m.draw))
          if (Ht.f[1, t] > 5000 | Ht.m[1, t] > 5000){
            Ht.f[, t] <- rbinom(n.age.cats.f, Nt.f, c(hunt.fawn.draw, hunt.juv.f.draw,
                                                      hunt.f.draw))
            Ht.m[, t] <- rbinom(n.age.cats.m, Nt.m, c(hunt.fawn.draw, hunt.juv.m.draw,
                                                      hunt.m.draw))
          } 
        } 
      }
      
      
      ### Switch for infect parameter (need to be able to turn this on, not hold at 0)
      if(!exists("infect")){
        infect = 0 # holder until switched with detection
      }
      
      ### Actions
      if(infect != 0){
        if(Action_young_bucks == 1){
          ### Target yearling bucks to reduce density demography [8 - 10%]
          # Overwriting the season's harvest after bump harvest numbers
          Ht.m[2, t] <- round(Ht.m[2, t]*1.1)
        }
        if(Action_lib_harvest == 1){ # produces warnings bc different lengths but functionally fine
          ### Liberalize harvest season 
          totalpopn <- sum(Nt.f + Nt.m)
          totalhunted <- sum(Ht.f[, t] + Ht.m[, t])
          huntingtarget <- (15000/143758) * n0 # make 4x bigger 
          prop.popn.to.take <- ((totalhunted + huntingtarget)/totalpopn )
          
          # want to pull a number btw 10 - 20% 
          huntertake <- sample(10:20, 1)*.01
          if(prop.popn.to.take >  huntertake){ # check to make sure not harvesting more of popn than random draw of hunter ability
            if(totalhunted/totalpopn > huntertake){
              huntingtarget <- 2 # effectively not increasing harvest 
            }else{
              #If harvesting more than draw, modify take to hit hunter draw
              modify.take <- ((huntertake * totalpopn) - totalhunted)*(1/huntingtarget)
              huntingtarget <- huntingtarget *modify.take
            }
          }
          
          #get prob for splitting up antlerless
          antlerless.probs <- c(Nt.m[1]/(Nt.m[1] + sum(Nt.f)), sum(Nt.f)/(Nt.m[1] + sum(Nt.f)))
          
          
          split <- as.vector(table(sample(1:2, size = huntingtarget, replace = T, prob = antlerless.probs)))
          
          if(huntingtarget == 2){ ## making sure no weird errors happen
            split <- c(1,1)
          }    
          if(length(split) == 0){ ## another weird error
            split <- c(1,1)
          }
          if(length(split) == 1){ ## another weird error
            split[2] <- c(1)
          }
          
          proportions.f <- Nt.f / sum(Nt.f)  # Calculate proportion
          Ht.f[, t] <- Ht.f[, t] + round(proportions.f * split[2])  # Calculate harvested numbers
            if(split[2] > Nt.m[1]){ #stops overharvesting male fawns
              split[2] <- Nt.m[1]
            }
          Ht.m[1, t] <- Ht.m[1, t] + split[1] # change from Ht.f[1, t] + split[2]
          
        }
      } # close first action loop 
      
      
      
      
      # those hunted in the I class overall based on the total hunted, the total that
      # are susceptible/infected and the relative hunting risk of S v. I can result in
      # a divide by 0 and NA.  this can also result in more hunting of a category than
      # are available.
      
      hunted.i.f <- round((rel.risk * Iall.f * Ht.f[, t]) /
                            (St.f[, t] + rel.risk * Iall.f))
      hunted.i.m <- round((rel.risk * Iall.m * Ht.m[, t]) /
                            (St.m[, t] + rel.risk * Iall.m))
      
      hunted.i.f[which(is.na(hunted.i.f))] <- 0
      hunted.i.m[which(is.na(hunted.i.m))] <- 0
      
      hunted.i.f[Iall.f < hunted.i.f] <- Iall.f[Iall.f < hunted.i.f]
      hunted.i.m[Iall.m < hunted.i.m] <- Iall.m[Iall.m < hunted.i.m]
      
      # subtracting out those hunted in the S class
      St.f[, t] <- St.f[, t] - (Ht.f[, t] - hunted.i.f)
      St.m[, t] <- St.m[, t] - (Ht.m[, t] - hunted.i.m)
      
      # allocate those deaths across the 10 I categories
      It.f[, t, ] <- allocate_deaths(hunted.i.f, It.f[, t, ])
      It.m[, t, ] <- allocate_deaths(hunted.i.m, It.m[, t, ])
    
      
      
      ### Surveillance starts here.. 
      # bio check 10% or less of total harvest - randomly select 10% of available harvest
      # no. individuals available for sampling
      samp <- (sum(Ht.f[, t]) + sum(Ht.m[, t]))*.1
      # split males/females 50/50?
      
      ### Check if number of individuals available for sampling is more or less than target
      # target = nosampled (set in params)
      if(samp < nosampled){
        availableSamples <- samp
      } else{
        availableSamples <- nosampled
      }
      
      samp.f <- round(availableSamples*.5)
      samp.m <- round(availableSamples*.5)
      
      # randomly draw samples
      Samp.f[,t ] <- randomsampling(matrix = Ht.f[, t], target_sum = samp.f, emptysample = Samp.f[,t ])
      Samp.m[,t ] <- randomsampling(matrix = Ht.m[, t], target_sum = samp.m, emptysample = Samp.m[,t ])
      
      # see how many of the samples have infected individuals..
      time <- t
      I.Samp.f <- infectedSamples(sample = Samp.f[,t ], infected = hunted.i.f, 
                                  removed = Ht.f[, t], trackInfect = I.Samp.f, time = time)
      I.Samp.m <- infectedSamples(sample = Samp.m[,t ], infected = hunted.i.m, 
                                  removed = Ht.m[, t], trackInfect = I.Samp.m, time = time)
      
      ## Check to see if found any infected individuals
      infect <- sum(I.Samp.f) + sum(I.Samp.m)
      
    
      
      ## !!! IF find infected individual then 
      
      if(infect != 0 & !is.na(infect)){

        if(Action_sharpshooting == 1){
          
          ### Change rel.risk parameter only for sharpshooting
          rel.risk <- 1.5 

          
            ### Sharp shooting  
            # Pull 300 per year until population drops to 5 deer/square mile 
            # In the simulation that would be ?? 
            # According to spreadsheet, target is approximately 15/square mile
            # Assuming we have that with calibration
            # would stop pulling when reduce population by 2/3s 
            
            ### Set up assumes that have more females than males.. 
            stopvalue <- (sum(It.f[1:n.age.cats.m, 2, ] + It.m[, 2, ]) + sum(It.f[(n.age.cats.m+1):n.age.cats.f, 2, ]) + 
                            sum(St.f[1:n.age.cats.m, 2] + St.m[, 2]  + St.f[(n.age.cats.m+1):n.age.cats.f, 2]))*1/3
            #stopvalue = starting popn * 1/3 (note assume don't have a burn in, if so change from 2)
            
            currentvalue <- (sum(It.f[1:n.age.cats.m, t, ] + It.m[, t, ]) + sum(It.f[(n.age.cats.m+1):n.age.cats.f, t, ]) + 
                               sum(St.f[1:n.age.cats.m, t] + St.m[, t]  + St.f[(n.age.cats.m+1):n.age.cats.f, t]))
            
            if(currentvalue > stopvalue){
              checksize = currentvalue - stopvalue
              if(checksize > 300){
                setsize = 300
              } else{
                setsize = checksize
              }
              ## Need to pull 300 individuals but this should be random?
              ## Can repurpose code used to pull hunting
              
              #What should the split between males and females be? 
              split <- as.vector(table(sample(1:2, size = setsize, replace = T)))
              
              # ## Redo this, right now, harvesting negative individuals... 
              # Sharp.f[, t] <- rbinom(n.age.cats, split[1], c(rep(1/n.age.cats, n.age.cats)))
              # Sharp.m[, t] <- rbinom(n.age.cats, split[2], c(rep(1/n.age.cats, n.age.cats)))
              
              Iall.f <- rowSums(It.f[, t, ])  # total # infected females
              Iall.m <- rowSums(It.m[, t, ])  # total # infected males
              Nt.f <- St.f[, t] + Iall.f  # total population of females
              Nt.m <- St.m[, t] + Iall.m  # total population of males
              
              proportions.f <- Nt.f / sum(Nt.f)  # Calculate proportions
              proportions.m <- Nt.m / sum(Nt.m) 
              Sharp.f[, t] <- round(proportions.f * split[1])  # Calculate harvested numbers
              Sharp.m[, t] <- round(proportions.m * split[2])  
              ### Now sharpshooters are harvesting based on what ages in popn are most available
              
              
              # those hunted in the I class overall based on the total hunted, the total that
              # are susceptible/infected and the relative hunting risk of S v. I can result in
              # a divide by 0 and NA.  this can also result in more hunting of a category than
              # are available.
              
              Sharped.i.f <- round((rel.risk * Iall.f * Sharp.f[, t]) /
                                     (St.f[, t] + rel.risk * Iall.f))
              Sharped.i.m <- round((rel.risk * Iall.m * Sharp.m[, t]) /
                                     (St.m[, t] + rel.risk * Iall.m))
              
              Sharped.i.f[which(is.na(Sharped.i.f))] <- 0
              Sharped.i.m[which(is.na(Sharped.i.m))] <- 0
              
              ### Can't remove negative individuals
              Sharped.i.f[which(Sharped.i.f<0)] <- 0
              Sharped.i.m[which(Sharped.i.m<0)] <- 0
              
              Sharped.i.f[Iall.f < Sharped.i.f] <- Iall.f[Iall.f < Sharped.i.f]
              Sharped.i.m[Iall.m < Sharped.i.m] <- Iall.m[Iall.m < Sharped.i.m]
              
              # subtracting out those removed by sharpshooters in the S class
              St.f[, t] <- St.f[, t] - (Sharp.f[, t] - Sharped.i.f)
              St.m[, t] <- St.m[, t] - (Sharp.m[, t] - Sharped.i.m)
              
              # allocate those deaths across the 10 I categories
              It.f[, t, ] <- allocate_deaths(Sharped.i.f, It.f[, t, ])
              It.m[, t, ] <- allocate_deaths(Sharped.i.m, It.m[, t, ])
              
              ### Change rel.risk parameter back to normal
              rel.risk <- 1
              
            }
          } 
          
          ### Breakup winter aggregation with harassment = reduce transmission by 25%
          ## **** Note this is only 2 strategies - include? Makes life harder..
          # beta.f = beta.f*.75 ; beta.m = beta.m*.75
        }
      } # close second action loop
    
    # Disease mortality stochastic movement of individuals from I1 to I2 disease
    # induced mortality here by advancing all I's and only a proportion of the 10th
    # category remains
    I.f.move <- matrix(rbinom(n.age.cats.f * 10, size = It.f[, t, ], prob = p),
                       nrow = n.age.cats.f)
    I.m.move <- matrix(rbinom(n.age.cats.m * 10, size = It.m[, t, ], prob = p),
                       nrow = n.age.cats.m)
    
    # store info on those that die directly from disease
    CWDt.f[, t] <- I.f.move[, 10]
    CWDt.m[, t] <- I.m.move[, 10]
    
    # move the I individuals forward in their categories
    It.f[, t, 1] <- It.f[, t, 1] - I.f.move[, 1]
    It.f[, t, 2:10] <- It.f[, t, 2:10] - I.f.move[, 2:10] + I.f.move[, 1:9]
    
    It.m[, t, 1] <- It.m[, t, 1] - I.m.move[, 1]
    It.m[, t, 2:10] <- It.m[, t, 2:10] - I.m.move[, 2:10] + I.m.move[, 1:9]
    
    # Direct transmission considering all I's are equal
    # Iall <- sum(It.f[, t, ] + It.m[, t, ])
    # Nall <- sum(St.f[, t] + St.m[, t]) + Iall
    
    ### Need to subset to add and then recombine? Might have been a better way to do this 
    # Making indexing dynamic for varying age class senarios 
    # Note, females must be greater than or equal to number of male age classes
    if(n.age.cats.m == n.age.cats.f){
      Iall <- sum(It.f[, t, ] + It.m[, t, ])
      Nall <- sum(St.f[, t] + St.m[, t]) + Iall
    } else if(n.age.cats.f - n.age.cats.m == 1){
      Iall <- sum(It.f[1:n.age.cats.m, t, ] + It.m[, t, ]) + sum(It.f[(n.age.cats.m+1), t, ])
      Nall <- sum(St.f[1:n.age.cats.m, t] + St.m[, t]) +  sum(St.f[(n.age.cats.m+1), t])  + Iall
    } else if(n.age.cats.f - n.age.cats.m > 1){
      Iall <- sum(It.f[1:n.age.cats.m, t, ] + It.m[, t, ]) + sum(It.f[(n.age.cats.m+1):n.age.cats.f, t, ])
      Nall <- sum(St.f[1:n.age.cats.m, t] + St.m[, t]) +  sum(St.f[(n.age.cats.m+1):n.age.cats.f, t])  + Iall
    }
    
    foi.f <- 1 - exp(-beta.f * Iall/Nall^theta)
    foi.m <- 1 - exp(-beta.m * Iall/Nall^theta)
    
    transmission.f <- rbinom(n.age.cats.f, St.f[, t], foi.f)
    transmission.m <- rbinom(n.age.cats.m, St.m[, t], foi.m)
    
    St.f[, t] <- St.f[, t] - transmission.f
    St.m[, t] <- St.m[, t] - transmission.m
    
    # update with the new infections
    # from transmission
    It.f[, t, 1] <- transmission.f + It.f[, t, 1]
    It.m[, t, 1] <- transmission.m + It.m[, t, 1] 
    
    ### from individuals arriving from out of state
    split <- as.vector(table(sample(1:2, size = arrival[t], replace = T)))
    ### Need to make this conditional or get error when nothing new to add
    if(sum(split) != 0){
      if(length(split) == 1){ split[2] <- 0} # otherwise have issues
      
      # ini.f.prev <- c(ini.fawn.prev, ini.juv.prev,
      #                 rep(ini.ad.f.prev, (n.age.cats - 2)))
      # # initial male prevalence
      # ini.m.prev <- c(ini.fawn.prev, ini.juv.prev,
      #                 rep(ini.ad.m.prev, (n.age.cats - 2)))
      
      ## Need to divide each split number into ago categories
      new.inf.f <- as.vector(rmultinom(1, size = split[1], prob = ini.f.prev))
      new.inf.m <- as.vector(rmultinom(1, size = split[2], prob = ini.m.prev))
      It.f[, t, 1] <- new.inf.f
      It.m[, t, 1] <- new.inf.m
    }
    # Environmental transmission happens last
    envcases.f <- rbinom(n.age.cats.f, St.f[, t], env.foi)
    envcases.m <- rbinom(n.age.cats.m, St.m[, t], env.foi)
    
    St.f[, t] <- St.f[, t] - envcases.f
    St.m[, t] <- St.m[, t] - envcases.m
    
    It.f[, t, 1] <- It.f[, t, 1] + envcases.f
    It.m[, t, 1] <- It.m[, t, 1] + envcases.m
  }
  # group the output
  counts <- list(St.f = St.f, St.m = St.m, I1t.f = It.f[, , 1],
                 I1t.m = It.m[, , 1], I2t.f = It.f[, , 2], I2t.m = It.m[, , 2],
                 I3t.f = It.f[, , 3], I3t.m = It.m[, , 3], I4t.f = It.f[, , 4],
                 I4t.m = It.m[, , 4], I5t.f = It.f[, , 5], I5t.m = It.m[ , , 5],
                 I6t.f = It.f[, , 6], I6t.m = It.m[, , 6], I7t.f = It.f[, , 7],
                 I7t.m = It.m[, , 7], I8t.f = It.f[, , 8], I8t.m = It.m[, , 8],
                 I9t.f = It.f[, , 9], I9t.m = It.m[, , 9], I10t.f = It.f[, , 10],
                 I10t.m = It.m[, , 10])
  
  deaths <- list(Ht.f = Ht.f, Ht.m = Ht.m, Dt.f = Dt.f, Dt.m = Dt.m, CWDt.f = CWDt.f,
                 CWDt.m = CWDt.m)
  
  mysurveillance <- list(I.Samp.f = I.Samp.f, I.Samp.m = I.Samp.m)
  
  sharpshooting <- list(Sharp.f = Sharp.f, Sharp.m = Sharp.m)
  
  
  # convert the output to long form
  counts.long <- melt(counts) %>%
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, - 1)), disease = "no")
  counts.long$disease[str_sub(counts.long$category, 1, 1) == "I"] <- "yes"
  counts.long$disease <- as.factor(counts.long$disease)
  
  deaths.long <- melt(deaths) %>%
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, - 1)))
  
  mysurveillance.long <- melt(mysurveillance) %>%
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, - 1)))
  
  sharpshooting.long <- melt(sharpshooting) %>%
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, - 1)))
  
  
  output <- list(counts = counts.long, deaths = deaths.long, survillance = mysurveillance.long,
                 sharpshooting = sharpshooting.long, f.R0 = f.R0, m.R0 = m.R0)
}

### Comment this all out when done
## Making graphs to check each action for single run of model

### Pick up here - make graph function for single run of model (use outputs..)
# arrival_line <- min(which(arrival > 0))/12


 # noact_prev <- plot_prev(output)
 # noact_abund <- plot_abundnace(output)
 # noact_harv <- plot_stoch_harvest(output, harvesttype= 1)
 # noact_harv_juv <- plot_stoch_harvest(output, harvesttype= 2)
 # noact_harv_antlerless <- plot_stoch_harvest(output, harvesttype= 3)

# prev <- plot_prev(output)
# abund <- plot_abundnace(output)
# harv <- plot_stoch_harvest(output, harvesttype= 1)
# harv_juv <- plot_stoch_harvest(output, harvesttype= 2)
# harv_antlerless <- plot_stoch_harvest(output, harvesttype= 3)
# 
# plot_sharpshooting(output)

