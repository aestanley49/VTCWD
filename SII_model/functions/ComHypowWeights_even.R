### The basis for this is that Hypos are all equally weighted, don't need to worry about weights and SD later

# source("CWD_model_Cross/functions/est_beta_params.R")
# source("CWD_model_Cross/functions/allocate_deaths.R")
# source("CWD_model_Cross/functions/MyFunctions.R")
# source("CWD_model_Cross/cwd_stoch_model_cal_act.r")
# source("CWD_model_Cross/ArrivalVectors.R")
# source("CWD_model_Cross/functions/cwd_stoch_wrapper_mod.r")
source("./SII_model/functions/CreateStableAgeStructure_func.R")
source("./SII_model/functions/popnmodel_func.R")
source("./SII_model/functions/sims_wrapper_func.R")
source("./SII_model/functions/randomsampling_func.R")

ComHypowWeights_even <- function(selectstrat = "SQ"){
  
  ### Set the parameters based on the selected straetgy 
  if(selectstrat == "SQ" | selectstrat == "SA" | selectstrat == "Ho"){
    params <- list(fawn.an.sur = 0.7, fawn.an.sur.var= 0.005, juv.f.an.sur = 0.8, juv.f.an.sur.var=0.005, juv.m.an.sur = .8, juv.m.an.sur.var=0.005, ad.an.f.sur = 0.87, ad.an.f.sur.var=0.005, ad.an.m.sur = 0.75, ad.an.m.sur.var=0.005,
                   
                   juv.repro = 2.6, ad.repro = 2.8,   repro.var = 0.005,
                   
                   hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1, hunt.mort.ad.f = 0.2, hunt.mort.ad.m = 0.45, hunt.var = 0.005,
                   
                   # turn stochasticity of vital rates on and off.. 
                   switch_stoch_survival_rates = 1, # (on = 1)
                   switch_stoch_hunt_rates = 1, # (on = 1)
                   switch_stoch_repro_rates = 1, # (on = 1)
                   
                   # temp set these to 0 for calibration
                   switchdiseasedy_on = 1, # turn disease dynamics on or off (on = 1)
                   arrival_input = c(0,0,0,0,0,0,0, 0, 0, 0), # currently set for 10 years 
                   
                   i.rate = .5,
                   env = 0, shedrate_I = .5, shedrate_H = .5, shedrate_ND = .5, expdecayconstant = .55, beta = .4,
                   
                   #model 1 (only one option for environmental transmission right now, turn on [1] or off [0])
                   switch_envirtrans = 0, # this selects set up for linear relationship between env and Fe
                   gamma = .01,
                   
                   
                   nyears = 200, n0 = 143000, K = 250000,
                   
                   ### Action parameters
                   nosampled = 5,
                   Action_sharpshooting = 0, 
                   Action_lib_harvest = 0, 
                   Action_young_bucks = 0,
                   
                   
                   ### Farm Years - comodel 
                   alpha = .01, beta1 = 2, # prevalence 
                   beta2 = -.8, # double fencing reduces by 80% 
                   beta3 = .1, # movement 
                   ## turn actions on or off
                   action1 = 1, # double fencing (indicator)
                   action2 = 1 # movement or export of captive individuals (action: restriction)
                   
    )
    if(selectstrat == "Ho"){
      params$nosampled <- 230 ## Need to overwrite number of samples in surviellance strategy tied to Ho
    }
  }else{ # other strategies have hunting actions turned on
    params <- list(fawn.an.sur = 0.7, fawn.an.sur.var= 0.005, juv.f.an.sur = 0.8, juv.f.an.sur.var=0.005, juv.m.an.sur = .8, juv.m.an.sur.var=0.005, ad.an.f.sur = 0.87, ad.an.f.sur.var=0.005, ad.an.m.sur = 0.75, ad.an.m.sur.var=0.005,
                   
                   juv.repro = 2.6, ad.repro = 2.8,   repro.var = 0.005,
                   
                   hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1, hunt.mort.ad.f = 0.2, hunt.mort.ad.m = 0.45, hunt.var = 0.005,
                   
                   # turn stochasticity of vital rates on and off.. 
                   switch_stoch_survival_rates = 1, # (on = 1)
                   switch_stoch_hunt_rates = 1, # (on = 1)
                   switch_stoch_repro_rates = 1, # (on = 1)
                   
                   # temp set these to 0 for calibration
                   switchdiseasedy_on = 1, # turn disease dynamics on or off (on = 1)
                   arrival_input = c(0,0,0,0,0,0,0, 0, 0, 0), # currently set for 10 years 
                   
                   i.rate = .5,
                   env = 0, shedrate_I = .5, shedrate_H = .5, shedrate_ND = .5, expdecayconstant = .55, beta = .4,
                   
                   #model 1 (only one option for environmental transmission right now, turn on [1] or off [0])
                   switch_envirtrans = 0, # this selects set up for linear relationship between env and Fe
                   gamma = .01,
                   
                   
                   nyears = 200, n0 = 143000, K = 250000,
                   
                   ### Action parameters
                   nosampled = 125,
                   Action_sharpshooting = 1, 
                   Action_lib_harvest = 1, 
                   Action_young_bucks = 1,
                   
                   
                   ### Farm Years - comodel 
                   alpha = .01, beta1 = 2, # prevalence 
                   beta2 = -.8, # double fencing reduces by 80% 
                   beta3 = .1, # movement 
                   ## turn actions on or off
                   action1 = 1, # double fencing (indicator)
                   action2 = 1 # movement or export of captive individuals (action: restriction)
    )
    if(selectstrat == "PareR" | selectstrat == "SK"){
      params$nosampled <- 460 ## Need to overwrite number of samples in surviellance strategy tied to PareR and Kitchen Sink
    }
  }
  
  
  simsout1 <- cwd_stoch_wrapper_arrvAVG(params, nsims = 50, n.years = 50, strat = selectstrat, hypothesis = "H1")
  simsout2 <- cwd_stoch_wrapper_arrvAVG(params, nsims = 50, n.years = 50, strat = selectstrat, hypothesis = "H2")
  simsout3 <- cwd_stoch_wrapper_arrvAVG(params, nsims = 50, n.years = 50, strat = selectstrat, hypothesis = "H3")
  simsout4 <- cwd_stoch_wrapper_arrvAVG(params, nsims = 50, n.years = 50, strat = selectstrat, hypothesis = "H4")
  
  simsout4$counts <- simsout1$counts %>% 
    full_join(simsout2$counts, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) %>% 
    full_join(simsout3$counts, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) %>% 
    full_join(simsout4$counts, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) 
  
  simsout4$harvest <- simsout1$harvest %>% 
    full_join(simsout2$harvest, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) %>% 
    full_join(simsout3$harvest, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) %>% 
    full_join(simsout4$harvest, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) 
  
  simsout4$envir <- simsout1$envir %>% 
    full_join(simsout2$envir, by = c("V1", "year","sim")) %>% 
    full_join(simsout3$envir, by = c("V1", "year","sim")) %>% 
    full_join(simsout4$envir, by = c("V1", "year","sim")) 
  
  simsout4$farmyearsout <- simsout1$farmyearsout %>% 
    full_join(simsout2$farmyearsout, by = c("V1", "year","sim")) %>% 
    full_join(simsout3$farmyearsout, by = c("V1", "year","sim")) %>% 
    full_join(simsout4$farmyearsout, by = c("V1", "year","sim")) 
  
  simsout4$detectionsamp <- simsout1$detectionsamp %>% 
    full_join(simsout2$detectionsamp, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) %>% 
    full_join(simsout3$detectionsamp, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) %>% 
    full_join(simsout4$detectionsamp, by = c("fawn.S.N", "juv.f.S.N", "juv.m.S.N", "ad.f.S.N", "ad.m.S.N", "fawn.E.N", "juv.f.E.N", "juv.m.E.N",
                                      "ad.f.E.N", "ad.m.E.N", "juv.f.I.N", "juv.m.I.N", "ad.f.I.N", "ad.m.I.N", "year", "sim")) 
  
  
  return(simsout4)
}

