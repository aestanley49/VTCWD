
source("CWD_model_Cross/functions/est_beta_params.R")
source("CWD_model_Cross/functions/allocate_deaths.R")
source("CWD_model_Cross/functions/MyFunctions.R")
source("CWD_model_Cross/cwd_stoch_model_cal_act.r")
source("CWD_model_Cross/ArrivalVectors.R")
source("CWD_model_Cross/functions/cwd_stoch_wrapper_mod.r")

ComHypowSelectH <- function(selectstrat = "SQ", selecthypo = "H1"){
  
  ### Set the parameters based on the selected straetgy 
  if(selectstrat == "SQ" | selectstrat == "SA" | selectstrat == "Ho"){
    params <- list(fawn.an.sur = 0.7, juv.an.sur = 0.8, ad.an.f.sur = 0.82, ad.an.m.sur = 0.8,
                   fawn.repro = 0.06, juv.repro = 1.3, ad.repro = 1.4,
                   hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1, hunt.mort.ad.f = 0.12, hunt.mort.ad.m = 0.50,
                   ini.fawn.prev = 0.01, ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
                   n.age.cats = 12,  n.age.cats.m = 10, n.age.cats.f = 15,
                   p = 0.27, env.foi = 0,  beta.f = 0.028,  beta.m = 0.028,
                   theta = 0.9, n0 = 10000, n.years = 25, rel.risk = 1.0,
                   repro.var = 0.005, fawn.sur.var = 0.005, sur.var = 0.005, hunt.var = 0.0005, juv.sur.var = 0.005, ad.f.sur.var = 0.005, ad.m.sur.var = 0.005, juv.repro.var = 0.005, ad.repro.var = 0.005,
                   WSI = 1,
                   arrival_input = c(0,0,0,0,0,0,0, 0, 0, 0), 
                   Action_young_bucks = 0, Action_lib_harvest = 0, Action_targetrm = 0, Action_sharpshooting = 0,
                   nosampled = 5
    )
    if(selectstrat == "Ho"){
      params$nosampled <- 230 ## Need to overwrite number of samples in surviellance strategy tied to Ho
    }
  }else{ # other strategies have hunting actions turned on
    params <- list(fawn.an.sur = 0.7, juv.an.sur = 0.8, ad.an.f.sur = 0.82, ad.an.m.sur = 0.8,
                   fawn.repro = 0.06, juv.repro = 1.3, ad.repro = 1.4,
                   hunt.mort.fawn = 0.01, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.1, hunt.mort.ad.f = 0.12, hunt.mort.ad.m = 0.50,
                   ini.fawn.prev = 0.01, ini.juv.prev = 0.03, ini.ad.f.prev = 0.04,  ini.ad.m.prev = 0.04,
                   n.age.cats = 12,  n.age.cats.m = 10, n.age.cats.f = 15,
                   p = 0.27, env.foi = 0,  beta.f = 0.028,  beta.m = 0.028,
                   theta = 0.9, n0 = 10000, n.years = 25, rel.risk = 1.0,
                   repro.var = 0.005, fawn.sur.var = 0.005, sur.var = 0.005, hunt.var = 0.0005, juv.sur.var = 0.005, ad.f.sur.var = 0.005, ad.m.sur.var = 0.005, juv.repro.var = 0.005, ad.repro.var = 0.005,
                   WSI = 1,
                   arrival_input = c(0,0,0,0,0,0,0, 0, 0, 0), 
                   Action_young_bucks = 1, Action_lib_harvest = 1, Action_targetrm = 1, Action_sharpshooting = 1,
                   nosampled = 125
    )
    if(selectstrat == "PareR" | selectstrat == "SK"){
      params$nosampled <- 460 ## Need to overwrite number of samples in surviellance strategy tied to PareR and Kitchen Sink
    }
  }
  
  simsout4 <- cwd_stoch_wrapper(params, nsims = 50, n.years = 25, strat = selectstrat, hypothesis = selecthypo)
  
  
  return(simsout4)
}
