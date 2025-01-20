
### ### Numbers for consequence table 

# libraries and functions

library(foreach)
library(doParallel)
library(tidyverse)
library(reshape2)
library(patchwork)

source("./SII_model/functions/CreateStableAgeStructure_func.R")
source("./SII_model/functions/popnmodel_func.R")
source("./SII_model/functions/sims_wrapper_func.R")
source("./SII_model/functions/randomsampling_func.R")
source("./SII_model/functions/ComHypowSelectH.R")
source("./SII_model/functions/ComHypowWeights_even.R")

### ### ### Plug and chug method - 
# Have to specify each strategy 

### Populatin MA
## this is suppose to be an average over 25 years (from year 25 - 50)

Trial <- ComHypowSelectH(selectstrat = "SQ", selecthypo = "H1")
## mean of means
PopnMA <- Trial$counts %>% 
  filter(year > 26) %>% 
  rowwise %>% 
  mutate(total_pop = sum(c_across(fawn.S.N:ad.m.I.N), na.rm= TRUE)) %>% 
  ungroup() %>% 
  group_by(sim) %>% 
  summarise(mean_sim_pop =mean(total_pop, na.rm=TRUE)) %>%  #get simulation average
  ungroup() %>% 
  summarise(mean_pop_strat_hypo =mean(mean_sim_pop, na.rm=TRUE))

## Overall mean 
PopnMA <- Trial$counts %>% 
  filter(year > 26) %>% 
  rowwise %>% 
  mutate(total_pop = sum(c_across(fawn.S.N:ad.m.I.N), na.rm= TRUE)) %>% 
  ungroup() %>% 
  summarise(mean_sim_pop =mean(total_pop, na.rm=TRUE))


### Prevalence MA
PrevMA <- Trial$counts %>% 
  filter(year == 50) %>% 
  mutate(noNOTinfected = sum(c_across(fawn.S.N:ad.m.E.N))) %>% ### There is some other way to sum across columns.. 
  mutate(noinfected = sum(c_across(juv.f.I.N:ad.m.I.N ))) %>% 
  mutate(prev = noinfected/sum(c_across(fawn.S.N:ad.m.I.N ))) %>% 
  summarise(mean_sim_prev =mean(prev, na.rm=TRUE))


## Can do the same thing with even weights approach 

### Farm Years MA
FarmYears <- Trial$farmyearsout

FarmYearsMA <- FarmYears %>% 
  mutate(op.years = V1 * year) %>% 
  group_by(sim) %>% 
  summarise(totalnoyrs = sum(op.years))

mean(FarmYearsMA$totalnoyrs)


### time til detection 

year_detected <- Trial$detectionsamp %>% 
  group_by(year, sim) %>% 
  mutate(noinfected_sampled = sum(c_across(juv.f.I.N:ad.m.I.N ))) %>% 
  mutate(yn_detected = case_when(noinfected_sampled >= 1 ~ 1, noinfected_sampled <= 1 ~  0)) %>% 
  mutate(yr_when_detect = yn_detected * year) %>% 
  ungroup() %>% 
  group_by(sim) %>% 
  filter({tmp = min(yr_when_detect);
  yr_when_detect > 0}) %>% 
  ## repeating code to find min value but have to take the 0s out first or 0 is min
  filter({tmp = min(yr_when_detect); 
  yr_when_detect == tmp}) %>% 
  ungroup() %>% 
  summarise(mean_yr_detect = mean(yr_when_detect))
  

### ### ### With even weights on hypotheses... 

Even <- ComHypowWeights_even(selectstrat = "SQ")


## Popn MA - Overall mean 
PopnMA <- Even$counts %>% 
  filter(year > 26) %>% 
  rowwise %>% 
  mutate(total_pop = sum(c_across(fawn.S.N:ad.m.I.N), na.rm= TRUE)) %>% 
  ungroup() %>% 
  summarise(mean_sim_pop =mean(total_pop, na.rm=TRUE))


### Prevalence MA
PrevMA <- Even$counts %>% 
  filter(year == 50) %>% 
  mutate(noNOTinfected = sum(c_across(fawn.S.N:ad.m.E.N))) %>% ### There is some other way to sum across columns.. 
  mutate(noinfected = sum(c_across(juv.f.I.N:ad.m.I.N ))) %>% 
  mutate(prev = noinfected/sum(c_across(fawn.S.N:ad.m.I.N ))) %>% 
  summarise(mean_sim_prev =mean(prev, na.rm=TRUE))


## Can do the same thing with even weights approach 

### Farm Years MA
FarmYears <- Even$farmyearsout

FarmYearsMA <- FarmYears %>% 
  mutate(op.years = V1 * year) %>% 
  group_by(sim) %>% 
  summarise(totalnoyrs = sum(op.years))

mean(FarmYearsMA$totalnoyrs)


### time til detection 

year_detected <- Even$detectionsamp %>% 
  group_by(year, sim) %>% 
  mutate(noinfected_sampled = sum(c_across(juv.f.I.N:ad.m.I.N ))) %>% 
  mutate(yn_detected = case_when(noinfected_sampled >= 1 ~ 1, noinfected_sampled <= 1 ~  0)) %>% 
  mutate(yr_when_detect = yn_detected * year) %>% 
  ungroup() %>% 
  group_by(sim) %>% 
  filter({tmp = min(yr_when_detect);
  yr_when_detect > 0}) %>% 
  ## repeating code to find min value but have to take the 0s out first or 0 is min
  filter({tmp = min(yr_when_detect); 
  yr_when_detect == tmp}) %>% 
  ungroup() %>% 
  summarise(mean_yr_detect = mean(yr_when_detect))



### ### ### Run all at once (parrellize..)


set.seed(666)
setcores <- detectCores() - 1
registerDoParallel(cores=setcores)
setstrats = c("SQ", "Ho", "SK", "PareR", "PventR", "SA")
tempout <- foreach(i = setstrats, .packages = c("dplyr", "tidyr", "stringr", "reshape2")) %dopar% {
  ## Need to load function into loop, might be worth looking into .export          
  source("./SII_model/functions/ComHypowWeights_even.R")
  simsout <- ComHypowWeights_opt2(selectstrat = i)
  simsout
}




### Get counts dataframe  
for(i in 1:length(setstrats)){
  hold <- tempout[[i]]$counts %>% 
    mutate(strat = setstrats[i])
  if(i == 1){
    count_df <- hold
  }else{
    count_df <- count_df %>% full_join(hold)
  }
}

## Popn MA - Overall mean 
# this takes a minute.. 
PopnMA <- count_df %>% 
  filter(year > 26) %>% 
  rowwise %>% 
  mutate(total_pop = sum(c_across(fawn.S.N:ad.m.I.N), na.rm= TRUE)) %>% 
  ungroup() %>% #need to drop rowwise
  group_by(strat) %>% 
  summarise(mean_sim_pop =mean(total_pop, na.rm=TRUE))


### Prevalence MA
PrevMA <- count_df %>% 
  filter(year == 50) %>% 
  group_by(strat, year, sim) %>% 
  mutate(noNOTinfected = sum(c_across(fawn.S.N:ad.m.E.N))) %>% ### There is some other way to sum across columns.. 
  mutate(noinfected = sum(c_across(juv.f.I.N:ad.m.I.N ))) %>% 
  mutate(prev = noinfected/sum(c_across(fawn.S.N:ad.m.I.N ))) %>% 
  ungroup() %>% #need to drop year and sim
  group_by(strat) %>% 
  summarise(mean_sim_prev =mean(prev, na.rm=TRUE))


### Get surveillance dataframe for detection 
for(i in 1:length(setstrats)){
  hold <- tempout[[i]]$detectionsamp %>% 
    mutate(strat = setstrats[i])
  if(i == 1){
    detect_df <- hold
  }else{
    detect_df <- detect_df %>% full_join(hold)
  }
}


year_detected <- detect_df %>% 
  group_by(strat, year, sim) %>% 
  mutate(noinfected_sampled = sum(c_across(juv.f.I.N:ad.m.I.N ))) %>% 
  mutate(yn_detected = case_when(noinfected_sampled >= 1 ~ 1, noinfected_sampled <= 1 ~  0)) %>% 
  mutate(yr_when_detect = yn_detected * year) %>% 
  ungroup() %>% 
  group_by(strat, sim) %>% 
  filter({tmp = min(yr_when_detect);
  yr_when_detect > 0}) %>% 
  ## repeating code to find min value but have to take the 0s out first or 0 is min
  filter({tmp = min(yr_when_detect); 
  yr_when_detect == tmp}) %>% 
  ungroup() %>% 
  group_by(strat) %>% 
  summarise(mean_yr_detect = mean(yr_when_detect))




### Get farmyears
for(i in 1:length(setstrats)){
  hold <- tempout[[i]]$farmyearsout %>% 
    mutate(strat = setstrats[i])
  if(i == 1){
    farmyrs_df <- hold
  }else{
    farmyrs_df <- farmyrs_df %>% full_join(hold)
  }
}


### Farm Years MA

FarmYearsMA <- farmyrs_df %>% 
  mutate(op.years = V1 * year) %>% 
  group_by(strat, sim) %>% 
  summarise(totalnoyrs = sum(op.years)) %>% 
  ungroup() %>% 
  group_by(strat) %>% 
  summarise(mean_farmyrs = mean(totalnoyrs))



