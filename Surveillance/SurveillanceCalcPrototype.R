###### Surveillance data 

## Load libraries 
library(tidyverse)
library(readxl)

## Read in data
countytbl <- read_excel("Surveillance/Tablesfrom2022DeerHarvestReport.xlsx", sheet = "Table012")

PercentSeason <- read_excel("Surveillance/Tablesfrom2022DeerHarvestReport.xlsx", sheet = "PercentSeason")

## Drop extra columns from dataframe 

countytbl <- countytbl %>% select(-Total, -HarvestperSqMile)
countytbl_long <- countytbl %>% 
  pivot_longer(cols = !c(County, Town), names_to = "Season", values_to = "Count", values_drop_na = T)

# join percents to main df by Season
fulldf <- countytbl_long %>% full_join(PercentSeason, by = "Season")

## Assume have 1/3 of total from riffle season 
fulldf <- fulldf %>% 
  mutate(Count = case_when(Season == "Regular" ~ Count*1/3, Season != "Regular" ~ Count))

set.seed(123)
fulldf <- fulldf %>%
  add_column(AdultBuckCount = NA, AdultDoeCount = NA, MaleFawnCount = NA, FemaleFawnCount = NA)
for(i in 1:nrow(fulldf)){
  fulldf[i, 9:12] <- t(rmultinom(1, as.numeric(fulldf[i,4]) , prob = fulldf[i,5:8]))
}



## From this, create 2 different data sets: (aggregated at county level but could do other boundaries)
# 1) All harvested indiviudals
# 2) All sampled inidividuals (just youth and open riffle (which is regular here??))


## Cost points for each sample
# low - Indiana - 16/sample
# High own lab - MI - 189/sample
# High no lab - Ohio - 200/sample

Indiv_county_ALL <- fulldf %>% 
  group_by(County) %>% 
  summarise(AdultBuckCount = sum(AdultBuckCount), 
            AdultDoeCount = sum(AdultDoeCount),
            MaleFawnCount = sum(MaleFawnCount),
            FemaleFawnCount = sum(FemaleFawnCount)) %>% 
  rowwise() %>%
  mutate(ContyTotal = sum(AdultBuckCount, AdultDoeCount, MaleFawnCount, FemaleFawnCount)) %>% 
  mutate(Indiana = 16) %>% 
  mutate(Michigan = 189) %>% 
  mutate(Ohio = 200)

Indiv_county_CheckStations <- fulldf %>% 
  filter(Season == "Youth" | Season == "Regular") %>% 
  group_by(County) %>% 
  summarise(AdultBuckCount = sum(AdultBuckCount), 
            AdultDoeCount = sum(AdultDoeCount),
            MaleFawnCount = sum(MaleFawnCount),
            FemaleFawnCount = sum(FemaleFawnCount)) %>% 
  rowwise() %>%
  mutate(ContyTotal = sum(AdultBuckCount, AdultDoeCount, MaleFawnCount, FemaleFawnCount)) %>% 
  mutate(Indiana = 16) %>% 
  mutate(Michigan = 189) %>% 
  mutate(Ohio = 200)





### ------------ Stratified Random Sample

SRSno.sampleCalcfunc <- function(prevalence = .01, confidence = .99){
  log(-(confidence - 1))/(log(1 - prevalence))
}


SRS_for_dfs_func <- function(Indiv_county_ALL){
Indiv_county_ALL_SRS <- Indiv_county_ALL %>% 
  ## Add prevalence and alphas testing for
  mutate(Prevalence_001 = .001) %>% 
  mutate(Prevalence_01 = .01) %>% 
  mutate(Prevalence_1 = .1) %>% 
  mutate(Alpha_01 = .01) %>% 
  mutate(Alpha_05 = .05) %>% 
  mutate(Alpha_1 = .1) %>% 
  pivot_longer(cols = c("Alpha_01", "Alpha_05", "Alpha_1"), 
               names_to = "Alpha", values_to = "Alpha_level") %>% 
  pivot_longer(cols = c("Prevalence_001", "Prevalence_01", "Prevalence_1"), 
               names_to = "Prevalence", values_to = "Prevalence_level") %>% 
  ### Now calc number of samples needed 
  mutate(NoSamples = log(-((1-Alpha_level) - 1))/(log(1 - Prevalence_level))) %>% 
  mutate(NoSamplesperCounty = NoSamples/nrow(Indiv_county_ALL)) %>% 
  ## Now calculate costs
  mutate(Cost_Indiana = NoSamplesperCounty *Indiana) %>% 
  mutate(Cost_Michigan = NoSamplesperCounty *Michigan) %>% 
  mutate(Cost_Ohio = NoSamplesperCounty *Ohio) %>% 
  ## Do we have enough samples in the county (Y/N)
  mutate(TargetMeet = if_else((NoSamplesperCounty - ContyTotal)  > 0, "False", "True"))
return(Indiv_county_ALL_SRS)
}


Indiv_county_CheckStations_SRS <- SRS_for_dfs_func(Indiv_county_CheckStations)
Indiv_county_ALL_SRS <- SRS_for_dfs_func(Indiv_county_ALL)


ggplot(Indiv_county_CheckStations_SRS) + 
  geom_point(aes(x = Prevalence_level, y = Cost_Indiana, color = TargetMeet, shape = Alpha))+
  facet_wrap(~County)
ggplot(Indiv_county_CheckStations_SRS) + 
  geom_point(aes(x = Prevalence_level, y = Cost_Michigan, color = TargetMeet, shape = Alpha))+
  facet_wrap(~County)
ggplot(Indiv_county_CheckStations_SRS) + 
  geom_point(aes(x = Prevalence_level, y = Cost_Ohio, color = TargetMeet, shape = Alpha))+
  facet_wrap(~County)

## edit df to compare costs on one graph
Indiv_county_CheckStations_SRS_long <- Indiv_county_CheckStations_SRS %>% 
  pivot_longer(cols = c("Cost_Ohio", "Cost_Michigan", "Cost_Indiana"), 
               names_to = "State_Cost", values_to = "Cost")
ggplot(Indiv_county_CheckStations_SRS_long) + 
  geom_point(aes(x = Prevalence_level, y = Cost, color = State_Cost, shape = Alpha))






### ------------- Weighted surveillance

## Function to calculate the number of points given
Target_points_func <- function(testsensitivity = 1, designprevalence, alpha){
  totalPoints <- (-log(alpha))/(testsensitivity*designprevalence)
  return(totalPoints)
}

## For 3 senarios, find total number of points needed and then see if have enough samples in bag
# Use numbers from GA weighted plan 

## Create table with GA numbers 
GAtablenums <- data.frame(
  Sample = c("AdultMale", "AdultFemale", "YearlingMale", "YearlingFemale"),
  Points = c(3, 1.5, 1, 1)
)

TargetPoints <- Target_points_func(designprevalence = .01, alpha = .01) 
## This is for the whole state, need to make county specific
# Could first just divide evenly
# later would want to distribute based on risk


## Senario 1 - prioritize the males for sampling (cost efficient)
PrioritizeMalesFunc <- function(){
  
}
  
## Senario 2 - random sampling for all strata 

## Calculate total number of points available in each county, 
# make column to indicate if that is greater than or less than total needed for senario
# record total number of samples and calculate cost according to Ohio, MI, Inidiana 


cost_success_func <- function(df = Indiv_county_CheckStations, TargetPoints,
                              designprevalence = setprev, alpha = setalpha){
checkme <-  df %>% mutate(AM_points = AdultBuckCount * 3,
                AF_points = AdultDoeCount *1.5,
                YM_points = MaleFawnCount * 1, 
                YF_points = FemaleFawnCount*1) %>% 
  mutate(AvgNoPoints = TargetPoints/nrow(df)) %>% 
    rowwise() %>% 
    mutate(totalPoints = sum(across(c(AM_points, AF_points, YM_points, YF_points)))) %>% 
    mutate(diffPoint = AvgNoPoints - totalPoints) %>% 
    mutate(TargetMeet = if_else(diffPoint > 0, "False", "True")) 
##Need to add row dividing the total amount of points evenly 

##okay so if we have enough points, prioritize the males for sampling (more cost efficient)

## !!!!!!!!!!! Current issue: 
# Have instance where only have 1 observation of young males
# This causes issues where just want to record the 1 and move to next one (current set up)
# because then when don't want to record any values (already have points), adding points don't need.. 
# also note, data is simulated, so can't target spef values 
# !!!!!!!!!
# Could change all counts to + 1 that way don't run into issue?
# This resolved by just making approach greedy - if have more points, then record that sample then stop

for(i in 1:nrow(checkme)){

  if(checkme$TargetMeet[i] == "True"){
    # Select all individuals in AdultMale group * point value and see if that meets target, 
    #record no individuals in new column 
    # if point value is exceeded, remove as many samples as needed until reach rounded point 
    
    ## Iteratively add aduilt bucks
    for(j in 1:max(checkme$AdultBuckCount[i])){
      bucktemp <- j*3
      if(bucktemp == checkme$AvgNoPoints[i]){
        checkme[i,18] <- j
        break
      } else if(bucktemp > checkme$AvgNoPoints[i]){
        checkme[i,18] <- j- 1
        break
      } else { # If don't hit cap
        checkme[i,18] <- j
      }
      
    }
    
    
    ## Iteratively add adult does 
    for(j in 1:max(checkme$AdultDoeCount[i])){
      doetemp <- j*1.5
      holdpoints <- checkme[i,18]*3 + doetemp ## FIX THIS
      if(holdpoints == checkme$AvgNoPoints[i]){
        checkme[i,19] <- j
        break
      } else if(holdpoints > checkme$AvgNoPoints[i]){
        checkme[i,19] <- j - 1
        break
      }
    }
    
    ## Iteratively add yearlings males 
    for(j in 1:max(checkme$MaleFawnCount[i])){
      YearMtemp <- j*1
      holdpoints <- checkme[i,18]*3 + checkme[i,19]*1.5 + YearMtemp
      if(holdpoints == checkme$AvgNoPoints[i]){
        checkme[i,20] <- j
        break
      } else if(holdpoints > checkme$AvgNoPoints[i]){

        checkme[i,20] <- j - 1  ## If only have one J and holdpoints isn't big enough, it gets stuck
        break 
      } else {
        checkme[i,20] <- j
      }
    }
    
    ## Iteratively add yearlings females 
    for(j in 1:max(checkme$FemaleFawnCount[i])){
      YearFtemp <- j*1
      holdpoints <- checkme[i,18]*3 + checkme[i,19]*1.5 + checkme[i,20]*1 + YearFtemp
      if(holdpoints == checkme$AvgNoPoints[i]){
        checkme[i,21] <- j
      } else if(holdpoints > checkme$AvgNoPoints[i]){
          checkme[i,21] <- j - 1
          break
      
      }
    }

      
    
    ## Get total number of samples 
    checkme[i,22] <- checkme[i,18] + checkme[i,19] + checkme[i,20] + checkme[i,21]
  } #close if targetmeet == True loop 
  # So if don't have enough points to meet target, ID how many have what is needed
  else{ ## If maxing out number of samples, record 0 
    ## Note, this is currently not an issue but if it was, 
      #could find the difference between the max number of samples needed
        # and number of samples in hand
   # max(checkme[,22], na.rm = T)
    
    checkme[i,21] <- checkme[i,20] <- checkme[i,19] <- checkme[i,18] <- 0
    checkme[i,22] <- checkme[i,18] + checkme[i,19] + checkme[i,20] + checkme[i,21]
  }
  
} # close row for loop

### Need extra rows calculated to compare across simulations later 
# First rename newey added columns 
colnames(checkme)[18] <- "Samp_Adult_Buck"; colnames(checkme)[19] <- "Samp_Adult_Doe" 
colnames(checkme)[20] <- "Samp_Young_Buck" ; colnames(checkme)[21] <- "Samp_Young_Doe"
colnames(checkme)[22] <- "Total_Samples_allStrata"

checkme2 <- checkme %>% 
  mutate(Prevalence = designprevalence) %>% 
  mutate(Alpha = alpha) %>% 
  ## Now calculate costs
  mutate(Cost_Indiana = Total_Samples_allStrata *Indiana) %>% 
  mutate(Cost_Michigan = Total_Samples_allStrata *Michigan) %>% 
  mutate(Cost_Ohio = Total_Samples_allStrata *Ohio)
  

return(checkme2)

} # close function


## Trying to back calculate the total number of samples taken for each senario


### Set up simulations
EstSampfunc <- function(setprev = .01, setalpha = .01, setdf = Indiv_county_CheckStations){
  TargetPoints <- Target_points_func(designprevalence = setprev, alpha = setalpha) 
  tempdf <- cost_success_func(df = setdf, TargetPoints, designprevalence = setprev, alpha = setalpha)
  return(sum(tempdf[,22]))
}

## Check this is working
ugh <- EstSampfunc(setprev = .01, setalpha = .001, setdf = Indiv_county_CheckStations)

### Create vectors of alpha and prevalence 
vecalpha <- c(.01, .05, .1)
vecprev <- c(.001, .01, .1)


checkstations_NoSampWeight <- data.frame()
setcounter <- 1
  for(i in vecalpha){
    for(j in vecprev){
      HoldNo <- EstSampfunc(setprev = j, setalpha = i, setdf = Indiv_county_CheckStations)
      IDCharString <- paste0("Est", "_", "alpha", i, "prev", j)
      ## Create counter that depends on i & j
      setcounter
      checkstations_NoSampWeight[setcounter, 1] <- IDCharString
      checkstations_NoSampWeight[setcounter, 2] <- HoldNo
      setcounter <- setcounter + 1
    }
  }


All_NoSampWeight <- data.frame()
setcounter <- 1
for(i in vecalpha){
  for(j in vecprev){
    HoldNo <- EstSampfunc(setprev = j, setalpha = i, setdf = Indiv_county_ALL)
    IDCharString <- paste0("Est", "_", "alpha", i, "prev", j)
    ## Create counter that depends on i & j
    setcounter
    All_NoSampWeight[setcounter, 1] <- IDCharString
    All_NoSampWeight[setcounter, 2] <- HoldNo
    setcounter <- setcounter + 1
  }
}


SamplesforWeighted <- cbind(checkstations_NoSampWeight, All_NoSampWeight[,2])
colnames(SamplesforWeighted) <- c("labels", "Checkstations_samples", "All_samples")




### Find way to stack dfs 

Indiv_county_CheckStations_sim_stacked <- data.frame()
for(i in vecalpha){
  for(j in vecprev){
    TargetPoints <- Target_points_func(designprevalence = j, alpha = i) 
    tempdf <- cost_success_func(df = Indiv_county_CheckStations, 
                                TargetPoints, designprevalence = j, alpha = i)
    
    if(i == .01 & j == .001){ ## Note - these aren't dynamic, if list changes need to change
      Indiv_county_CheckStations_sim_stacked <- tempdf
    } else {
      Indiv_county_CheckStations_sim_stacked <- rbind(Indiv_county_CheckStations_sim_stacked, tempdf)
    }
  }
}


### Put together visualization comparison 


ggplot(Indiv_county_CheckStations_sim_stacked) + 
  geom_point(aes(x = Prevalence, y = Cost_Indiana, color = TargetMeet, shape = as.factor(Alpha)))+
  facet_wrap(~County)
ggplot(Indiv_county_CheckStations_sim_stacked) + 
  geom_point(aes(x = Prevalence, y = Cost_Michigan, color = TargetMeet, shape = as.factor(Alpha)))+
  facet_wrap(~County)
ggplot(Indiv_county_CheckStations_sim_stacked) + 
  geom_point(aes(x = Prevalence, y = Cost_Ohio, color = TargetMeet, shape = as.factor(Alpha)))+
  facet_wrap(~County)

## edit df to compare costs on one graph
Indiv_county_CheckStations_sim_stacked_long <- Indiv_county_CheckStations_sim_stacked %>% 
  pivot_longer(cols = c("Cost_Ohio", "Cost_Michigan", "Cost_Indiana"), 
               names_to = "State_Cost", values_to = "Cost")
ggplot(Indiv_county_CheckStations_sim_stacked_long) + 
  geom_point(aes(x = Prevalence, y = Cost, color = State_Cost, shape = as.factor(Alpha)))



### Do some back calculating to figure out prevalence can detect with set no samples and alpha

backcalcprevFunc <- function(samp = sethere, confidence = 1 - setalpha){
  exp(log(-(confidence - 1))/samp)
}

## From no samples available, back calc prevalence for varying alphas

Indiv_county_CheckStations

### Create vectors of alpha and prevalence 
vecalpha <- c(.01, .05, .1)


## Loop isn't working 
checkstations_NoSampWeight <- data.frame()
setcounter <- 1
for(i in vecalpha){
  for(k in 1:(Indiv_county_CheckStations$ContyTotal)){
    HoldNo <- backcalcprevFunc(samp = k, confidence = 1 - i)
    Indiv_county_CheckStations$setalpha <- i 
    Indiv_county_CheckStations$resultingprev <- HoldNo
  }
}

Prev_Indiv_county_CheckStations <- Indiv_county_CheckStations %>% 
  mutate(Alpha_01 = .01) %>% 
  mutate(Alpha_05 = .05) %>%
  mutate(Alpha_1 = .1) %>% 
  pivot_longer(cols = c("Alpha_01", "Alpha_05", "Alpha_1"), 
               names_to = "Alpha", values_to = "Alpha_level") %>% 
  ### Now calc prev from number of samples 
  # extra step to make sure not too wonky 
  mutate(confidencecalc = 1 - Alpha_level) %>% 
  mutate(EstPrev = exp(log(-(confidencecalc - 1))/ContyTotal) )

### How are prevalence values so high? Is this right? 



### -------------  -------------  ------------- Scrap


### ### ### double check results:

## No samples needed if we only sampled adult males (3 points year)
CheckSampfunc <- function(setprev = .01, setalpha = .01, AdultMalesOnly = 1){
  TargetPoints <- Target_points_func(designprevalence = setprev, alpha = setalpha)
  if(AdultMalesOnly == 1)
    checkSamp <- round(TargetPoints/3)
  else{
    checkSamp <- round(TargetPoints/1)
  }
  return(checkSamp)
}

Check_Weight_AdultBucks <- data.frame()
setcounter <- 1
for(i in vecalpha){
  for(j in vecprev){
    HoldNo <- CheckSampfunc(setprev = j, setalpha = i, AdultMalesOnly = 1)
    IDCharString <- paste0("Est", "_", "alpha", i, "prev", j)
    ## Create counter that depends on i & j
    setcounter
    Check_Weight_AdultBucks[setcounter, 1] <- IDCharString
    Check_Weight_AdultBucks[setcounter, 2] <- HoldNo
    setcounter <- setcounter + 1
  }
}


### Manual check 
vecalpha <- c(.01, .05, .1)
vecprev <- c(.001, .01, .1)
#dfs: 
Indiv_county_CheckStations_SRS <- SRS_for_dfs_func(Indiv_county_CheckStations)
Indiv_county_ALL_SRS <- SRS_for_dfs_func(Indiv_county_ALL)
#Est_alpha0.01prev0.001 - almost all are 0 in all county df (besides 1)
#Est_alpha0.05prev0.001 
TargetPoints <- Target_points_func(designprevalence = .1, alpha = .1) 
tempdf <- cost_success_func(df = Indiv_county_ALL, TargetPoints)
