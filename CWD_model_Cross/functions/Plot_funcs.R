#' Plotting functions for processing output from stoch_wrapper 
#' 
library(foreach)
library(doParallel)
library(tidyverse)
library(reshape2)
source("CWD_model_Cross/functions/ComHyposFunc.R")

# stupid little function to round to .5 (used in detection piece)
round_to_half <- function(x) {
  integer_part <- floor(x)
  fractional_part <- x - integer_part
  if (fractional_part >= 0.25 & fractional_part < 0.75) {
    return(integer_part + 0.5)
  } else if (fractional_part <= 0.25) {
    return((round(x * 2) / 2)+ .5 ) 
  } else {
    return((round(x * 2) / 2)- .5 )
  }
}

## vectorized... !!! Need to check use of function in harvest code
round_to_half <- function(x) {
  fractional_part <- x - floor(x)
  if_else(fractional_part >= 0.25 & fractional_part < 0.75,
          floor(x) + 0.5,
          if_else(fractional_part <= 0.25,
                  (round(x * 2) / 2) + 0.5,
                  (round(x * 2) / 2) - 0.5))
}




### There is probably some plot could put together that combines multiple simulations/hypothesis outputs
# strat choices: Ho, SA, PventR, PareR, NoAA, SK, SQ

### three plot outputs possible: harvest, abundance, prevalence 

# Note, this will take a looong time to run
### Averages across all hypotheses ------------
CompareAltsFunc <- function(setstrats = c("SQ", "Ho", "SK", "NoAA"), plottype = "harvest"){
  ### Create empty df
  # #df <- empty.dump()
  # df <- data.frame(
  #   sim = c(1),
  #   year = c(1), 
  #   harvest = c(0)
  # )

  # for(i in setstrats){
  #   simsout <- ComHypowWeights_opt2(selectstrat = i)
  # 
  #   ### If reporting the whole harvest (all age & sex classes)
  #   if(harvesttype == 1){
  #     ### All harvest
  #     mod_deaths_harv <- simsout$deaths %>%
  #       filter(month %% 12 == 7) %>%
  #       filter(population > 0) %>%
  #       filter(category == "Ht.f" |category == "Ht.m") %>%
  #       group_by(sim, year) %>%
  #       dplyr::summarise(harvest = sum(population)) %>%
  #       mutate(strategy = i)
  #   }
  # 
  #     df <- df %>% full_join(mod_deaths_harv)
  #     df2 <- df %>% filter(!is.na(strategy)) # remove holder row
  # }

  
  ### Set up with a foreach loop to make run faster
  
  setcores <- detectCores() - 1
  registerDoParallel(cores=setcores)
  tempout <- foreach(i = setstrats, .packages = c("dplyr", 
            "tidyr", "stringr", "reshape2")) %dopar% {
        
          ## Need to load function into loop, might be worth looking into .export          
              source("CWD_model_Cross/functions/ComHyposFunc.R")
              
    simsout <- ComHypowWeights_opt2(selectstrat = i)
 
  ###   Set up dataframes
    ### If reporting the whole harvest (all age & sex classes)
    if(plottype == "harvest"){
      ### All harvest
      mod_deaths_harv <- simsout$deaths %>% 
        filter(month %% 12 == 7) %>%
        filter(population > 0) %>% 
        filter(category == "Ht.f" |category == "Ht.m") %>% 
        group_by(sim, year) %>% 
        dplyr::summarise(harvest = sum(population)) %>% 
        mutate(strategy = i)
      df2 <- mod_deaths_harv %>% filter(!is.na(strategy)) # remove holder row
      df <- simsout$survillance %>% 
        mutate(strategy = i)
      
      returnme <- list(df, df2)
    } else if(plottype == "abundance"){
      
      ### Get counts data 
      # summarize by year and disease
      dat.counts <- simsout$counts
      dat.counts$age.cat <- "adult"
      dat.counts$age.cat[dat.counts$age == 1] <- "fawn"
      
      dat.sum <- dat.counts %>%
        filter(month %% 12 == 11) %>%
        group_by(year, sim) %>%
        dplyr::summarize(n = sum(population)) %>%
        arrange(sim, year) %>% 
        mutate(strategy = i)
      ### get detection data
      df <- simsout$survillance %>% 
        mutate(strategy = i)
      ### output
      returnme <- list(df, dat.sum)
    }else if(plottype == "prevalence") { ## This is prevalence 
      ### prevalence 
      dat <- simsout$counts
      dat.sum <- dat %>%
        filter(month %% 12 == 1) %>%
        group_by(year, sim, disease) %>%
        dplyr::summarize(n = sum(population)) %>%
        spread(key = disease, value = n) %>%
        mutate(prev = yes/ (no + yes)) %>%
        arrange(sim, year) %>% 
        mutate(strategy = i)
      ### get detection data
      df <- simsout$survillance %>% 
        mutate(strategy = i)
      ### output
      returnme <- list(df, dat.sum)
    }
 }
  
  
  ### This code works regardless of selected plot type
  # need to pull dataframes togehter out of list 
  for(i in 1:length(setstrats)){
    hold <- tempout[[i]][[2]] 
    if(i == 1){
      df <- hold
    }else{
      df <- df %>% full_join(hold)
    }
  }
  ### For detection
  for(i in 1:length(setstrats)){
    hold <- tempout[[i]][[1]] 
    if(i == 1){
      dect <- hold
    }else{
      dect <- dect %>% full_join(hold)
    }
  }

### Pull togther dataframes  plots  ------ Harvest
  if(plottype == "harvest"){  
  ## Get average for each strategy (across all hypotheses)
  mean_sim <- df %>% 
    group_by(year, strategy) %>% 
    summarise(harvest_sd = sd(harvest),
              harvest = mean(harvest)) %>% 
    mutate(upper_CI = harvest+harvest_sd, 
           lower_CI = harvest-harvest_sd)
  
  ## Set up dots for arrival
  ### Find average year disease was detected
  detected_prev_sims <- dect %>%
    filter(population > 0) %>%
    group_by(strategy, sim, year, month) %>%
    summarise(detectedprev_count = n()) %>% 
    group_by(strategy, sim) %>%
    summarise(year = min(year))%>% 
    left_join(df, by = c("year", "strategy", "sim"))
  
  detected_prev_mean <- dect %>%
    filter(population > 0) %>%
    group_by(strategy, sim, year, month) %>%
    summarise(detectedprev_count = n()) %>% 
    group_by(strategy, sim) %>%
    summarise(year = min(year))%>% 
    group_by(strategy) %>% 
    summarise(year = round_to_half(mean(year))) %>% 
    left_join(mean_sim, by = c("year", "strategy")) 
  
  
  # Start constructing the plot
 p <-  ggplot(data = df, aes(x = year, y = harvest))  +
    geom_line(aes(group = interaction(sim, strategy), color = strategy), size = 0.5, alpha = 0.5)+
    geom_line(data = mean_sim, aes(x = year, y = harvest, color = as.factor(strategy)), size = 1) +
   geom_point(data = detected_prev_sims, aes(x = year, y = harvest), color = "grey", size = 1) +
   geom_point(data = detected_prev_mean, aes(x = year, y = harvest, color = as.factor(strategy))) +
    stat_summary(aes(color = as.factor(strategy)), geom = "ribbon", fun.data = mean_cl_boot, 
                 conf.int = 0.95, alpha = 0.0, linetype = "dashed", size = 0.75) +
   theme(panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         axis.text=element_text(size=16), axis.title=element_text(size=20)) + 
   xlab("Year") + ylab("Prevalence") + 
   theme_classic(base_size = 14) +
   theme_classic()
  }
  
  

  
  ### Pull together dataframes  plots  ------ Abundance
  if(plottype == "abundance"){  
    dat.mean <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(avg = mean(n, na.rm = T)) %>% 
      mutate(year = round_to_half(year))
    
    dat.errors <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                       hi = quantile(n, 0.975, na.rm=T))
    
    df <- df %>% mutate(year = round_to_half(year))
    
    ## Set up dots for arrival
    ### Find average year disease was detected
    detected_prev_sims <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      left_join(df, by = c("year", "strategy", "sim"))
    
    detected_prev_mean <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      group_by(strategy) %>% 
      summarise(year = round_to_half(mean(year))) %>% 
      left_join(dat.mean, by = c("year", "strategy"))
    
    
    p <- ggplot(data = df, aes(x = year, y = n)) + 
      geom_line(aes(group = interaction(sim, strategy), color = strategy), size = 0.5, alpha = 0.5) +
      geom_line(data = dat.mean, aes(x = year, y = avg, group = NULL, color = as.factor(strategy)),
                linewidth = 1) +
      geom_line(data = dat.errors, aes(x = year, y = lo, group = NULL, color = as.factor(strategy)), 
                linetype = "dashed", linewidth = .75) +
      geom_line(data = dat.errors, aes(x = year, y = hi, group = NULL, color = as.factor(strategy)),
                linetype = "dashed", linewidth = .75) +
      geom_point(data = detected_prev_sims, aes(x = year, y = n), color = "grey", size = 1) +
      geom_point(data = detected_prev_mean, aes(x = year, y = avg, color = as.factor(strategy))) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text=element_text(size=16), axis.title=element_text(size=20)) + 
      xlab("Year") + ylab("Abundance") + 
      theme_classic(base_size = 14) +
      theme_classic()
 
 
    }
  
  
  ### Pull together dataframes  plots  ------ prevalence
  if(plottype == "prevalence"){  
    dat.mean <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(avg.prev = mean(prev, na.rm = T)) %>% 
      mutate(year = round_to_half(year))
    
    dat.errors <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(lo = quantile(prev, 0.025, na.rm=T),
                       hi = quantile(prev, 0.975, na.rm=T))
    
    df <- df %>% mutate(year = round_to_half(year))
    
    ## Set up dots for arrival
    ### Find average year disease was detected
    detected_prev_sims <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      left_join(df, by = c("year", "strategy", "sim"))
    
    detected_prev_mean <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      group_by(strategy) %>% 
      summarise(year = round_to_half(mean(year))) %>% 
      left_join(dat.mean, by = c("year", "strategy"))
    
    
    p <- ggplot(data = df, aes(x = year, y = prev)) + 
      geom_line(aes(group = interaction(sim, strategy), color = strategy), size = 0.5, alpha = 0.5) +
      geom_line(data = dat.mean, aes(x = year, y = avg.prev, group = NULL, color = as.factor(strategy)),
                linewidth = 1) +
      geom_line(data = dat.errors, aes(x = year, y = lo, group = NULL, color = as.factor(strategy)), 
                linetype = "dashed", linewidth = .75) +
      geom_line(data = dat.errors, aes(x = year, y = hi, group = NULL, color = as.factor(strategy)),
                linetype = "dashed", linewidth = .75) +
      geom_point(data = detected_prev_sims, aes(x = year, y = prev), color = "grey", size = 1) +
      geom_point(data = detected_prev_mean, aes(x = year, y = avg.prev, color = as.factor(strategy))) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text=element_text(size=16), axis.title=element_text(size=20)) + 
      xlab("Year") + ylab("Prevalence") + 
      theme_classic(base_size = 14) +
      theme_classic()
    
    
  }
  
  
  
  
  return(p)
}

test <- CompareAltsFunc(setstrats = c("SQ", "Ho", "SK", "NoAA", "PareR", "PventR", "SA"), plottype = "abundance")


### Select which hypothesis want to run multiple strategies on... -------

CompareAltsFuncwSetH <- function(setstrats = c("SQ", "Ho", "SK", "NoAA"), plottype = "harvest", SetH = "H1"){

  
  
  ### Set up with a foreach loop to make run faster
  
  setcores <- detectCores() - 1
  registerDoParallel(cores=setcores)
  tempout <- foreach(i = setstrats, .packages = c("dplyr", 
                                                  "tidyr", "stringr", "reshape2")) %dopar% {
                                                    
                                                    ## Need to load function into loop, might be worth looking into .export          
                                                    source("CWD_model_Cross/functions/ComHyposFunc.R")
                                                    
                                                    simsout <- ComHypowSelectH(selectstrat = i, selecthypo = SetH)
                                                    
                                                    ###   Set up dataframes
                                                    ### If reporting the whole harvest (all age & sex classes)
                                                    if(plottype == "harvest"){
                                                      ### All harvest
                                                      mod_deaths_harv <- simsout$deaths %>% 
                                                        filter(month %% 12 == 7) %>%
                                                        filter(population > 0) %>% 
                                                        filter(category == "Ht.f" |category == "Ht.m") %>% 
                                                        group_by(sim, year) %>% 
                                                        dplyr::summarise(harvest = sum(population)) %>% 
                                                        mutate(strategy = i)
                                                      df2 <- mod_deaths_harv %>% filter(!is.na(strategy)) # remove holder row
                                                      df <- simsout$survillance %>% 
                                                        mutate(strategy = i)
                                                      
                                                      returnme <- list(df, df2)
                                                    } else if(plottype == "abundance"){
                                                      
                                                      ### Get counts data 
                                                      # summarize by year and disease
                                                      dat.counts <- simsout$counts
                                                      dat.counts$age.cat <- "adult"
                                                      dat.counts$age.cat[dat.counts$age == 1] <- "fawn"
                                                      
                                                      dat.sum <- dat.counts %>%
                                                        filter(month %% 12 == 11) %>%
                                                        group_by(year, sim) %>%
                                                        dplyr::summarize(n = sum(population)) %>%
                                                        arrange(sim, year) %>% 
                                                        mutate(strategy = i)
                                                      ### get detection data
                                                      df <- simsout$survillance %>% 
                                                        mutate(strategy = i)
                                                      ### output
                                                      returnme <- list(df, dat.sum)
                                                    }else if(plottype == "prevalence") { ## This is prevalence 
                                                      ### prevalence 
                                                      dat <- simsout$counts
                                                      dat.sum <- dat %>%
                                                        filter(month %% 12 == 1) %>%
                                                        group_by(year, sim, disease) %>%
                                                        dplyr::summarize(n = sum(population)) %>%
                                                        spread(key = disease, value = n) %>%
                                                        mutate(prev = yes/ (no + yes)) %>%
                                                        arrange(sim, year) %>% 
                                                        mutate(strategy = i)
                                                      ### get detection data
                                                      df <- simsout$survillance %>% 
                                                        mutate(strategy = i)
                                                      ### output
                                                      returnme <- list(df, dat.sum)
                                                    }
                                                  }
  
  
  ### This code works regardless of selected plot type
  # need to pull dataframes togehter out of list 
  for(i in 1:length(setstrats)){
    hold <- tempout[[i]][[2]] 
    if(i == 1){
      df <- hold
    }else{
      df <- df %>% full_join(hold)
    }
  }
  ### For detection
  for(i in 1:length(setstrats)){
    hold <- tempout[[i]][[1]] 
    if(i == 1){
      dect <- hold
    }else{
      dect <- dect %>% full_join(hold)
    }
  }
  
  ### Pull togther dataframes  plots  ------ Harvest
  if(plottype == "harvest"){  
    ## Get average for each strategy (across all hypotheses)
    mean_sim <- df %>% 
      group_by(year, strategy) %>% 
      summarise(harvest_sd = sd(harvest),
                harvest = mean(harvest)) %>% 
      mutate(upper_CI = harvest+harvest_sd, 
             lower_CI = harvest-harvest_sd)
    
    ## Set up dots for arrival
    ### Find average year disease was detected
    detected_prev_sims <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      left_join(df, by = c("year", "strategy", "sim"))
    
    detected_prev_mean <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      group_by(strategy) %>% 
      summarise(year = round_to_half(mean(year))) %>% 
      left_join(mean_sim, by = c("year", "strategy")) 
    
    
    # Start constructing the plot
    p <-  ggplot(data = df, aes(x = year, y = harvest))  +
      geom_line(aes(group = interaction(sim, strategy), color = strategy), size = 0.5, alpha = 0.5)+
      geom_line(data = mean_sim, aes(x = year, y = harvest, color = as.factor(strategy)), size = 1) +
      geom_point(data = detected_prev_sims, aes(x = year, y = harvest, fill = as.factor(strategy)), shape = 21, color = "black", size = 3) +
      geom_point(data = detected_prev_mean, aes(x = year, y = harvest, color = as.factor(strategy))) +
      stat_summary(aes(color = as.factor(strategy)), geom = "ribbon", fun.data = mean_cl_boot, 
                   conf.int = 0.95, alpha = 0.0, linetype = "dashed", size = 0.75) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text=element_text(size=16), axis.title=element_text(size=20)) + 
      xlab("Year") + ylab("Prevalence") + 
      theme_classic(base_size = 14) +
      theme_classic()
  }
  
  
  
  
  ### Pull together dataframes  plots  ------ Abundance
  if(plottype == "abundance"){  
    dat.mean <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(avg = mean(n, na.rm = T)) %>% 
      mutate(year = round_to_half(year))
    
    dat.errors <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                       hi = quantile(n, 0.975, na.rm=T))
    
    df <- df %>% mutate(year = round_to_half(year))
    
    ## Set up dots for arrival
    ### Find average year disease was detected
    detected_prev_sims <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      left_join(df, by = c("year", "strategy", "sim"))
    
    detected_prev_mean <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      group_by(strategy) %>% 
      summarise(year = round_to_half(mean(year))) %>% 
      left_join(dat.mean, by = c("year", "strategy"))
    
    
    p <- ggplot(data = df, aes(x = year, y = n)) + 
      geom_line(aes(group = interaction(sim, strategy), color = strategy), size = 0.5, alpha = 0.5) +
      geom_line(data = dat.mean, aes(x = year, y = avg, group = NULL, color = as.factor(strategy)),
                linewidth = 1) +
      geom_line(data = dat.errors, aes(x = year, y = lo, group = NULL, color = as.factor(strategy)), 
                linetype = "dashed", linewidth = .75) +
      geom_line(data = dat.errors, aes(x = year, y = hi, group = NULL, color = as.factor(strategy)),
                linetype = "dashed", linewidth = .75) +
      geom_point(data = detected_prev_sims, aes(x = year, y = n, fill = as.factor(strategy)), shape = 21, color = "black", size = 3) +
      geom_point(data = detected_prev_mean, aes(x = year, y = avg, color = as.factor(strategy))) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text=element_text(size=16), axis.title=element_text(size=20)) + 
      xlab("Year") + ylab("Abundance") + 
      theme_classic(base_size = 14) +
      theme_classic()
    
    
  }
  
  
  ### Pull together dataframes  plots  ------ prevalence
  if(plottype == "prevalence"){  
    dat.mean <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(avg.prev = mean(prev, na.rm = T)) %>% 
      mutate(year = round_to_half(year))
    
    dat.errors <- df %>%
      group_by(strategy, year) %>%
      dplyr::summarize(lo = quantile(prev, 0.025, na.rm=T),
                       hi = quantile(prev, 0.975, na.rm=T))
    
    df <- df %>% mutate(year = round_to_half(year))
    
    ## Set up dots for arrival
    ### Find average year disease was detected
    detected_prev_sims <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      left_join(df, by = c("year", "strategy", "sim"))
    
    detected_prev_mean <- dect %>%
      filter(population > 0) %>%
      group_by(strategy, sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(strategy, sim) %>%
      summarise(year = min(year))%>% 
      group_by(strategy) %>% 
      summarise(year = round_to_half(mean(year))) %>% 
      left_join(dat.mean, by = c("year", "strategy"))
    
    
    p <- ggplot(data = df, aes(x = year, y = prev)) + 
      geom_line(aes(group = interaction(sim, strategy), color = strategy), size = 0.5, alpha = 0.5) +
      geom_line(data = dat.mean, aes(x = year, y = avg.prev, group = NULL, color = as.factor(strategy)),
                linewidth = 1) +
      geom_line(data = dat.errors, aes(x = year, y = lo, group = NULL, color = as.factor(strategy)), 
                linetype = "dashed", linewidth = .75) +
      geom_line(data = dat.errors, aes(x = year, y = hi, group = NULL, color = as.factor(strategy)),
                linetype = "dashed", linewidth = .75) +
      geom_point(data = detected_prev_sims, aes(x = year, y = prev, fill = as.factor(strategy)), shape = 21, color = "black", size = 3) +
      geom_point(data = detected_prev_mean, aes(x = year, y = avg.prev, color = as.factor(strategy))) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text=element_text(size=16), axis.title=element_text(size=20)) + 
      xlab("Year") + ylab("Prevalence") + 
      theme_classic(base_size = 14) +
      theme_classic()
    
    
  }
  
  
  
  
  return(p)
}






plot_stoch_harvest <- function(dat, all.lines, error.bars, detectbar, harvesttype, detect){
  if(missing(dat)==TRUE) warning("missing data to plot")
  if(missing(all.lines)){all.lines = TRUE}
  if(missing(harvesttype)){harvesttype = 1} # default to reporting all harvest 
  if(missing(detect)){detect = TRUE}
  
  ### Find average year disease was detected
  detected_prev <- dat$survillance %>%
    filter(population > 0) %>%
    group_by(sim, year, month) %>%
    summarise(detectedprev_count = n()) %>% 
    group_by(sim) %>%
    summarise(min_year = min(year))
  
  Avg_min_year = mean(detected_prev$min_year) # this is the average value 
  
  
  ### If reporting the whole harvest (all age & sex classes)
  if(harvesttype == 1){
    ### All harvest
    mod_deaths_harv <- dat$deaths %>% 
      filter(month %% 12 == 7) %>%
      filter(population > 0) %>% 
      filter(category == "Ht.f" |category == "Ht.m") %>% 
      group_by(sim, year) %>% 
      dplyr::summarise(harvest = sum(population)) 
    
  }else if (harvesttype == 2 | harvesttype == 3){
    
    mod_deaths_harv <- dat$deaths %>% 
      filter(month %% 12 == 7) %>%
      filter(population > 0) %>% 
      filter(category == "Ht.f" |category == "Ht.m") %>% 
      mutate(IDfawn = case_when(age == 1 ~ 1, age != 1 ~  0)) %>% 
      mutate(IDjuv = case_when(age == 2 ~ 1, age != 2 ~  0)) %>% 
      mutate(IDadult = case_when(age > 2 ~ 1, age < 2 ~  0)) %>% 
      mutate(IDf = case_when(sex == "f" ~ 1, sex != "f" ~  0)) %>% 
      mutate(IDm = case_when(sex == "m" ~ 1, sex != "m" ~  0))
    
    if(harvesttype == 2){ 
      ### Juvinille males
      mod_deaths_harv <- mod_deaths_harv %>% 
        mutate(FawnCount = population * IDfawn) %>% 
        mutate(JuvMCount = population * IDjuv* IDm) %>% 
        mutate(AFCount = population * IDadult * IDf) %>% 
        mutate(AMCount = population * IDadult * IDm)  %>% 
        mutate(age_class = ifelse(age == 1, "fawn",
                                  ifelse(age == 2, "juv", "adult"))) %>% 
        group_by(sim, year, age_class) %>% 
        dplyr::summarise(dead_fawns = sum(FawnCount),
                         dead_juv = sum(JuvMCount),
                         dead_AF = sum(AFCount),
                         dead_AM = sum(AMCount), 
                         .groups = 'keep') %>% 
        mutate_at(vars(dead_fawns, dead_juv, dead_AF, dead_AM), replace_na, 0) %>% 
        
        filter(age_class == "juv") %>% 
        mutate(harvest = dead_juv)
      
    }else if (harvesttype == 3){
      ### Antlerless 
      mod_deaths_harv <- mod_deaths_harv %>% 
        mutate(FawnCount = population * IDfawn) %>% 
        mutate(JuvMCount = population * IDjuv* IDm) %>% 
        mutate(AFCount = population * IDadult * IDf) %>% 
        mutate(AMCount = population * IDadult * IDm)  %>% 
        mutate(age_class = ifelse(age == 1, "fawn",
                                  ifelse(age == 2, "juv", "adult"))) %>% 
        group_by(sim, year, age_class) %>% 
        dplyr::summarise(dead_fawns = sum(FawnCount),
                         dead_juv = sum(JuvMCount),
                         dead_AF = sum(AFCount),
                         dead_AM = sum(AMCount), 
                         .groups = 'keep') %>% 
        mutate_at(vars(dead_fawns, dead_juv, dead_AF, dead_AM), replace_na, 0) %>% 
        ungroup() %>% 
        select(-age_class, -dead_AM) %>% 
        pivot_longer(cols = c("dead_fawns", "dead_juv", "dead_AF"),
                     names_to = "class", values_to =  "popn") %>% 
        group_by(sim, year) %>% 
        summarise(harvest = sum(popn))
      
    }
  }  
  
  # calculate mean, lo and hi percentiles.
  # calculate the mean
  # dat.mean <- dat.sum %>%
  # group_by(year) %>%
  # dplyr::summarize(avg = mean(n, na.rm = T))
  # dat.errors <- dat.sum %>%
  # group_by(year) %>%
  # dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
  #                  hi = quantile(n, 0.975, na.rm=T))
  
  mean_sim <- mod_deaths_harv %>% 
    group_by(year) %>% 
    summarise(harvest = mean(harvest),
              harvest_sd = sd(harvest)) %>% 
    mutate(upper_CI = harvest+harvest_sd, 
           lower_CI = harvest-harvest_sd) %>% 
    mutate(year = round_to_half(year))
  
  
  # Start constructing the plot
  
  p <- ggplot(data = mod_deaths_harv,
              aes(x = year, y = harvest)) +
    geom_line(aes(group = sim), color = "grey", size = 0.5) +
    geom_line(data = mean_sim, color = "black", size = 1) +
    stat_summary(geom="ribbon", fun.data = mean_cl_boot, 
                 conf.int=0.95, alpha = 0.0, linetype = "dashed", color = "black", size = .75)
  
  if(missing(detectbar) == FALSE){
    p <- p+ geom_vline(xintercept = Avg_min_year, color = "red") # add line for detection
  }
  if(detect == TRUE){
    
    mod_deaths_harv <- mod_deaths_harv %>% 
      mutate(year = round_to_half(year))
    
    dect <- dat$survillance
    
    
    detected_prev_sims <- dect %>%
      filter(population > 0) %>%
      group_by(sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(sim) %>%
      summarise(year = min(year))%>% 
      left_join(mod_deaths_harv, by = c("year", "sim"))
    
    detected_prev_mean <- dect %>%
      filter(population > 0) %>%
      group_by(sim, year, month) %>%
      summarise(detectedprev_count = n()) %>% 
      group_by(sim) %>%
      summarise(year = min(year))%>% 
      summarise(year = round_to_half(mean(year))) %>% 
      left_join(mean_sim, by = c("year"))
    
    p <- p + 
      geom_point(data = detected_prev_sims, aes(x = year, y = harvest, group = sim),
                 color = "darkslateblue", size = 1) +
      geom_point(data = detected_prev_mean, aes(x = year, y = harvest, group = NULL), color = "mediumblue")
    
  }
  
  if(harvesttype == 1){
    p <- p +  labs(x = "time (years)", 
                   y = "harvest", title = "Total Harvest per Year")+
      theme_classic() 
  }else if(harvesttype == 2){
    p <- p +  labs(x = "time (years)", 
                   y = "harvest", title = "Total Juvinille Male Harvest per Year")+
      theme_classic() 
    
  }else if(harvesttype == 3){    
    p <- p +  labs(x = "time (years)", 
                   y = "harvest", title = "Total Antlerless Harvest per Year")+
      theme_classic() 
  }
  
  return(p)
  
}

# plot_stoch_harvest(simsout, all.lines, error.bars, detectbar, harvesttype= 3)



plot_stoch_abundance <- function(dat, all.lines, error.bars, detect){
  if(missing(dat)==TRUE) warning("missing data to plot")
  if(missing(all.lines)){all.lines = TRUE}
  if(missing(detect)){detect = TRUE}
  # summarize by year and disease
  dat.counts <- dat$counts
  dat.counts$age.cat <- "adult"
  dat.counts$age.cat[dat.counts$age == 1] <- "fawn"
  
  dat.sum <- dat.counts %>%
    filter(month %% 12 == 11) %>%
    group_by(year, sim) %>%
    dplyr::summarize(n = sum(population)) %>%
    arrange(sim, year)
  
  # calculate mean, lo and hi percentiles.
  # calculate the mean
  dat.mean <- dat.sum %>%
    group_by(year) %>%
    dplyr::summarize(avg = mean(n, na.rm = T)) %>% 
    mutate(year = round_to_half(year))
  
  dat.errors <- dat.sum %>%
    group_by(year) %>%
    dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                     hi = quantile(n, 0.975, na.rm=T))
  
  dat.sum <- dat.sum%>% 
    mutate(year = round_to_half(year))
  
  dect <- dat$survillance
  
  
  detected_prev_sims <- dect %>%
    filter(population > 0) %>%
    group_by(sim, year, month) %>%
    summarise(detectedprev_count = n()) %>% 
    group_by(sim) %>%
    summarise(year = min(year))%>% 
    left_join(dat.sum, by = c("year", "sim"))
  
  detected_prev_mean <- dect %>%
    filter(population > 0) %>%
    group_by(sim, year, month) %>%
    summarise(detectedprev_count = n()) %>% 
    group_by(sim) %>%
    summarise(year = min(year))%>% 
    summarise(year = round_to_half(mean(year))) %>% 
    left_join(dat.mean, by = c("year"))
  
  
  # Start constructing the plot
  if(all.lines == TRUE){
    p <- ggplot(data = dat.sum, aes(x = year, y = n, group = sim)) +
      geom_line(color = "grey") +
      geom_line(data = dat.mean, aes(x = year, y = avg, group = NULL),
                linewidth = 1, color="black")
  }
  if(all.lines == FALSE){
    p <- ggplot(data = dat.mean, aes(x = year, y = avg, group = NULL)) +
      geom_line(linewidth = 1)
  }
  ### Not Currently Working... !!!!
  if(missing(error.bars) == FALSE){
    # plot the error bars
    p <- p + geom_line(data = dat.errors, aes(x = year, y = lo, group = NULL),
                       linetype = "dashed", color = "black") +
      geom_line(data = dat.errors, aes(x = year, y = hi, group = NULL),
                linetype = "dashed", color = "black")
  }
  if(detect == TRUE){
    p <- p + 
      geom_point(data = detected_prev_sims, aes(x = year, y = n, group = sim),
                 color = "darkslateblue", size = 1) +
      geom_point(data = detected_prev_mean, aes(x = year, y = avg, group = NULL), color = "mediumblue")
      
  }
  p <- p + xlab("Year") + ylab("Abundance") + theme_classic(base_size = 14) +
    theme_classic()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text=element_text(size=16), axis.title=element_text(size=20))
  return(p)
  
  
}

# plot_stoch_abundance(simsout, all.lines)


# plot prevalence over time. 
plot_stoch_prev_single <- function(dat, all.lines, error.bars){
  if(missing(dat)==TRUE) warning("missing data to plot")
  if(missing(all.lines)){all.lines = TRUE}
  # summarize by year and disease
  dat.sum <- dat %>%
    filter(month %% 12 == 1) %>%
    group_by(year, sim, disease) %>%
    dplyr::summarize(n = sum(population)) %>%
    spread(key = disease, value = n) %>%
    mutate(prev = yes/ (no + yes)) %>%
    arrange(sim, year)
  # calculate mean, lo and hi percentiles.
  dat.mean <- dat.sum %>%
    group_by(year) %>%
    dplyr::summarize(avg.prev = mean(prev))
  # Start constructing the plot
  if(all.lines == TRUE){
    p <- ggplot(data = dat.sum, aes(x = year, y = prev, group = sim)) +
      geom_line(color = "grey") +
      geom_line(data = dat.mean, aes(x = year, y = avg.prev, group = NULL),
                linewidth = 1, color="black")
  }
  if(all.lines == FALSE){
    p <- ggplot(data = dat.mean, aes(x = year, y = avg.prev, group = NULL)) +
      geom_line(linewidth = 1)
  }
  if(missing(error.bars) == FALSE){
    # calculate mean, lo and hi percentiles.
    dat.mean <- dat.sum %>%
      group_by(year) %>%
      dplyr::summarize(avg.prev = mean(prev),
                       lo = quantile(prev, error.bars[1]),
                       hi = quantile(prev, error.bars[2]))
    # plot the error bars
    p <- p + geom_line(data = dat.mean, aes(x = year, y = lo, group = NULL),
                       linetype = "dashed", color = "black") +
      geom_line(data = dat.mean, aes(x = year, y = hi, group = NULL),
                linetype = "dashed", color = "black")
  }
  p <- p + xlab("Year") + ylab("Prevalence") + theme_classic(base_size = 14) +
    theme_classic()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
  p
}

# plot_stoch_prev_single(simsout$counts, all.lines)
