#' Plotting functions for processing output from stoch_wrapper 
#' 
library(foreach)
library(doParallel)
library(tidyverse)
library(reshape2)
## loaded in forloop
# source("CWD_model_Cross/functions/ComHypowSelectH.R") # functions for setting parameters with Set Hypothesis

# stupid little function to round to .5 (used in detection piece)
## vectorized... !!! Need to check use of function in harvest code
round_to_half <- function(x) {
  fractional_part <- x - floor(x)
  if_else(fractional_part >= 0.25 & fractional_part < 0.75,
          floor(x) + 0.5,
          if_else(fractional_part <= 0.25,
                  (round(x * 2) / 2) + 0.5,
                  (round(x * 2) / 2) - 0.5))
}



### Select which hypothesis want to run multiple strategies on... -------

CompareAltsFuncwSetH <- function(setstrats = c("SQ", "Ho", "SK", "NoAA"), 
                                 plottype = "harvest", SetH = "H1"){
  
  
  
  ### Set up with a foreach loop to make run faster
  
  setcores <- detectCores() - 1
  registerDoParallel(cores=setcores)
  tempout <- foreach(i = setstrats, .packages = c("dplyr", 
                                                  "tidyr", "stringr", "reshape2")) %dopar% {
                                                    
                                                    ## Need to load function into loop, might be worth looking into .export          
                                                    source("CWD_model_Cross/functions/ComHypowSelectH.R")
                           
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

