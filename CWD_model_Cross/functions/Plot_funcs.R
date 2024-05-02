#' Plotting functions for processing output from stoch_wrapper 
#' 


plot_stoch_harvest <- function(dat, all.lines, error.bars, detectbar, harvesttype){
  if(missing(dat)==TRUE) warning("missing data to plot")
  if(missing(all.lines)){all.lines = TRUE}
  if(missing(harvesttype)){harvesttype = 1} # default to reporting all harvest 
  
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
           lower_CI = harvest-harvest_sd)
  
  
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



plot_stoch_abundance <- function(dat, all.lines, error.bars){
  if(missing(dat)==TRUE) warning("missing data to plot")
  if(missing(all.lines)){all.lines = TRUE}
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
    dplyr::summarize(avg = mean(n, na.rm = T))
  
  dat.errors <- dat.sum %>%
    group_by(year) %>%
    dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                     hi = quantile(n, 0.975, na.rm=T))
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
