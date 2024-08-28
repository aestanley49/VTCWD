---
  title: "VT_Calibrate_Visual"
output: word_document
date: "2024-03-19"
---
  
  set up
```{r}

rm (list=ls())
library(Hmisc)
library(tidyverse)
library(reshape2)

source("functions/est_beta_params.R")
source("functions/allocate_deaths.R")
source("cwd_stoch_model_calibrate.r")
source("functions/cwd_stoch_wrapper.r")
```

Run sims
```{r params}

#' Randomly allocates deaths in the stochastic CWD model
#'
#' @param deaths A vector of how many die in each age category
#' @param pop A matrix of age.categories = rows and number of I categories = columns
#'
#' @return A new matrix of the # of individuals in each I and age category
#' @export

allocate_deaths <- function(deaths, pop) {
  condition1 <- which(deaths > 0)
  cats <- seq(1, 10, 1)
  
  for (i in condition1) {
    # vector of column locations of length equal to all individuals possible
    c2 <- rep(cats, pop[i, ])
    if (length(c2) == 1) {
      pop[i, c2] <- pop[i, c2] - deaths[i]
    } else {
      # sample these
      c3 <- plyr::count(sample(c2, deaths[i], replace = F))
      # remove those that died
      pop[i, c3$x] <- pop[i, c3$x] - c3$freq
    }
  }
  return(pop)
}

est_beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


cwd_stoch_wrapper <- function(params, nsims) {
  
  if(missing(nsims) == T) warning('nsims not provided')
  if(missing(params) == T) warning('params not provided')
  
  #pre-allocate the output vectors
  counts.sims <- vector("list", nsims)
  deaths.sims <- vector("list", nsims)
  
  for(i in 1:nsims){
    outa <- cwd_stoch_model(params)
    counts.sims[[i]] <- outa$counts
    deaths.sims[[i]] <- outa$deaths
  }
  
  # organize the output into a long data.frame
  counts <- melt(counts.sims, id = c("age", "month", "population", "category",
                                     "year", "sex", "disease")) %>% 
    dplyr::rename(sim = L1)
  
  deaths <- melt(deaths.sims, id = c("age", "month", "population", "category",
                                     "year", "sex")) %>% 
    dplyr::rename(sim = L1)
  
  out <- list(counts = counts, deaths = deaths, f.R0 = outa$f.R0, m.R0 = outa$m.R0)
}

## trying to get 2.5 sex ratio.. can' almost reach 2. Popn is declining, but still have .75 male survival 
# params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.82,
# ad.an.m.sur = 0.1,
# fawn.repro = 0.06, juv.repro = 1.7, ad.repro = 1.8,
# hunt.mort.fawn = 0.015, hunt.mort.juv.f = 0.1, hunt.mort.juv.m = 0.2,
# hunt.mort.ad.f = 0.2, hunt.mort.ad.m = 0.95,
# ini.fawn.prev = 0.00,
# ini.juv.prev = 0.00, ini.ad.f.prev = 0.00,  ini.ad.m.prev = 0.00,
# n.age.cats = 12,  p = 0.27, env.foi = 0,  beta.f = 0.08,  beta.m = 0.08,
# theta = 1, n0 = 1000, n.years = 50, rel.risk = 1.0,
# repro.var = 0.005, fawn.sur.var = 0.005, sur.var = 0.005, hunt.var = 0.005,
# #my added var params
# juv.sur.var = 0.005, ad.f.sur.var = 0.05, ad.m.sur.var = 0.05,
# juv.repro.var = 0.005, ad.repro.var = 0.005
# )


cwd_stoch_model <- function(params) {
  
  # write the list objects to the local environment
  for (v in 1:length(params)) assign(names(params)[v], params[[v]])
  
  ######### CREATE INITIAL CONDITIONS##########
  months <- seq(1, n.years * 12)  # monthly timestep
  hunt.mo <- rep(0, n.years * 12)  # months in where the hunt occurs
  hunt.mo[months%%12 == 7] <- 1  # hunt.mo==1 on Nov
  
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
  
  # group into a vector
  # initial female prevalence
  ini.f.prev <- c(ini.fawn.prev, ini.juv.prev,
                  rep(ini.ad.f.prev, (n.age.cats - 2)))
  # initial male prevalence
  ini.m.prev <- c(ini.fawn.prev, ini.juv.prev,
                  rep(ini.ad.m.prev, (n.age.cats - 2)))
  
  # Create the Leslie Matrix to start the population at stable age dist
  M <- matrix(rep(0, n.age.cats * 2 * n.age.cats * 2), nrow = n.age.cats * 2)
  
  # replace the -1 off-diagonal with the survival rates
  M[row(M) == (col(M) + 1)] <- c(juv.an.sur * (1 - hunt.mort.juv.f),
                                 rep(ad.an.f.sur *  (1 - hunt.mort.ad.f),
                                     n.age.cats - 2), 0,
                                 c(juv.an.sur * (1 - hunt.mort.juv.m),
                                   rep(ad.an.m.sur * (1 - hunt.mort.ad.m),
                                       n.age.cats - 2)))
  
  # if you want the top age category to continue to survive
  M[n.age.cats, n.age.cats] <- ad.an.f.sur * (1 - hunt.mort.ad.f)
  M[n.age.cats * 2, n.age.cats * 2] <- ad.an.m.sur * (1 - hunt.mort.ad.m)
  
  # insert the fecundity vector prebirth census
  M[1, 1:n.age.cats] <- c(0, juv.repro, rep(ad.repro, n.age.cats - 2)) * 0.5 *
    fawn.an.sur * (1 - hunt.mort.fawn)
  M[n.age.cats + 1, 1:n.age.cats] <- M[1, 1:n.age.cats]
  
  # pre-allocate the output matrices
  tmp <- matrix(0, nrow = n.age.cats, ncol = n.years * 12)
  St.f <- tmp  # susceptible female vector
  St.m <- tmp  # suceptible male vector
  # infectious categories
  It.m <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 10))  # females
  It.f <- array(rep(tmp), dim = c(n.age.cats, n.years * 12, 10))  # males
  
  # tracking the # hunted
  Ht.f <- tmp
  Ht.m <- tmp
  # natural deaths
  Dt.f <- tmp
  Dt.m <- tmp
  # disease deaths
  CWDt.f <- tmp
  CWDt.m <- tmp
  
  ########STOCHASTIC WINTER CONDITIONS#########
  # Annual Parameter Draws monthly stochastic survival rates
  #winter.type.draws <- rbinom(n.years, 1, sev.win.prob)
  #winter.type <- rep(winter.type.draws, each=12, times=n.years)
  #############################################
  
  # Intializing with the stable age distribution.
  St.f[, 1] <- round(popbio::stable.stage(M)[1:n.age.cats] * n0 * (1 - ini.f.prev))
  St.m[, 1] <- round(popbio::stable.stage(M)[(n.age.cats + 1):(n.age.cats * 2)] *
                       n0 * (1 - ini.m.prev))
  
  ### ### initializing with wrong sex ratio...
  # sum(St.f[, 1])/sum(St.m[, 1]) 1.5, need 2.5 
  
  if(sum(St.f[,1]) <= 0) {
    warning("These parameters result in a stable age structure with no surviving 
            females.")
  } 
  
  # randomly allocating infecteds across ages and categories.
  It.m[, 1, 1:10] <- rbinom(n.age.cats * 10, round(popbio::stable.stage(M)[1:n.age.cats] *
                                                     n0/10), ini.m.prev)
  It.f[, 1, 1:10] <- rbinom(n.age.cats * 10, round(popbio::stable.stage(M)[1:n.age.cats] *
                                                     n0/10), ini.f.prev)
  
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
    
    # group into a vector
    Sur.f <- c(fawn.sur.draw, juv.sur.draw, rep(ad.f.sur.draw, n.age.cats - 2))
    Sur.m <- c(fawn.sur.draw, juv.sur.draw, rep(ad.m.sur.draw, n.age.cats - 2))
    
    # stochastic hunting survival rates
    hunt.fawn.draw <- rbeta(1, hunt.fawn.b$alpha, hunt.fawn.b$beta, ncp = 0)
    hunt.juv.f.draw <- rbeta(1, hunt.juv.f.b$alpha, hunt.juv.f.b$beta, ncp = 0)
    hunt.juv.m.draw <- rbeta(1, hunt.juv.m.b$alpha, hunt.juv.m.b$beta, ncp = 0)
    hunt.f.draw <- rbeta(n.age.cats - 2, hunt.f.b$alpha, hunt.f.b$beta, ncp = 0)
    hunt.m.draw <- rbeta(n.age.cats - 2, hunt.m.b$alpha, hunt.m.b$beta, ncp = 0)
    
    # on birthdays add in recruits and age everyone by one year also on birthdays do
    # the random parameter draws births happen in June, model starts in May
    if (t%%12 == 2) {
      
      # the last age category remains in place and doesn't die 
      ## *** Changed this to assume 75% mortality and round... 
      St.f[2:(n.age.cats - 1), t] <- St.f[1:(n.age.cats - 2), t - 1]
      St.f[n.age.cats, t] <- round((St.f[n.age.cats, t - 1] +
                                      St.f[(n.age.cats - 1), t - 1]))
      St.m[2:(n.age.cats - 1), t] <- St.m[1:(n.age.cats - 2), t - 1]
      St.m[n.age.cats, t] <- round((St.m[n.age.cats, t - 1] +
                                      St.m[(n.age.cats - 1), t - 1]))
      It.f[2:(n.age.cats - 1), t, ] <- It.f[1:(n.age.cats - 2), t - 1, ]
      It.f[n.age.cats, t, ] <- round((It.f[n.age.cats, t - 1, ] +
                                        It.f[(n.age.cats - 1), t - 1, ]))
      It.m[2:(n.age.cats - 1), t, ] <- It.m[1:(n.age.cats - 2), t - 1, ]
      It.m[n.age.cats, t, ] <- round((It.m[n.age.cats, t - 1, ] +
                                        It.m[(n.age.cats - 1), t - 1, ]))
      
      # reproduction
      I_juv <- sum(It.f[2, t - 1, ])
      I_adults <- sum(It.f[3:n.age.cats, t - 1, ])
      
      fawns_born <- rbinom(1, (St.f[2, t - 1] + I_juv), juv.preg.draw) *2  +
        rbinom(1, (sum(St.f[3:n.age.cats, t - 1]) + I_adults), ad.preg.draw) *2 
      
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
    nat.s.f <- rbinom(n.age.cats, St.f[, t], (1 - Sur.f))
    nat.s.m <- rbinom(n.age.cats, St.m[, t], (1 - Sur.m))
    
    St.f[, t] <- St.f[, t] - nat.s.f
    St.m[, t] <- St.m[, t] - nat.s.m
    # infecteds
    nat.i.f <- matrix(rbinom(length(It.f[, t, ]), size = It.f[, t, ],
                             prob = (1 - Sur.f)), nrow = 12)
    nat.i.m <- matrix(rbinom(length(It.m[, t, ]), size = It.m[, t, ],
                             prob = (1 - Sur.m)), nrow = 12)
    
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
      Ht.f[, t] <- rbinom(n.age.cats, Nt.f, c(hunt.fawn.draw, hunt.juv.f.draw,
                                              hunt.f.draw))
      
      Ht.m[, t] <- rbinom(n.age.cats, Nt.m, c(hunt.fawn.draw, hunt.juv.m.draw,
                                              hunt.m.draw))
      
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
    }
    
    # Disease mortality stochastic movement of individuals from I1 to I2 disease
    # induced mortality here by advancing all I's and only a proportion of the 10th
    # category remains
    I.f.move <- matrix(rbinom(n.age.cats * 10, size = It.f[, t, ], prob = p),
                       nrow = 12)
    I.m.move <- matrix(rbinom(n.age.cats * 10, size = It.m[, t, ], prob = p),
                       nrow = 12)
    
    # store info on those that die directly from disease
    CWDt.f[, t] <- I.f.move[, 10]
    CWDt.m[, t] <- I.m.move[, 10]
    
    # move the I individuals forward in their categories
    It.f[, t, 1] <- It.f[, t, 1] - I.f.move[, 1]
    It.f[, t, 2:10] <- It.f[, t, 2:10] - I.f.move[, 2:10] + I.f.move[, 1:9]
    
    It.m[, t, 1] <- It.m[, t, 1] - I.m.move[, 1]
    It.m[, t, 2:10] <- It.m[, t, 2:10] - I.m.move[, 2:10] + I.m.move[, 1:9]
    
    # Direct transmission considering all I's are equal
    Iall <- sum(It.f[, t, ] + It.m[, t, ])
    Nall <- sum(St.f[, t] + St.m[, t]) + Iall
    
    foi.f <- 1 - exp(-beta.f * Iall/Nall^theta)
    foi.m <- 1 - exp(-beta.m * Iall/Nall^theta)
    
    transmission.f <- rbinom(n.age.cats, St.f[, t], foi.f)
    transmission.m <- rbinom(n.age.cats, St.m[, t], foi.m)
    
    St.f[, t] <- St.f[, t] - transmission.f
    St.m[, t] <- St.m[, t] - transmission.m
    
    # update with the new infections
    It.f[, t, 1] <- transmission.f + It.f[, t, 1]
    It.m[, t, 1] <- transmission.m + It.m[, t, 1]
    
    # Environmental transmission happens last
    envcases.f <- rbinom(n.age.cats, St.f[, t], env.foi)
    envcases.m <- rbinom(n.age.cats, St.m[, t], env.foi)
    
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
  
  # convert the output to long form
  counts.long <- melt(counts) %>%
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, - 1)), disease = "no")
  counts.long$disease[str_sub(counts.long$category, 1, 1) == "I"] <- "yes"
  counts.long$disease <- as.factor(counts.long$disease)
  
  deaths.long <- melt(deaths) %>%
    dplyr::rename(age = Var1, month = Var2, population = value, category = L1) %>%
    mutate(year = (month - 1)/12, sex = as.factor(str_sub(category, - 1)))
  
  output <- list(counts = counts.long, deaths = deaths.long, f.R0 = f.R0, 
                 m.R0 = m.R0)
}



params <- list(fawn.an.sur = 0.6, juv.an.sur = 0.8, ad.an.f.sur = 0.82,
               ad.an.m.sur = 0.8,
               fawn.repro = 0.0, juv.repro = 1.7, ad.repro = 1.8,
               hunt.mort.fawn = 0.015, hunt.mort.juv.f = 0.25, hunt.mort.juv.m = 0.25,
               hunt.mort.ad.f = 0.08, hunt.mort.ad.m = 0.75,
               ini.fawn.prev = 0.00,
               ini.juv.prev = 0.00, ini.ad.f.prev = 0.00,  ini.ad.m.prev = 0.00,
               n.age.cats = 12,  p = 0.27, env.foi = 0,  beta.f = 0.08,  beta.m = 0.08,
               theta = 1, n0 = 143000, n.years = 30, rel.risk = 1.0,
               repro.var = 0.005, fawn.sur.var = 0.05, sur.var = 0.005, hunt.var = 0.005,
               #my added var params
               juv.sur.var = 0.005, ad.f.sur.var = 0.05, ad.m.sur.var = 0.05,
               juv.repro.var = 0.005, ad.repro.var = 0.005
)

simsout <- cwd_stoch_wrapper(params, nsims = 100)


####
simsout.count <- simsout$counts
simsout.count$age.cat <- "adult"
simsout.count$age.cat[simsout.count$age == 1] <- "fawn"

sexratio <- simsout.count %>% 
  filter(age.cat == "adult") %>%
  filter(month %% 12 == 9) %>%
  dplyr::select(-c(category, age, disease)) %>% 
  pivot_wider(names_from = c(sex), values_from = population, values_fn = ~sum(.x, na.rm = TRUE)) %>% 
  mutate(sexratio = f/m)

mean_sim <- sexratio %>% 
  group_by(year) %>% 
  summarise(sexratio = mean(sexratio),
            sexratio_sd = sd(sexratio)) %>% 
  mutate(upper_CI = sexratio+sexratio_sd, 
         lower_CI = sexratio-sexratio_sd)

ggplot(data = sexratio,
       aes(x = year, y = sexratio)) +
  geom_line(aes(group = sim), color = "grey", size = 0.5) +
  geom_line(data = mean_sim, color = "black", size = 1) +
  stat_summary(geom="ribbon", fun.data = mean_cl_boot, 
               conf.int=0.95, alpha = 0.0, linetype = "dashed", color = "black", size = .75)+
  geom_hline(yintercept = 2.5, color = "red") +
  labs(x = "time (years)", 
       y = "sex ratio", title = "Sex Ratio Over Time")+
  theme_classic()


#overall abundance
dat.sum <- simsout.count %>%
  filter(month %% 12 == 9) %>%
  group_by(year, sim) %>%
  dplyr::summarize(n = sum(population)) %>%
  arrange(sim, year)

# calculate the mean
dat.mean <- dat.sum %>%
  group_by(year) %>%
  dplyr::summarize(avg = mean(n, na.rm = T))

dat.errors <- dat.sum %>%
  group_by(year) %>%
  dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                   hi = quantile(n, 0.975, na.rm=T))

ggplot(data = dat.sum, aes(x = year, y = n)) +
  geom_line(color = "grey", aes(group = sim)) +
  geom_line(data = dat.mean, aes(x = year, y = avg), size = 1.5) +
  geom_hline(yintercept = 143758, color = "red") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20))


#buck abundance
dat.sum.m <- simsout.count %>%
  filter(month %% 12 == 9) %>%
  filter(age.cat == "adult") %>%
  filter(sex=="m") %>%
  group_by(year, sim) %>%
  dplyr::summarize(n = sum(population)) %>%
  arrange(sim, year)

# calculate the mean
dat.mean.m <- dat.sum.m %>%
  group_by(year) %>%
  dplyr::summarize(avg = mean(n, na.rm = T))

dat.errors.m <- dat.sum.m %>%
  group_by(year) %>%
  dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                   hi = quantile(n, 0.975, na.rm=T))

ggplot(data = dat.sum.m, aes(x = year, y = n)) +
  geom_line(color = "grey", aes(group = sim)) +
  geom_line(data = dat.mean.m, aes(x = year, y = avg), size = 1.5) +
  geom_hline(yintercept = 27729, color = "red") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20))


#doe abundance
dat.sum.f <- simsout.count %>%
  filter(month %% 12 == 9) %>%
  filter(age.cat == "adult") %>%
  filter(sex=="f") %>%
  group_by(year, sim) %>%
  dplyr::summarize(n = sum(population)) %>%
  arrange(sim, year)

# calculate the mean
dat.mean.f <- dat.sum.f %>%
  group_by(year) %>%
  dplyr::summarize(avg = mean(n, na.rm = T))

dat.errors.f <- dat.sum.f %>%
  group_by(year) %>%
  dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                   hi = quantile(n, 0.975, na.rm=T))

ggplot(data = dat.sum.f, aes(x = year, y = n)) +
  geom_line(color = "grey", aes(group = sim)) +
  geom_line(data = dat.mean.f, aes(x = year, y = avg), size = 1.5) +
  geom_hline(yintercept = 70710, color = "red") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20))


#fawn abundance
dat.sum.fawn <- simsout.count %>%
  filter(month %% 12 == 9) %>%
  filter(age.cat == "fawn") %>%
  group_by(year, sim) %>%
  dplyr::summarize(n = sum(population)) %>%
  arrange(sim, year)

# calculate the mean
dat.mean.fawn <- dat.sum.fawn %>%
  group_by(year) %>%
  dplyr::summarize(avg = mean(n, na.rm = T))

dat.errors.fawn <- dat.sum.fawn %>%
  group_by(year) %>%
  dplyr::summarize(lo = quantile(n, 0.025, na.rm=T),
                   hi = quantile(n, 0.975, na.rm=T))

ggplot(data = dat.sum.fawn, aes(x = year, y = n)) +
  geom_line(color = "grey", aes(group = sim)) +
  geom_line(data = dat.mean.fawn, aes(x = year, y = avg), size = 1.5) +
  geom_hline(yintercept = 45000, color = "red") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20))


# Doe:Fawn ratio
# Highest I can get this is 1.2 and then the abundance goes exponential 
 

doe_fawn_ratio <- simsout$counts %>% 
  mutate(IDfawn = case_when(age == 1 ~ "fawns", age != 1 ~  "notfawn")) %>% 
  filter(!(sex == "m" & IDfawn == "notfawn")) %>% 
  filter(month %% 12 == 9) %>%
  dplyr::select(-c(category, age, sex, disease)) %>% 
  pivot_wider(names_from = c(IDfawn), values_from = population, values_fn = ~sum(.x, na.rm = TRUE)) %>% 
  mutate(ratio = notfawn/fawns)

mean_sim <- doe_fawn_ratio %>% 
  group_by(year) %>% 
  summarise(ratio = mean(ratio),
            ratio_sd = sd(ratio)) %>% 
  mutate(upper_CI = ratio+ratio_sd, 
         lower_CI = ratio-ratio_sd)

ggplot(data = doe_fawn_ratio,
       aes(x = year, y = ratio)) +
  geom_line(aes(group = sim), color = "grey", size = 0.5) +
  geom_line(data = mean_sim, color = "black", size = 1) +
  stat_summary(geom="ribbon", fun.data = mean_cl_boot, 
               conf.int=0.95, alpha = 0.0, linetype = "dashed", color = "black", size = .75)+
  geom_hline(yintercept = 1.5, color = "red") + ## This will depend on pop starting size
  labs(x = "time (years)", 
       y = "Doe:Fawn Ratio", title = "Doe:Fawn Ratio Over Time")+
  theme_classic()

