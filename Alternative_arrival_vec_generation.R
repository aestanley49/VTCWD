
### Must arrive in the 1st year

set.seed(133)
initial <- ArrivalVecsFunc(n.years = 25, nsims = 5000)


mod <- initial
indices_to_remove <- c()

for(i in 1:nrow(mod)){
  if(unlist(mod$Vec[i])[1] != 1){
    indices_to_remove <- c(indices_to_remove, i)
  }
}

# Removing rows based on indices obtained
mod <- mod[-indices_to_remove, ]

## How many sims do we want? 50? 
newarrivalvec <- mod %>% 
  group_by(Strategy, Hypo) %>% 
  sample_n(size = 15, replace = FALSE)  ## Change to 50 when have that in each group
 
# replicating... 
newarrivalvec <- rbind(newarrivalvec, newarrivalvec, newarrivalvec, newarrivalvec) 

ArrivalVec <-   newarrivalvec %>% 
  mutate(sim = row_number())

# save as .r Object to load in later
save(ArrivalVec, file = "CWD_model_Cross/arrivalvec_yr1.RData")

### Arrives in 10th year

mod10 <- initial

## first remove where arrives in first 9 years
for(i in 1:nrow(mod10)){
  if(sum(unlist(mod10$Vec[i])[1:9]) == 1){
    indices_to_remove <- c(indices_to_remove, i)
  }
}
mod10 <- mod10[-indices_to_remove, ]

### then remove rows if not there in 10th year

for(i in 1:nrow(mod10)){
  if(unlist(mod10$Vec[i])[1] != 1){
    indices_to_remove <- c(indices_to_remove, i)
  }
}

# Removing rows based on indices obtained
mod10 <- mod10[-indices_to_remove, ]

newarrivalvec_10 <- mod10 %>% 
  group_by(Strategy, Hypo) %>% 
  sample_n(size = 1, replace = FALSE) %>%  ## Change to 50 when have that in each group
  mutate(sim = row_number())


### Arrives after year 20
# This is basically going to be Hypothesis 3... 

### Average vector for each strategy - hypo combo 

set.seed(133)
initial500 <- ArrivalVecsFunc(n.years = 25, nsims = 500)

stablevec <- function(setstrat, sethypo){
  selectstrat <- initial500 %>% 
    filter(Strategy == setstrat & Hypo == sethypo)
  
  hold <- set_n.years + 4
  
  selectstrat[,5:hold] <- NA
  
  for(i in 1:nrow(selectstrat)){
    tempvec <- as.vector(unlist(selectstrat$Vec[i]))
    selectstrat[i,5:hold] <- t(tempvec)
  }
  
  avgvec <- selectstrat %>%
    mutate(summe = rowSums(select(selectstrat, '...5':'...29'))) %>%
    select(-c(Strategy, Hypo, Vec, sim)) 
  
  noindiv <- round(sum(avgvec$summe)/500)
  hold <- (sum(avgvec$summe)/500)
  
  ugh <- avgvec %>% 
    summarise(across(, ~sum(.)/hold))
  
  probsgen <- as.vector(unlist(ugh[,c(1:25)]))
  
  if(hold == 0){
    returnme <- rep(0, 25)
  }else{
    returnme <- as.vector(t(rmultinom(n = 1, size = noindiv, prob = probsgen)))
  }
  return(returnme)
}

Strategy = c("Ho", "SA", "PventR", "PareR", "NoAA", "SK", "SQ")
Hypo = c("H1", "H2", "H3", "H4")

avgarrivalvec <- tibble(
  Strategy = c("hold"),
  Hypo = c("hold"),
  Vec = list(c())
)

counter <- 1
for(j in Strategy){
  for(k in Hypo){
    avgarrivalvec[counter,1] <- j
    avgarrivalvec[counter,2] <- k
    avgarrivalvec$Vec[counter] <- list(stablevec(setstrat = j, sethypo = k))
    counter <- counter + 1
  }
}

## replicate 50 times
for(i in 1:50){
  avgarrivalvec$sim <- i
  if( i == 1){
    newname <- avgarrivalvec
  }else{
    newname <- rbind(newname, avgarrivalvec) 
  }
}

ArrivalVec <- newname

save(ArrivalVec, file = "CWD_model_Cross/arrivalvec_AVG.RData")




