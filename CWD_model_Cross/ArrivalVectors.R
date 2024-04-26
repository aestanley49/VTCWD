### Generate and organize arrival vectors 


## H0 - no arrival ---------------------------
noCWD <- rep(0, 25)

## H1 - carcass ---------------------------

carcassinfectsdeer <- function(rangelow, rangehigh){
  # draw no infected individuals brought back
  set.seed(1)
  nocarcasses <- sample(rangelow:rangehigh, 1)
  carcassinfects <- empty.dump()
  for(i in 1:nocarcasses){
    #make arrival of infections across different years
    seed <- 1 + i
    set.seed(seed)
    temp <- rbinom(25, 1, .01)
    temp <- as.data.frame(t(temp))
    carcassinfects <- rbind(carcassinfects, temp)
  }
  newcarcassinfects <- colSums(carcassinfects)
}


Ho_H1 <- carcassinfectsdeer(rangelow = 5, rangehigh = 10)
SA_H1 <- carcassinfectsdeer(rangelow = 5, rangehigh = 10)
PventR_H1 <- carcassinfectsdeer(rangelow = 10, rangehigh = 20)
PareR_H1 <- carcassinfectsdeer(rangelow = 8, rangehigh = 15)
NoAA_H1 <- carcassinfectsdeer(rangelow = 8, rangehigh = 15)
SK_H1 <- carcassinfectsdeer(rangelow = 5, rangehigh = 10)
SQ_H1 <- carcassinfectsdeer(rangelow = 10, rangehigh = 20)



## H2 - captive ---------------------------

# Note, these are all 0s... 

interpolate_prob_ests <- function(x1, y1, x2, y2, seqlength = 25){
  X <- seq(1, seqlength)
  y = ((y1 - y2)/x1 * X + x2) / 100
  return(y)
}
noinfectdeer <- function(endx = 25, endy, startx = 0, starty, seqlength){
  vec <- c()
  x1 = endx
  y1 = endy
  x2 = startx
  y2 = starty
  seqlength = seqlength
  setprob <- interpolate_prob_ests(x1, y1, x2, y2)
  for(i in setprob){
    set.seed(1000) # need same vec if have the same starting values
    infect <- rbinom(1, 2, prob = i)
    vec <- c(vec, infect)
  }
  return(vec)
}
set.seed(1)
noinfectdeer(25, endy = 25, 0, starty = 23, seqlength = 25)

Ho_H2 <- noinfectdeer(25, endy = 10, 0, starty = 5, seqlength = 25)
SA_H2 <- noinfectdeer(25, endy = 0, 0, starty = 0, seqlength = 25)
PventR_H2 <- noinfectdeer(25, endy = 10, 0, starty = 5, seqlength = 25)
PareR_H2 <- noinfectdeer(25, endy = 10, 0, starty = 1, seqlength = 25)
NoAA_H2 <- noinfectdeer(25, endy = 10, 0, starty = 5, seqlength = 25)
SK_H2 <- noinfectdeer(25, endy = 0, 0, starty = 0, seqlength = 25)
SQ_H2 <- noinfectdeer(25, endy = 10, 0, starty = 5, seqlength = 25)

# for(j in 1:10000000){
# vec <- c()
# for(i in setprob){
#   set.seed(j) # need same vec if have the same starting values
#   infect <- rbinom(1, 2, prob = i)
#   vec <- c(vec, infect)
#   hold <- sum(vec[])
#   break(hold != 0)
# }
# }
### This will always be vectors of 0... 


## H3 - migration ---------------------------

H3 <- c(rep(0, 20), 1, 2, 4, 16, 256)

## H4 - combined ---------------------------

Ho_H4 <- as.vector(Ho_H1) + Ho_H2 + H3
SA_H4 <- as.vector(SA_H1) + SA_H2 + H3
PventR_H4 <- as.vector(PventR_H1) + PventR_H2 + H3
PareR_H4 <- as.vector(PareR_H1) + PareR_H2 + H3
NoAA_H4 <- as.vector(NoAA_H1) + NoAA_H2 + H3
SK_H4 <- as.vector(SK_H1) + SK_H2 + H3
SQ_H4 <- as.vector(SQ_H1) + SQ_H2 + H3

## Put together datatable of arrival vecs ---------------------------

#easiset thing will be to iterate through a list... set up when have iteration code nailed down.
arrivalvec <- list(Ho_H1 = c(as.vector(Ho_H1)), Ho_H2 = c(Ho_H2), H3 = c(H3), H4 = c(Ho_H4))
