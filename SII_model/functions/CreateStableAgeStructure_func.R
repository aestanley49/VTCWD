# Create the Leslie Matrix to start the population at stable age dist

### Need the following parameters but they are loaded into larger popn model
# juv.f.an.sur, hunt.mort.juv.f, ad.an.f.sur, 
# hunt.mort.ad.f, ad.an.m.sur, hunt.mort.ad.m, juv.m.an.sur, 
# hunt.mort.juv.m, juv.repro, ad.repro, fawn.an.sur, hunt.mort.fawn, 
# n0
  
returnstartingpopn <- function(){
  
  M <- matrix(rep(0, 5*5), nrow = 5) #fawn, juv f, adult f, juv m, adult m
  
  # replace the -1 off-diagonal with the survival rates
  M[row(M) == (col(M) + 1)] <- c(juv.f.an.sur * (1 - hunt.mort.juv.f) * .5 , # 
                                 (ad.an.f.sur *  (1 - hunt.mort.ad.f)), 0,
                                 (ad.an.m.sur * (1 - hunt.mort.ad.m)))
  M[4,1] <-                     (juv.m.an.sur * (1 - hunt.mort.juv.m) * .5) # don't have male and female fawns so this is wonky
  
  
  
  # if you want the top age category to continue to survive
  M[3, 3] <- ad.an.f.sur * (1 - hunt.mort.ad.f)
  M[5, 5] <- ad.an.m.sur * (1 - hunt.mort.ad.m)
  
  # insert the fecundity vector prebirth census - fawns don't have sex so only doing this once
  M[1, 1:3] <- c(0, juv.repro, ad.repro) * fawn.an.sur * (1 - hunt.mort.fawn)
  
  
  ## Code source - Matthipoulos chapter 6 
  det(M - diag(5)) # if non-zero, model has a unique equilibrium 
  
  e <- eigen(M)
  
  de <- e$values[1]
  
  abs(de)
  
  ev <- e$vectors[,1]
  
  startingprop <- round(abs(100*ev/sum(ev)),0)/100
  
  ## Nick's estimate - fawns = 30%, does = 50%, bucks = 20% 
  startingpopn <- startingprop * n0
  
  # Order: fawns, juv females, adult females, juv males, adult males
  return(startingpopn)
}
