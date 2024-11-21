# Create the Leslie Matrix to start the population at stable age dist

### Need the following parameters but they are loaded into larger popn model
# juv.f.an.sur, hunt.mort.juv.f, ad.an.f.sur, 
# hunt.mort.ad.f, ad.an.m.sur, hunt.mort.ad.m, juv.m.an.sur, 
# hunt.mort.juv.m, juv.repro, ad.repro, fawn.an.sur, hunt.mort.fawn, 
# n0
  
returnstartingpopn <- function(params){
  
  for (v in 1:length(params)){
    assign(names(params)[v], params[[v]]) 
  }
  
  M4 <-matrix(c(0,  (juv.f.an.sur ^ (5/12)  * juv.repro ), (ad.an.f.sur^ (5/12) * ad.repro),     0,     0,
                (0.5*(fawn.an.sur^ (5/12)  *(1 - hunt.mort.fawn))),  0,       0,       0,     0,
                0, (juv.f.an.sur ^ (5/12) * (1 - hunt.mort.juv.f)), (ad.an.f.sur ^ (5/12) * (1 - hunt.mort.ad.f)),     0,     0,
                (0.5*(fawn.an.sur ^ (5/12)*(1 - hunt.mort.fawn))),  0,       0,       0,     0,
                0,       0,       0,       (juv.m.an.sur ^ (5/12) * (1 - hunt.mort.juv.m)),   (ad.an.m.sur ^ (5/12) * (1 - hunt.mort.ad.m))),
              nrow=5,ncol=5,byrow = T)
  prop4<-popbio::stable.stage(M4) 
  
  ## Nick's estimate - fawns = 30%, does = 50%, bucks = 20% 
  startingpopn <- prop4 * n0
  
  # Order: fawns, juv females, adult females, juv males, adult males
  return(startingpopn)
}
