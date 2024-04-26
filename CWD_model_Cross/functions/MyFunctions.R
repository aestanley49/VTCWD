randomsampling <- function(matrix, target_sum, emptysample) {
  while (sum(emptysample) < target_sum) {
    # Randomly select a non-zero element
    non_zero_indices <- which(matrix > 0)
    selected_index <- sample(non_zero_indices, 1, replace = TRUE)
    
    # Move 1 from the selected element to the new matrix
    matrix[selected_index] <- matrix[selected_index] - 1
    emptysample[selected_index] <- 1 + emptysample[selected_index]
    
    # Update the sum
    if (sum(emptysample) == target_sum) {
      break
    }
  }
  return(emptysample)
}

# here <- randomsampling(matrix = og, target_sum = 8, emptysample = Samp.f[,t ])

## calculate probability that an individual will have CWD from surveillance sample.. 

# identify wether have sampled age classes that could have infected individuals
infectedSamples <- function(sample = Samp.m[,t ], infected = hunted.i.m, 
                            removed = Ht.m[, t], trackInfect = I.Samp.m){
  binarySamp <- sample
  selected_index <- which(binarySamp > 0) 
  binarySamp[selected_index] <- 1
  ProbInfected <- infected/removed * binarySamp # probability of getting infected individual 
  for(i in ProbInfected[!is.na(ProbInfected) & ProbInfected > 0]){ # for every age class where we could have sampled an infected individual 
    getindex <- which(ProbInfected == i )
    Noindividualschoosing <- sample[getindex]
    hold <- rbinom(Noindividualschoosing, size = 1, prob = i)
    trackInfect[getindex, t] <- sum(hold) # keep track of how many infected we find across time and age class
  }
  return(trackInfect)
}
# infectedSamples(sample = Samp.m[,t ], infected = hunted.i.m, 
#                             removed = Ht.m[, t], trackInfect = I.Samp.m)
#   

