
### Basic: doesn't include arrival hypos or select specific strategies 
sims_wrapper <- function(params, nsims) {
  
  if(missing(nsims) == T) warning('nsims not provided')
  if(missing(params) == T) warning('params not provided')
  
  #pre-allocate the output vectors
  counts_sims <- vector("list", nsims)
  harvest_sims <- vector("list", nsims)
  environemnt_sims <- vector("list", nsims)
  diseaseexplore_sims <- vector("list", nsims)
  checkvitals_sims <- vector("list", nsims)
  farmyears_sims <- vector("list", nsims)
  
  for(i in 1:nsims){
    pop <- cwd_stoch_model(params)
    counts_sims[[i]] <- as.data.frame(pop$counts)
    harvest_sims[[i]] <- as.data.frame(pop$removed)
    environemnt_sims[[i]] <- as.data.frame(pop$environemnt)
    diseaseexplore_sims[[i]] <- as.data.frame(pop$diseaseexplore)
    checkvitals_sims[[i]] <- as.data.frame(pop$checkvitals)
    farmyears_sims[[i]] <- as.data.frame(pop$farmyears)
  }
  
  # Reshape the list into a single data frame

  counts <- counts_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise

  harvest <- harvest_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  envir <- environemnt_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  disease <- diseaseexplore_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  vitalscheck <- checkvitals_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise

  farmyearsout <- farmyears_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  out <- list(counts = counts, harvest = harvest, envir = envir, 
              disease = disease, vitalscheck = vitalscheck, farmyearsout = farmyearsout)
}



cwd_stoch_wrapper_arrvAVG <- function(params, nsims, strat, hypothesis, n.years) {
  
  if(missing(nsims) == T) warning('nsims not provided')
  if(missing(params) == T) warning('params not provided')
  
  nsims = nsims
  load("CWD_model_Cross/arrivalvec_AVG.RData")
  ArrivalVec <- ArrivalVec %>% 
    filter(Strategy == strat & Hypo == hypothesis)
  
  #pre-allocate the output vectors
  counts_sims <- vector("list", nsims)
  harvest_sims <- vector("list", nsims)
  environemnt_sims <- vector("list", nsims)
  diseaseexplore_sims <- vector("list", nsims)
  checkvitals_sims <- vector("list", nsims)
  farmyears_sims <- vector("list", nsims)
  detectionsamp_sims <- vector("list", nsims)
  
  for(i in 1:nsims){
    vechere <- ArrivalVec[which(ArrivalVec$sim == i),]
    params$arrival_input <- unlist(vechere$Vec)
    pop <- cwd_stoch_model(params)
    counts_sims[[i]] <- as.data.frame(pop$counts)
    harvest_sims[[i]] <- as.data.frame(pop$removed)
    environemnt_sims[[i]] <- as.data.frame(pop$environemnt)
    diseaseexplore_sims[[i]] <- as.data.frame(pop$diseaseexplore)
    checkvitals_sims[[i]] <- as.data.frame(pop$checkvitals)
    farmyears_sims[[i]] <- as.data.frame(pop$farmyears)
    detectionsamp_sims[[i]] <- as.data.frame(pop$detectionsamp)
  }
  
  # Reshape the list into a single data frame
  
  counts <- counts_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  harvest <- harvest_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  envir <- environemnt_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  disease <- diseaseexplore_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  vitalscheck <- checkvitals_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  farmyearsout <- farmyears_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  detectionsamp <- detectionsamp_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise
  
  
  out <- list(counts = counts, harvest = harvest, envir = envir, disease = disease, 
              detectionsamp = detectionsamp, vitalscheck = vitalscheck, farmyearsout = farmyearsout)
}
