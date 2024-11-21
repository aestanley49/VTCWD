sims_wrapper <- function(params, nsims) {
  
  if(missing(nsims) == T) warning('nsims not provided')
  if(missing(params) == T) warning('params not provided')
  
  #pre-allocate the output vectors
  counts_sims <- vector("list", nsims)
  harvest_sims <- vector("list", nsims)
  
  for(i in 1:nsims){
    pop <- cwd_stoch_model(params)
    counts_sims[[i]] <- as.data.frame(pop$counts)
    harvest_sims[[i]] <- as.data.frame(pop$removed)
  }
  
  # Reshape the list into a single data frame

  counts <- counts_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise

  harvest <- harvest_sims %>%
    purrr::imap(~mutate(.x, year = row_number(), sim = .y)) %>%  # Add 'year' and 'sim' columns
    purrr::reduce(bind_rows)  # Combine them row-wise

  
  out <- list(counts = counts, harvest = harvest)
}