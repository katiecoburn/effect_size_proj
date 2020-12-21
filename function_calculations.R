trad_es_calc <- function(delta, msw, mswb, mswab, m, n, sigma_error_pop, sigma_block_pop, 
                         sigma_int_pop, sigma_total_pop, ...){
  
  sigma_error <- msw
  sigma_block <- (mswb - msw)/(2*n)
  sigma_int <- (mswab - msw)/n
  sigma_total <- sigma_error + sigma_block + sigma_int
  
  N <- 2 * m * n
  
  d2 <- delta/sqrt(sigma_total_pop)
  
  v2 <- (2/(m*n)) + (d2^2)/(2 * (N - 2))
  
  return(c(d2 = d2, v2 = v2))
  
}

par_est_calc <- function(delta, msw, mswb, mswab, m, n, sigma_error_pop, sigma_block_pop, 
                         sigma_int_pop, sigma_total_pop, ...){
  # this is the version in excel spreadsheet karina sent on 11/18/20
  
  sigma_error <- msw
  sigma_block <- (mswb - msw)/(2*n)
  sigma_int <- (mswab - msw)/n
  sigma_total <- sigma_error + sigma_block + sigma_int
  
  N <- 2 * m * n
  
  d1 <- sqrt(((N - 2 * m) * sigma_error + (m - 1) * (sigma_error + 2 * n * sigma_block) + (m - 1) * (sigma_error + n * sigma_int))/((N - 2) * sigma_total)) * delta/sqrt(sigma_total_pop)
  
  v1_term_1 <- (2/(m * n)) * ((sigma_error_pop + n * sigma_int_pop)/sigma_total_pop)
  v1_term_2_num <- ((N - 2 * m) * sigma_error^2 + (m - 1) * (sigma_error + 2 * n * sigma_block)^2 + (m - 1) * (sigma_error + n * sigma_int)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * sigma_error + (m - 1) * (sigma_error + 2 * n * sigma_block) + (m - 1) * (sigma_error + n * sigma_int))^2
  
  v1 <- v1_term_1 + v1_term_2_num/v1_term_2_denom
  
  return(c(d1 = d1, v1 = v1))
  
}

