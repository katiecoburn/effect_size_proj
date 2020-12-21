trad_es_calc <- function(delta, msw, mswb, mswab, m, n, sigma_error_pop, sigma_block_pop, 
                         sigma_int_pop, sigma_total_pop, ...){
  
  N <- 2 * m * n
  
  sigma_error_est <- msw
  sigma_block_est <- (mswb - msw)/(2*n)
  sigma_int_est <- (mswab - msw)/n
  sigma_total_est <- ss_total/(N - 2)
  
  d2 <- delta/sqrt(sigma_total_est)
  
  v2 <- (2/(m*n)) + (d2^2)/(2 * (N - 2))
  
  return(c(d2 = d2, v2 = v2))
  
}

par_est_calc <- function(delta, msw, mswb, mswab, m, n, sigma_error_pop, sigma_block_pop, 
                         sigma_int_pop, sigma_total_pop, ...){
  # this is the version in excel spreadsheet karina sent on 11/18/20
  
  N <- 2 * m * n
  
  sigma_error_est <- msw
  sigma_block_est <- (mswb - msw)/(2*n)
  sigma_int_est <- (mswab - msw)/n
  sigma_total_est <- ss_total/(N - 2)
  
  d1 <- sqrt(((N - 2 * m) * sigma_error_pop + (m - 1) * (sigma_error_pop + 2 * n * sigma_block_pop) + (m - 1) * (sigma_error_pop + n * sigma_int_pop))/((N - 2) * sigma_total_pop)) * delta/sqrt(sigma_total_est)
  
  v1_term_1 <- (2/(m * n)) * ((sigma_error_pop + n * sigma_int_pop)/sigma_total_pop)
  v1_term_2_num <- ((N - 2 * m) * sigma_error_pop^2 + (m - 1) * (sigma_error_pop + 2 * n * sigma_block_pop)^2 + (m - 1) * (sigma_error_pop + n * sigma_int_pop)^2) * 0.40^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * sigma_error_pop + (m - 1) * (sigma_error_pop + 2 * n * sigma_block_pop) + (m - 1) * (sigma_error_pop + n * sigma_int_pop))^2
  
  v1 <- v1_term_1 + v1_term_2_num/v1_term_2_denom
  
  return(c(d1 = d1, v1 = v1))
  
}

