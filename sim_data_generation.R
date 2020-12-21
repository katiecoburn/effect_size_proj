sim_data_generation <- function(params, ...){
  
    betas <- rnorm(m, 0, sqrt((0.50 - var_inter)))
    tibble(betas = betas)

}