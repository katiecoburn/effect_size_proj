sigma_total <- 4.475 # fixed
sigma_error <- 2.12 # fixed
delta <- 0.40 # fixed
int_prop <- c(.10, .30) # specify levels to vary
n <- c(10, 50) # specify levels to vary
m <- c(2, 4)

sim_pars <- tibble( # creating a big tibble of parameter values/combinations
  delta = delta,
  expand_grid(n, int_prop, m) # anything that you're varying gets put in here
  ) %>% 
  mutate(N = 2 * m * n,
         sigma_int = sigma_total * int_prop,
         sigma_total = sigma_total,
         sigma_block = (sigma_total - sigma_error - sigma_int),
         sigma_error = sigma_error,
         msb = sigma_error + sigma_block * 2 * n,
         msab = sigma_error + n * sigma_int) %>% 
  select(delta, m, n, N, int_prop, sigma_total, sigma_int, 
         sigma_block, sigma_error, msb, msab) # organizing it so it's easier to read

# d1 and v1 are ours; d2 and v2 are traditional metrics
sim_results <- bind_cols(sim_pars, 
                         (sim_pars %>% pmap_dfr(par_est_calc)),
                         (sim_pars %>% pmap_dfr(trad_es_calc)))
