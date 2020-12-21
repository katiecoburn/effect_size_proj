source("simulation_code_dec_16.R")
source("function_calculations.R")

# 
# our_simulation_function(i_count = 2, j_count = 4, 
#                                 k_count = 100, var_error = 0.50, 
#                                 var_block = 0.40, var_inter = 0.10, reps = 10)

# j_count <- c(2, 4, 10)
# k_count <- c(10, 50, 100)
# var_inter <- c(0, .10, .30)
# j_count <- c(2, 10)
library(tidyverse)
j_count <- 2
k_count <- c(10, 50, 100, 500)
var_inter <- c(-1, 0, .10, .30)
sim_pars <- tibble(expand_grid(j_count, k_count, var_inter)) %>% 
  mutate(i_count = 2, var_error = 0.50, reps = 100)
# sim_pars %>% pmap_dfr(our_simulation_function)
sim_results <- (sim_pars %>% pmap_dfr(our_simulation_function))

sim_results <- tibble(d2 = sim_results$trad$d2, v2 = sim_results$trad$v2, d1 = sim_results$ours$d1, v1 = sim_results$ours$v1, m = sim_results$trad$j_count, n = sim_results$trad$k_count, var_block = sim_results$trad$var_block, var_inter = sim_results$trad$var_inter)


sim_results %>% 
  group_by(m, n, var_block) %>% 
  summarise_at(c("d1", "d2"), var) %>% 
  mutate(actual_var_d1 = d1, actual_var_d2 = d2) %>% 
  select(-c(d1, d2)) %>% 
  cbind( (sim_results %>% group_by(m, n, var_block) %>% summarise_at(c("v1", "v2"), mean) %>% ungroup() %>% select(v1, v2))) %>% 
  mutate(ratio_d1 = v1/actual_var_d1,
         ratio_d2 = v2/actual_var_d2) 

sim_results %>% 
  group_by(n, var_block, m) %>% 
  summarise_at(vars(d1, v1, d2, v2), mean)


write_csv(sim_results, "sim_results_updated_dec21.csv")
