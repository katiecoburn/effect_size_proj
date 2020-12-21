sim_results_10 <- read_csv("sim_results_updated_dec16.csv")
sim_results_10 %>% 
  group_by(m, n, var_inter) %>% 
  summarise_at(c("d1", "d2"), var) %>% 
  mutate(actual_var_d1 = d1, actual_var_d2 = d2) %>% 
  select(-c(d1, d2)) %>% 
  cbind( (sim_results_10 %>% group_by(m, n, var_inter) %>% summarise_at(c("v1", "v2"), mean) %>% ungroup() %>% select(v1, v2))) %>% 
  mutate(ratio_d1 = v1/actual_var_d1,
         ratio_d2 = v2/actual_var_d2) 

sim_resultsdec16 <- read_csv("sim_results_updated_dec16.csv")
sim_resultsdec16 %>% 
  group_by(n, var_inter, m) %>% 
  summarise_at(vars(d1, v1, d2, v2), mean)
