library(tidyverse)

basic_setup <- tibble(var_e = 0.5) %>% 
  crossing(
    tibble(arms = 2)
  ) %>% 
  crossing(
    # number of blocks per arm
    tibble(m = c(2, 4, 10))
  ) %>% 
  crossing(
    # number of units per block
    tibble(n = c(10, 50, 100))
  ) %>% 
  crossing(
    tibble(var_inter = c(0, 0.1, 0.3))
  ) %>% 
  mutate(N = arms*m*n) %>% 
  crossing(sim_index = 1:10) %>% 
  nest(
    params = c(var_e, arms, m, n, var_inter, N)
  )
# basic_setup$params[[1]]
# psuedo code
# basic_setup %>% 
#   mutate(
#     sim_data = map_df(params, ~ sim_function),
#     sim_estiamte = map_df
#   )

# basic_setup %>% 
#   mutate(
#     sim_data = lmap(params, ~sim_data_generation)
#     )
