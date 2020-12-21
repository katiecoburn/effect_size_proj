source("sim_setup_from_arend.R")
source("sim_data_generation_karina.R")
source("function_calculations.R")

sim_pars <- basic_setup$params[[1]]
reps <- 10000
results <- tibble(d1 = NA, v1 = NA, d2 = NA, v2 = NA, empir_var = NA, m = NA, n = NA, var_inter = NA)

for(i in 1:nrow(sim_pars)){
  
  params <- sim_pars[i,]
  
  for(j in 1:reps){
    
    data <- sim_data_generation(params)
    results <- rbind(results,
                     c(d1 = data$d1, v1 = data$v1, d2 = data$d2, v2 = data$v2, 
                       empir_var = var(data$Y), m = data$m, n = data$n, 
                       var_inter = data$var_inter))
    
    print(paste("cell: ", i, "rep: ", j))
    
  }
   
}

results <- results[-1,] %>% #drop na row
  tibble()

write_csv(results, file = "karinacoderesults.csv")
