our_simulation_function <- function(i_count, j_count, k_count, var_error,
                                    var_inter, reps, ...){
  print("cell start ")
  
  var_block <- 0.5 - var_inter
  if(var_inter < 0){
    var_inter <- 0
    var_block <- 0
    var_error <- 1
  }
  var_total <- 1
  
  grandmean <- 0
  treatment_effect <- 0.40
  # now will be half of the treatment effect (half of 0.4, for example)
  alpha_group1 <- 0.5 * treatment_effect
  alpha_group2 <- -0.5 * treatment_effect
  m <- j_count
  reps <- reps
  
  constant_group1 <- grandmean + alpha_group1
  constant_group2 <- grandmean + alpha_group2
  
  trad_es_results <- list(d2 = NA, v2 = NA)
  our_es_results <- list(d1 = NA, v1 = NA)
  true_mean <- list(delta = NA)
  interaction_term <- c()
  constant_term <- c()
  
  for(g in 1:reps){
    
    random_betas <- rnorm(n = (m - 1), mean = 0, sd = sqrt(var_block))
    last_beta <- -1 * sum(random_betas)
    
    beta <- c(random_betas, last_beta)
    
    # this (below) is the old version, generating m betas instead of m - 1
    # beta <- rnorm(n = j_count, mean = 0, sd = sqrt(var_block))
    random_inters <- rnorm(n = (m - 1), mean = 0, sd = sqrt(var_inter))
    interaction <- c(random_inters, -1*random_inters, -1*random_inters, random_inters)
    
    # interaction <- rnorm(n = (i_count*j_count), mean = 0, sd = sqrt(var_inter))
    noise <- rnorm(n = (i_count*j_count*k_count), mean = 0, sd = sqrt(var_error))
    block <- rep(seq(1, j_count), ((k_count*i_count*j_count)/j_count))
    group <- c(rep(c(1), (length(noise)/2)), rep(c(2), (length(noise)/2)))
    # For group, 1 is treatment; 2 is control
    
    data <- tibble(y = NA, block = block, group = group, noise = noise)
    
    for(i in 1:(i_count*j_count*k_count)){
      
      group_1_int_terms <- interaction[1:(length(interaction)/2)]
      group_2_int_terms <- interaction[(length(interaction)/2 + 1):length(interaction)]
      
      if(data$group[i] == 1){
         interaction_term[i] <- group_1_int_terms[data$block[i]]
         constant_term[i] <- constant_group1
      }
      if(data$group[i] == 2){
        interaction_term[i] <- group_2_int_terms[data$block[i]]
        constant_term[i] <- constant_group2
      }
      
      data$y[i] <- constant_term[i] + beta[data$block[i]] + data$noise[i] + interaction_term[i]
      
    }
    
    delta <- (data %>% filter(group == 1) %>% .$y %>% mean()) - (data %>% filter(group == 2) %>% .$y %>% mean())
    # This sets delta to the actual difference between means of the generated data (not to the parameter value)
    
    # trad_es_results[[g]] <- trad_es_calc(delta = delta, m = m, n = (i_count*j_count*k_count), sigma_total = 1, sigma_block = var_block, sigma_int = var_inter)
    # Now we do an anova on the generated data to get the mean squares:
    
    msw <- summary(aov(data$y ~ data$block * data$group))[[1]]$`Mean Sq`[4]
    mswb <- summary(aov(data$y ~ data$block * data$group))[[1]]$`Mean Sq`[1]
    mswab <- summary(aov(data$y ~ data$block * data$group))[[1]]$`Mean Sq`[3]
    total_sum_squares <- sum(summary(aov(data$y ~ data$block * data$group))[[1]]$`Sum Sq`)
    
    trad_es_results[[g]] <- trad_es_calc(delta = delta, m = m, n = (i_count*j_count*k_count), msw = msw,
                                         mswb = mswb, mswab = mswab, sigma_total_pop = 1, 
                                         sigma_error_pop = var_error, sigma_block_pop = var_block, 
                                         sigma_int_pop = var_inter, ss_total = total_sum_squares)
    our_es_results[[g]] <- par_est_calc(delta = delta, sigma_total_pop = 1, sigma_block_pop = var_block, 
                                        sigma_int_pop = var_inter, sigma_error_pop = var_error, m = m, 
                                        n = (i_count*j_count*k_count), msw = msw, mswb = mswb, 
                                        mswab = mswab, ss_total = total_sum_squares)
    }
  
  trad_es_results <- bind_rows(trad_es_results) %>% 
    cbind(j_count, k_count, var_inter, var_block)
  our_es_results <- bind_rows(our_es_results) %>% 
    cbind(j_count, k_count, var_inter, var_block)

  return(list(trad = trad_es_results, ours = our_es_results))
  
}