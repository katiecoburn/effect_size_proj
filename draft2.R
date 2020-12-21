library(tidyverse)
source("function_calculations.R")
grandmean <- 0
alpha <- 0.5
constant <- grandmean + alpha
i_count <- 2
j_count <- 2
k_count <- 100
m <- j_count
var_error <- 0.50
var_block <- 0.40
var_inter <- 0.10
var_total <- 1

beta <- rnorm(n = j_count, mean = 0, sd = sqrt(var_block))
interaction <- rnorm(n = i_count*j_count, mean = 0, sd = sqrt(var_inter))
noise <- rnorm(n = i_count*j_count*k_count, mean = 0, sd = sqrt(var_error))
group <- c(rep("T", k_count*j_count), rep("C", k_count*j_count))
block <- rep((c(rep(1, k_count), rep(2, k_count))), m)


data <- tibble(y = NA, block = block, group = group, noise = noise)

for(i in 1:(i_count*j_count*k_count)){
  if(data$group[i] == "T"){
    if(data$block[i] == 1){
      data$y[i] <- constant + beta[1] + interaction[1] + data$noise[i]
    }
    if(data$block[i] == 2){
      data$y[i] <- constant + beta[2] + interaction[2] + data$noise[i]
    }
  }
  if(data$group[i] == "C"){
    if(data$block[i] == 1){
      data$y[i] <- constant + beta[1] + interaction[3] + data$noise[i]
    }
    if(data$block[i] == 2){
      data$y[i] <- constant + beta[2] + interaction[4] + data$noise[i]
    }
  }
}
# data <- data %>% select(-noise)
# 
# summary(lm(y ~ block * group, data = data))
# # summary(aov(y ~ block * group, data = data))

effect_sizes <- tibble(diff = (data %>% filter(group == "T") %>% .$y) - (data %>% filter(group == "C") %>% .$y))

delta <- mean(effect_sizes$diff)

trad_es_calc(delta = delta, m = m, n = (i_count*j_count*k_count), sigma_total = 1)
par_est_calc(delta = delta, sigma_total = 1, sigma_block = var_block, 
             sigma_int = var_inter, msw = var_error, m = m, n = (i_count*j_count*k_count))
