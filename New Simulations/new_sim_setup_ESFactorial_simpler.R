library(tidyverse)

m <- 2
n <- 20
N <- 2*m*n
var_e <- 1
var_block <- 0.0
var_inter <- 0.0
var_total <- var_e + var_block + var_inter

output <- data.frame(d1 = NA, d2 = NA, v1 = NA, v2 = NA)

for (j in 1:1000){   
  
  ##############################################################
  #Start with basic case where sigma_block & sigma_inter are zero
  
  noise <- rnorm(N, 0, sqrt(var_e))
  group <- c(rep(1, N/2), rep(2, N/2))
  data <- tibble(noise = noise, group = group, Y = NA)
  
  for(i in 1:N){
    if(data$group[i] == 1){
      data$Y[i] = data$noise[i] + 0.40
    }else
      if(data$group[i] == 2){
        data$Y[i] = data$noise[i]
      }
  }
  
  y_trt <- data %>% filter(group == 1) %>% select(Y) %>% .$Y
  y_cnt <- data %>% filter(group == 2) %>% select(Y) %>% .$Y
  
  var_trt <- var(y_trt)
  var_cnt <- var(y_cnt)
  
  #calculate delta
  delta <- mean(y_trt) - mean(y_cnt)
  s <- sqrt(((n*m - 1) * var_trt + (n*m-1) * var_cnt)/(N - 2))
  
  d2 <- delta/s
  d1 <- (delta/s) * sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) * var_total)) 
  
  #Determine v1 and v2
  
  v2 <- (2/(m*n)) + (d2^2)/(2 * (N - 2))
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1 <- v1_term_1 + (v1_term_2_num/v1_term_2_denom)
  
  print(j)
  output[j, 1] <- d1
  output[j, 2] <- d2
  output[j, 3] <- v1
  output[j, 4] <- v2
}

mean(output$d1)
mean(output$d2)

mean(output$v1)
mean(output$v2)

