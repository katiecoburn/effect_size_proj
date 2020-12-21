sim_data_generation <- function(params, ...){
  
  grandmean <- 0
  alpha_g1 <- 0.5
  alpha_g2 <- 0.1
  
  beta <- rnorm(params$m, 0, sqrt((params$var_e - params$var_inter)))
  
  interaction <- rnorm(params$arms*params$m, 0, sqrt(params$var_inter))
  
  noise <- rnorm(params$N, 0, sqrt(params$var_e))

  b1_up= params$N/(params$arms*params$m)
  b2_lw= b1_up + 1
  b2_up= b1_up + params$n
  b3_lw= b2_up + 1
  b3_up= b2_up+ params$n 
  b4_lw= b3_up + 1
  b4_up= b3_up  + params$n
  
  Y_trt1= c()
  Y_trt2= c()
  Y_cnt1=c()
  Y_cnt2 = c()
  
  for(i in 1:b1_up) {
    
    trt_1 = grandmean + alpha_g1 + beta[1] + interaction[1] + noise[i]
    Y_trt1 = rbind(Y_trt1, trt_1)
  }
    
  for(i in b2_lw:b2_up) {
    trt_2 = grandmean + alpha_g1 + beta[2] + interaction[2] + noise[i]
    Y_trt2 = rbind(Y_trt2, trt_2)
  }
  
  for(i in b3_lw:b3_up) {
    cnt_1 = grandmean + alpha_g2 + beta[1] + interaction[3] + noise[i]
    Y_cnt1 = rbind(Y_cnt1, cnt_1)
  }
  
  for(i in b4_lw:b4_up) {
    cnt_2 = grandmean + alpha_g2 + beta[2] + interaction[4] + noise[i]
    Y_cnt2 = rbind(Y_cnt2, cnt_2)
  }

  Y <- c(Y_trt1, Y_trt2, Y_cnt1, Y_cnt2)
  
  group<- c(rep(1, params$N/2), rep(2, params$N/2))
  
  block<- rep(seq(1, params$m), (params$N/params$m))
  
  delta <- mean(c(Y_trt1, Y_trt2)) - mean(c(Y_cnt1, Y_cnt2))
  
   out<- list(m = params$m,
              n = params$n,
              var_inter = params$var_inter,
              beta = beta, 
             interaction = interaction, 
             noise = noise,
             Y_trt1 = Y_trt1, 
             Y_trt2 = Y_trt2, 
             Y_cnt1 = Y_cnt1,
             Y_cnt2 = Y_cnt2,
             group = group,
             block = block,
             Y = Y,
             delta= delta,
             d1 = par_est_calc(delta, 1, params$var_e - params$var_inter, params$var_inter, params$var_e, params$m, params$n)[1] ,
             d2=  trad_es_calc(delta, params$m, params$n, 1)[1],
             v1 = par_est_calc(delta, 1, params$var_e - params$var_inter, params$var_inter, params$var_e, params$m, params$n)[2] ,
             v2=  trad_es_calc(delta, params$m, params$n, 1)[2])
   
  return (out)

}
