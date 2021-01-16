
##############################################################
  #Start with basic case where sigma_block & sigma_inter are zero
  m=2
  n=10
  N=2*m*n
  arms = 2
  var_e= 1.0
  var_inter = 0
  var_block = 0
  var_total = var_e + var_inter + var_block
  
  var_y = (var_e + n*var_inter)/N
  Y<- rnorm(10000, 0.4, sqrt(var_y))
  
  df_ssw = N - (2*m)
  ssw<- var_e*rchisq(10000, df_ssw)
  
  #df_ssb = m -1
  #ssb<- (var_e + 2*n*var_block)*rchisq(10000, df_ssb)
  
  #df_ssab = m -1
  #ssab<- (var_e + n*var_inter)*rchisq(10000,  df_ssab)
  
  
  sst <- ssw 
  S<- sqrt(sst/(N-2))

  d2 <- Y/S
  
  d1<- (Y/S)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  mean(Y)
  mean(S)
  
  mean(d2)
  mean(d1)
  
  var(d1)
  var(d2)
  
  mean(v1)
  mean(v2)  


  ###########################################################
  ###########################################################
  #Move on to the second case where sigma_block=0.7 & sigma_inter=0
  m=2
  n=100
  N=2*m*n
  arms = 2
  var_e= 0.3
  var_inter = 0
  var_block = 0.7
  var_total = var_e + var_inter + var_block
  
  var_y = (var_e + n*var_inter)/N
  Y<- rnorm(10000, 0.4, sqrt(var_y))
  
  df_ssw = N - (2*m)
  ssw<- var_e*rchisq(10000, df_ssw)
  
  df_ssb = m -1
  ssb<- (var_e + 2*n*var_block)*rchisq(10000, df_ssb)
  
  #df_ssab = m -1
  #ssab<- (var_e + n*var_inter)*rchisq(10000,  df_ssab)
  
  
  sst <- ssw + ssb 
  S<- sqrt(sst/(N-2))
  
  d2 <- Y/S
  
  d1<- (Y/S)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  mean(Y)
  mean(S)
  
  mean(d2)
  mean(d1)
  
  var(d1)
  var(d2)
  
  mean(v1)
  mean(v2)  
  
  ####################################################
  ###########################################################
  ###########################################################
  #Move on to the third case where sigma_block=0.5 & sigma_inter=0.2 & sigma_e=0.3
  m=2
  n=100
  N=2*m*n
  arms = 2
  var_e=0.3
  var_inter = 0.2
  var_block = 0.5
  var_total = var_e + var_inter + var_block
  
  var_y = (var_e + n*var_inter)/N
  Y<- rnorm(10000, 0.4, sqrt(var_y))
  
  df_ssw = N - (2*m)
  ssw<- var_e*rchisq(10000, df_ssw)
  
  df_ssb = m -1
  ssb<- (var_e + 2*n*var_block)*rchisq(10000, df_ssb)
  
  df_ssab = m -1
  ssab<- (var_e + n*var_inter)*rchisq(10000,  df_ssab)
  
  
  sst <- ssw + ssb + ssab
  S<- sqrt(sst/(N-2))
  
  d2 <- Y/S
  
  d1<- (Y/S)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  mean(Y)
  mean(S)
  
  mean(d2)
  mean(d1)
  
  var(d1)
  var(d2)
  
  mean(v1)
  mean(v2)  
  
  
  
  
  