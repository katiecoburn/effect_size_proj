
##############################################################
  #Start with basic case where sigma_block & sigma_inter are zero
data_clean <- data.frame()
n_vector<- c(10, 100, 500) 
m_vector <- c(2, 4, 10)

for(i in 1:length(m_vector)){
  
  m= m_vector[i]
  
  for(j in 1:length(n_vector)) {
    
  n= n_vector[j] 
  N=2*m*n
  arms = 2
  var_e= 1.0
  var_inter = 0
  var_block = 0
  var_total = var_e + var_inter + var_block
  
  var_y = 2*(var_e + n*var_inter)/(m*n)
  Y<- rnorm(10000, 0.4, sqrt(var_y))
  
  df_ssw = N - (2*m)
  ssw<- var_e*rchisq(10000, df_ssw)
  
  df_ssb = m -1
  ssb<- (var_e + 2*n*var_block)*rchisq(10000, df_ssb)
  
  df_ssab = m -1
  ssab<- (var_e + n*var_inter)*rchisq(10000,  df_ssab)
  
  df_num<- ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  df_den <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2)
  
  df_total<- df_num/df_den
  
  sst <- ssw + ssb + ssab
  #sst<- ssw
  S<- sqrt(sst/(N-2))

  d2 <- Y/S
  
  d1<- (Y/S)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  c<- 1 - (3/(4*df_total - 1))
    
  g2<- d2*c
  
  g1<- d1*c
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  v2_g<- c^2*v2
  
  v1_g <- c^2*v1
  
  result<- cbind( m, n, 
                  var_block , 
                  var_inter , 
                  var_e,
                  var_total,
                  mean(Y),
                  mean(S),
                  df_total,
                  mean(d1),
                  mean(g1),
                  mean(d2),
                  mean(g2),
                  mean(v1_term_1),
                  mean(v1_term_2),
                  mean(v1),
                  mean(v1_g),
                  mean(v2_term_1),
                  mean(v2_term_2),
                  mean(v2),
                  mean(v2_g),
                  var(d1),
                  var(g1),
                  mean(v1)/var(d1),
                  mean(v1_g)/var(d1),
                  var(d2),
                  var(g2),
                  mean(v2)/var(d2),
                  mean(v2_g)/var(d2))
  
  data_clean<- rbind(data_clean, result)
  
  }}
  
  


  ###########################################################
  ###########################################################
  #Move on to the second case where sigma_block=0.7 & sigma_inter=0

data_clean <- data.frame()
n_vector<- c(10, 100, 500) 
m_vector <- c(2, 4, 10)

for(i in 1:length(m_vector)){
  
  m= m_vector[i]
  
  for(j in 1:length(n_vector)) {
    
  n= n_vector[j]
  N=2*m*n
  arms = 2
  var_e= 0.7
  var_inter = 0
  var_block = 0.3
  var_total = var_e + var_inter + var_block
  
  var_y = 2*(var_e + n*var_inter)/(m*n)
  Y<- rnorm(10000, 0.4, sqrt(var_y))
  
  df_ssw = N - (2*m)
  ssw<- var_e*rchisq(10000, df_ssw)
  
  df_ssb = m -1
  ssb<- (var_e + 2*n*var_block)*rchisq(10000, df_ssb)
  
  df_ssab = m -1
  ssab<- (var_e + n*var_inter)*rchisq(10000,  df_ssab)
  
  df_num<- ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  df_den <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2)
  
  df_total<- df_num/df_den
  
  sst <- ssw + ssb + ssab 
  
  #sst<- ssw + ssb
  S<- sqrt(sst/(N-2))
  
  d2 <- Y/S
  
  d1<- (Y/S)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  c<- 1 - (3/(4*df_total - 1))
  
  g2<- d2*c
  
  g1<- d1*c
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  v2_g<- c^2*v2
  
  v1_g <- c^2*v1
  
  result<- cbind( m, n, 
                  var_block , 
                  var_inter , 
                  var_e,
                  var_total,
                  mean(Y),
                  mean(S),
                  df_total,
                  mean(d1),
                  mean(g1),
                  mean(d2),
                  mean(g2),
                  mean(v1_term_1),
                  mean(v1_term_2),
                  mean(v1),
                  mean(v1_g),
                  mean(v2_term_1),
                  mean(v2_term_2),
                  mean(v2),
                  mean(v2_g),
                  var(d1),
                  var(g1),
                  mean(v1)/var(d1),
                  mean(v1_g)/var(d1),
                  var(d2),
                  var(g2),
                  mean(v2)/var(d2),
                  mean(v2_g)/var(d2))
  
  data_clean<- rbind(data_clean, result)
  
  }}
  
  
  ####################################################
  ###########################################################
  ###########################################################
  #Move on to the third case where sigma_block=0.5 & sigma_inter=0.2 & sigma_e=0.3
  
  data_clean <- data.frame()
  n_vector<- c(10, 100, 500) 
  m_vector <- c(2, 4, 10)
  
  for(i in 1:length(m_vector)){
  
  m= m_vector[i]
  
  for(j in 1:length(n_vector)) {
  
  n= n_vector[j]
  N=2*m*n
  arms = 2
  var_e=0.2
  var_inter = 0.7
  var_block = 0.1
  var_total = var_e + var_inter + var_block
  
  var_y =  2*(var_e + n*var_inter)/(m*n)
  Y<- rnorm(10000, 0.4, sqrt(var_y))
  
  df_ssw = N - (2*m)
  ssw<- var_e*rchisq(10000, df_ssw)
  
  df_ssb = m -1
  ssb<- (var_e + 2*n*var_block)*rchisq(10000, df_ssb)
  
  df_ssab = m -1
  ssab<- (var_e + n*var_inter)*rchisq(10000,  df_ssab)
  
  df_num<- ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  df_den <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2)
  
  df_total<- df_num/df_den
  
  sst <- ssw + ssb + ssab
  S<- sqrt(sst/(N-2))
  
  d2 <- Y/S
  
  d1<- (Y/S)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  c<- 1 - (3/(4*df_total - 1))
  
  g2<- d2*c
  
  g1<- d1*c
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
 
  v2_g<- c^2*v2
  
  v1_g <- c^2*v1
  
  result<- cbind( m, n, 
                  var_block , 
                  var_inter , 
                  var_e,
                  var_total,
                  mean(Y),
                  mean(S),
                 df_total,
                 mean(d1),
                 mean(g1),
                 mean(d2),
                 mean(g2),
                 mean(v1_term_1),
                 mean(v1_term_2),
                 mean(v1),
                 mean(v1_g),
                 mean(v2_term_1),
                 mean(v2_term_2),
                 mean(v2),
                 mean(v2_g),
                 var(d1),
                 var(g1),
                 mean(v1)/var(d1),
                 mean(v1_g)/var(d1),
                 var(d2),
                 var(g2),
                 mean(v2)/var(d2),
                 mean(v2_g)/var(d2))
  
  data_clean<- rbind(data_clean, result)
  
  } }
  