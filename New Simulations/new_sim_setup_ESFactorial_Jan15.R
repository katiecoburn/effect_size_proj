
df_clean <- data.frame()

for (j in 1:1000){   
  
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

grandmean <- 0
alpha_g1 <- (1/2)*0.4
alpha_g2 <- (-1/2)*0.4


rm_beta <- rnorm((m-1), 0, sqrt(var_block))

last_beta<- -1*sum(rm_beta)

beta<- c(rm_beta, last_beta)

rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))

interaction <- c(rm_interaction, -1*rm_interaction, -1*rm_interaction, rm_interaction)  

#rm_noise <- rnorm( (N - 2) , 0, sqrt(var_e))

#last_noise <- -1*sum(rm_noise)

#noise1<- c(rm_noise[1:199], last_noise)

#noise2 <- c(rm_noise[200:398], -last_noise)

#noise<- c(noise1, noise2)

#noise<- c(rm_noise,-last_noise, last_noise)

noise <- rnorm(N, 0, sqrt(var_e))

b1_up= N/(arms*m)
b2_lw= b1_up + 1
b2_up= b1_up + n
b3_lw= b2_up + 1
b3_up= b2_up + n 
b4_lw= b3_up + 1
b4_up= b3_up + n

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

y_trt<- c(Y_trt1, Y_trt2)
y_cnt<- c(Y_cnt1, Y_cnt2)

var_trt<- var(y_trt)
var_cnt<- var(y_cnt)

group<- c(rep(1, N/2), rep(2, N/2))

block<- rep(seq(1, m), (N/m))

#calculate delta
delta <- mean(y_trt) - mean(y_cnt)

#Estimate SS using formula
# s<- sqrt(((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2))  #This is: S= sqrt(SS_error/(N-2))

### Added KC 1/15
# SS as sum of chi-squared distributions

ssw <- var_e*rchisq(n = 1, df = (2 * m * n - 2 * m))
ssb <- (var_e + 2*n*var_block)*rchisq(n = 1, df = (m - 1))
ssab <- (var_e + n*var_inter)*rchisq(n = 1, df = (m - 1))

s <- sqrt((ssw + ssb + ssab)/(N - 2))

### Added KC 1/15

#this matches what we get from summary(lm1)$sigma, so we are good

#Determine d1 and d2
d2 <- delta/s

d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 

#g2 <- d2 * (1 - (3/(4*(N-2) -1)))

#g1 <- d1 * (1 - (3/(4*(N-2) -1)))

#d1 & d2 match and are "close" to 0.4

#Determine v1 and v2

v2_term_1 <- (2/(m*n))

v2_term_2 <- (d2^2)/(2 * (N - 2))

v2 <- v2_term_1  + v2_term_2

v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2

v1_term_2 <- v1_term_2_num/v1_term_2_denom

v1 <- v1_term_1 + v1_term_2

#v1 & v2 match

#All good so far

df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)

df_clean <- rbind(df_clean, df_clean_current)

names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")

}

mean(df_clean$Y)  #0.997
mean(df_clean$S)  #0.9949

mean(df_clean$d1)  #0.408
mean(df_clean$d2) #0.408

mean(df_clean$v1)  #0.103
mean(df_clean$v2)  #0.103

mean(df_clean$v1_1)  #0.103
mean(df_clean$v2_1)

mean(df_clean$v1_2)  #0.103
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  #0.07 (actual variance)

var(df_clean$d1)
var(df_clean$d2)

###########################################################
###########################################################
#Move on to the second case where sigma_block=0.7 & sigma_inter=0
df_clean <- data.frame()

for (j in 1:10000){   
m=2
n=100
N=2*m*n
arms = 2
var_e= 0.3
var_inter = 0
var_block = 0.7
var_total = var_e + var_inter + var_block

grandmean <- 0
alpha_g1 <- (1/2)*0.4
alpha_g2 <- (-1/2)*0.4


rm_beta <- rnorm((m-1), 0, sqrt(var_block))

last_beta<- -1*sum(rm_beta)

beta<- c(rm_beta, last_beta)

rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))

interaction <- c(rm_interaction, -1*rm_interaction, -1*rm_interaction, rm_interaction)  

noise <- rnorm( N , 0, sqrt(var_e))

#rm_noise <- rnorm( (N - 2) , 0, sqrt(var_e))

#last_noise <- -1*sum(rm_noise)

#noise1<- c(rm_noise[1:199], last_noise)

#noise2 <- c(rm_noise[200:398], -last_noise)

#noise<- c(noise1, noise2)


b1_up= N/(arms*m)
b2_lw= b1_up + 1
b2_up= b1_up + n
b3_lw= b2_up + 1
b3_up= b2_up + n 
b4_lw= b3_up + 1
b4_up= b3_up + n


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

y_trt<- c(Y_trt1, Y_trt2)
y_cnt<- c(Y_cnt1, Y_cnt2)

y_block1<- c(Y_trt1, Y_cnt1)
y_block2<- c(Y_trt2, Y_cnt2)

var_trt<- var(y_trt)
var_cnt<- var(y_cnt)

grandmeanpop<- mean(Y)

group<- c(rep(1, N/2), rep(2, N/2))

block<- rep(seq(1, m), (N/m))

#calculate delta
delta <- mean(y_trt) - mean(y_cnt)

#Calculate SS using formulas
# ss1 <- ((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2)  #SS Error/(N-2)
# 
# ss2<- (n*arms*((mean(y_block1) - grandmeanpop)^2  + (mean(y_block2) - grandmeanpop)^2))/(N-2) #SS Block/(N-2)
# 
# s<- sqrt(ss1 + ss2) # S = sqrt(SS/(N-2))

### Added KC 1/15
# SS as sum of chi-squared distributions

ssw <- var_e*rchisq(n = 1, df = (2 * m * n - 2 * m))
ssb <- (var_e + 2*n*var_block)*rchisq(n = 1, df = (m - 1))
ssab <- (var_e + n*var_inter)*rchisq(n = 1, df = (m - 1))

s <- sqrt((ssw + ssb + ssab)/(N - 2))

### Added KC 1/15


#Determine d1 and d2
d2 <- delta/s

d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 

#d1 & d2 are both close to 0.4, but d1 is closer. Great news

#Determine v1 and v2

v2_term_1 <- (2/(m*n))

v2_term_2 <- (d2^2)/(2 * (N - 2))

v2 <- v2_term_1  + v2_term_2

v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2

v1_term_2 <- v1_term_2_num/v1_term_2_denom

v1 <- v1_term_1 + v1_term_2

#v1 tends to be smaller than v2, this is good news too

df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)

df_clean <- rbind(df_clean, df_clean_current)

names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")

}

mean(df_clean$Y)  
mean(df_clean$S)  

mean(df_clean$d1)  
mean(df_clean$d2) 

mean(df_clean$v1)  
mean(df_clean$v2)  

mean(df_clean$v1_1)  
mean(df_clean$v2_1)

mean(df_clean$v1_2)  
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  

(mean(df_clean$v2))/(var(df_clean$d2)) 

var(df_clean$d1)
var(df_clean$d2)

####################################################
###########################################################
###########################################################
#Move on to the third case where sigma_block=0.5 & sigma_inter=0.2 & sigma_e=0.3

df_clean <- data.frame()

for (j in 1:10000){ 
m=2
n=100
N=2*m*n
arms = 2
var_e=0.3
var_inter = 0.2
var_block = 0.5
var_total = var_e + var_inter + var_block

grandmean <- 0
alpha_g1 <- (1/2)*0.4
alpha_g2 <- (-1/2)*0.4


rm_beta <- rnorm((m-1), 0, sqrt((var_block)))

last_beta<- -1*sum(rm_beta)

beta<- c(rm_beta, last_beta)

rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))

interaction <- c(rm_interaction, -1*rm_interaction,  rm_interaction, -1*rm_interaction)  

noise <- rnorm( N , 0, sqrt(var_e))

b1_up= N/(arms*m)
b2_lw= b1_up + 1
b2_up= b1_up + n
b3_lw= b2_up + 1
b3_up= b2_up + n 
b4_lw= b3_up + 1
b4_up= b3_up + n

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

y_trt<- c(Y_trt1, Y_trt2)
y_cnt<- c(Y_cnt1, Y_cnt2)

y_block1<- c(Y_trt1, Y_cnt1)
y_block2<- c(Y_trt2, Y_cnt2)

var_trt<- var(y_trt)
var_cnt<- var(y_cnt)

grandmeanpop<- mean(Y)

group<- c(rep(1, N/2), rep(2, N/2))

block<- rep(seq(1, m), (N/m))

#calculate delta
delta <- mean(y_trt) - mean(y_cnt)

# ss1 <- ((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2)  #SS_Error/(N-2)
# 
# ss2<- (n*arms*((mean(y_block1) - grandmeanpop)^2  + (mean(y_block2) - grandmeanpop)^2))/(N-2) #SS_Block/(N-2)
# 
# ss3<- (n*(  (mean(Y_trt1) - mean(y_trt) - mean(y_block1) + grandmeanpop)^2 
#           + (mean(Y_trt2)   - mean(y_trt) - mean(y_block2) + grandmeanpop)^2 
#           + (mean(Y_cnt1) - mean(y_cnt) - mean(y_block1) + grandmeanpop)^2 
#           + (mean(Y_cnt2) - mean(y_cnt) - mean(y_block2) + grandmeanpop)^2 ))/(N-2)  #SS_inter/(N-2)
# 
# s<- sqrt(ss1 + ss2 + ss3) # S = sqrt(SS/(N-2))

### Added KC 1/15
# SS as sum of chi-squared distributions

ssw <- var_e*rchisq(n = 1, df = (2 * m * n - 2 * m))
ssb <- (var_e + 2*n*var_block)*rchisq(n = 1, df = (m - 1))
ssab <- (var_e + n*var_inter)*rchisq(n = 1, df = (m - 1))

s <- sqrt((ssw + ssb + ssab)/(N - 2))

### Added KC 1/15


#Determine d1 and d2
d2 <- delta/s

d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 

#Determine v1 and v2

v2_term_1 <- (2/(m*n))

v2_term_2 <- (d2^2)/(2 * (N - 2))

v2 <- v2_term_1  + v2_term_2

v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2

v1_term_2 <- v1_term_2_num/v1_term_2_denom

v1 <- v1_term_1 + v1_term_2

df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)

df_clean <- rbind(df_clean, df_clean_current)

names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")

}


mean(df_clean$Y)  
mean(df_clean$S)  

mean(df_clean$d1)  
mean(df_clean$d2) 

mean(df_clean$v1)  
mean(df_clean$v2)  

mean(df_clean$v1_1)  
mean(df_clean$v2_1)

mean(df_clean$v1_2)  
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  

(mean(df_clean$v2))/(var(df_clean$d2)) 

var(df_clean$d1)
var(df_clean$d2)

################################################
###### Artificial case var_e=0 & var_block=1

df_clean <- data.frame()

for (j in 1:10000){   
  m=2
  n=10
  N=2*m*n
  arms = 2
  var_e= 0
  var_inter = 0
  var_block = 1.0
  var_total = var_e + var_inter + var_block
  
  grandmean <- 0
  alpha_g1 <- (1/2)*0.4
  alpha_g2 <- (-1/2)*0.4
  
  
  rm_beta <- rnorm((m-1), 0, sqrt(var_block))
  
  last_beta<- -1*sum(rm_beta)
  
  beta<- c(rm_beta, last_beta)
  
  rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))
  
  interaction <- c(rm_interaction, -1*rm_interaction, -1*rm_interaction, rm_interaction)  
  
  noise <- rnorm( N , 0, sqrt(var_e))
  
  #rm_noise <- rnorm( (N - 2) , 0, sqrt(var_e))
  
  #last_noise <- -1*sum(rm_noise)
  
  #noise1<- c(rm_noise[1:199], last_noise)
  
  #noise2 <- c(rm_noise[200:398], -last_noise)
  
  #noise<- c(noise1, noise2)
  
  
  b1_up= N/(arms*m)
  b2_lw= b1_up + 1
  b2_up= b1_up + n
  b3_lw= b2_up + 1
  b3_up= b2_up + n 
  b4_lw= b3_up + 1
  b4_up= b3_up + n
  
  
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
  
  y_trt<- c(Y_trt1, Y_trt2)
  y_cnt<- c(Y_cnt1, Y_cnt2)
  
  y_block1<- c(Y_trt1, Y_cnt1)
  y_block2<- c(Y_trt2, Y_cnt2)
  
  var_trt<- var(y_trt)
  var_cnt<- var(y_cnt)
  
  grandmeanpop<- mean(Y)
  
  group<- c(rep(1, N/2), rep(2, N/2))
  
  block<- rep(seq(1, m), (N/m))
  
  #calculate delta
  delta <- mean(y_trt) - mean(y_cnt)
  
  #Calculate SS using formulas
  # ss1 <- ((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2)  #SS Error/(N-2)
  # 
  # ss2<- (n*arms*((mean(y_block1) - grandmeanpop)^2  + (mean(y_block2) - grandmeanpop)^2))/(N-2) #SS Block/(N-2)
  # 
  # #s<- sqrt(ss1 + ss2) # S = sqrt(SS/(N-2))
  # 
  # s<- sqrt(ss2) 
  
  ### Added KC 1/15
  # SS as sum of chi-squared distributions
  
  ssw <- rchisq(n = 1, df = (2 * m * n - 2 * m))
  ssb <- rchisq(n = 1, df = (m - 1))
  ssab <- rchisq(n = 1, df = (m - 1))
  
  s <- sqrt((ssw + ssb + ssab)/(N - 2))
  
  ### Added KC 1/15
  
  #Determine d1 and d2
  d2 <- delta/s
  
  d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  #d1 & d2 are both close to 0.4, but d1 is closer. Great news
  
  #Determine v1 and v2
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  #v1 tends to be smaller than v2, this is good news too
  
  df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)
  
  df_clean <- rbind(df_clean, df_clean_current)
  
  names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")
  
}

mean(df_clean$Y)  
mean(df_clean$S)  

mean(df_clean$d1)  
mean(df_clean$d2) 

mean(df_clean$v1)  
mean(df_clean$v2)  

mean(df_clean$v1_1)  
mean(df_clean$v2_1)

mean(df_clean$v1_2)  
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  

(mean(df_clean$v2))/(var(df_clean$d2)) 

var(df_clean$d1)
var(df_clean$d2)

###############################################
##########Case d=0############

#### sigma_e=1.0
df_clean <- data.frame()

for (j in 1:10000){   
  
  ##############################################################
  #Start with basic case where sigma_block & sigma_inter are zero
  m=2
  n=500
  N=2*m*n
  arms = 2
  var_e= 1.0
  var_inter = 0
  var_block = 0
  var_total = var_e + var_inter + var_block
  
  grandmean <- 0
  alpha_g1 <- (1/2)*0.4
  alpha_g2 <- (1/2)*0.4
  
  
  rm_beta <- rnorm((m-1), 0, sqrt(var_block))
  
  last_beta<- -1*sum(rm_beta)
  
  beta<- c(rm_beta, last_beta)
  
  rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))
  
  interaction <- c(rm_interaction, -1*rm_interaction, -1*rm_interaction, rm_interaction)  
  
  #rm_noise <- rnorm( (N - 2) , 0, sqrt(var_e))
  
  #last_noise <- -1*sum(rm_noise)
  
  #noise1<- c(rm_noise[1:199], last_noise)
  
  #noise2 <- c(rm_noise[200:398], -last_noise)
  
  #noise<- c(noise1, noise2)
  
  #noise<- c(rm_noise,-last_noise, last_noise)
  
  noise <- rnorm(N, 0, sqrt(var_e))
  
  b1_up= N/(arms*m)
  b2_lw= b1_up + 1
  b2_up= b1_up + n
  b3_lw= b2_up + 1
  b3_up= b2_up + n 
  b4_lw= b3_up + 1
  b4_up= b3_up + n
  
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
  
  y_trt<- c(Y_trt1, Y_trt2)
  y_cnt<- c(Y_cnt1, Y_cnt2)
  
  var_trt<- var(y_trt)
  var_cnt<- var(y_cnt)
  
  group<- c(rep(1, N/2), rep(2, N/2))
  
  block<- rep(seq(1, m), (N/m))
  
  #calculate delta
  delta <- mean(y_trt) - mean(y_cnt)
  
  #Estimate SS using formula
  # s<- sqrt(((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2))  #This is: S= sqrt(SS_error/(N-2))
  
  ### Added KC 1/15
  # SS as sum of chi-squared distributions
  
  ssw <- rchisq(n = 1, df = (2 * m * n - 2 * m))
  ssb <- rchisq(n = 1, df = (m - 1))
  ssab <- rchisq(n = 1, df = (m - 1))
  
  s <- sqrt((ssw + ssb + ssab)/(N - 2))
  
  ### Added KC 1/15
  
  #this matches what we get from summary(lm1)$sigma, so we are good
  
  #Determine d1 and d2
  d2 <- delta/s
  
  d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  #g2 <- d2 * (1 - (3/(4*(N-2) -1)))
  
  #g1 <- d1 * (1 - (3/(4*(N-2) -1)))
  
  #d1 & d2 match and are "close" to 0.4
  
  #Determine v1 and v2
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  #v1 & v2 match
  
  #All good so far
  
  df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)
  
  df_clean <- rbind(df_clean, df_clean_current)
  
  names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")
  
}

mean(df_clean$Y)  #0.997
mean(df_clean$S)  #0.9949

mean(df_clean$d1)  #0.408
mean(df_clean$d2) #0.408

mean(df_clean$v1)  #0.103
mean(df_clean$v2)  #0.103

mean(df_clean$v1_1)  #0.103
mean(df_clean$v2_1)

mean(df_clean$v1_2)  #0.103
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  #0.07 (actual variance)

var(df_clean$d1)
var(df_clean$d2)


#######Case####################
#Move on to the second case where sigma_block=0.7 & sigma_inter=0
df_clean <- data.frame()

for (j in 1:10000){   
  m=2
  n=500
  N=2*m*n
  arms = 2
  var_e= 0.3
  var_inter = 0
  var_block = 0.7
  var_total = var_e + var_inter + var_block
  
  grandmean <- 0
  alpha_g1 <- (1/2)*0.4
  alpha_g2 <- (1/2)*0.4
  
  
  rm_beta <- rnorm((m-1), 0, sqrt(var_block))
  
  last_beta<- -1*sum(rm_beta)
  
  beta<- c(rm_beta, last_beta)
  
  rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))
  
  interaction <- c(rm_interaction, -1*rm_interaction, -1*rm_interaction, rm_interaction)  
  
  noise <- rnorm( N , 0, sqrt(var_e))
  
  #rm_noise <- rnorm( (N - 2) , 0, sqrt(var_e))
  
  #last_noise <- -1*sum(rm_noise)
  
  #noise1<- c(rm_noise[1:199], last_noise)
  
  #noise2 <- c(rm_noise[200:398], -last_noise)
  
  #noise<- c(noise1, noise2)
  
  
  b1_up= N/(arms*m)
  b2_lw= b1_up + 1
  b2_up= b1_up + n
  b3_lw= b2_up + 1
  b3_up= b2_up + n 
  b4_lw= b3_up + 1
  b4_up= b3_up + n
  
  
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
  
  y_trt<- c(Y_trt1, Y_trt2)
  y_cnt<- c(Y_cnt1, Y_cnt2)
  
  y_block1<- c(Y_trt1, Y_cnt1)
  y_block2<- c(Y_trt2, Y_cnt2)
  
  var_trt<- var(y_trt)
  var_cnt<- var(y_cnt)
  
  grandmeanpop<- mean(Y)
  
  group<- c(rep(1, N/2), rep(2, N/2))
  
  block<- rep(seq(1, m), (N/m))
  
  #calculate delta
  delta <- mean(y_trt) - mean(y_cnt)
  
  #Calculate SS using formulas
  # ss1 <- ((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2)  #SS Error/(N-2)
  # 
  # ss2<- (n*arms*((mean(y_block1) - grandmeanpop)^2  + (mean(y_block2) - grandmeanpop)^2))/(N-2) #SS Block/(N-2)
  # 
  # s<- sqrt(ss1 + ss2) # S = sqrt(SS/(N-2))
  
  ### Added KC 1/15
  # SS as sum of chi-squared distributions
  
  ssw <- rchisq(n = 1, df = (2 * m * n - 2 * m))
  ssb <- rchisq(n = 1, df = (m - 1))
  ssab <- rchisq(n = 1, df = (m - 1))
  
  s <- sqrt((ssw + ssb + ssab)/(N - 2))
  
  ### Added KC 1/15
  
  #Determine d1 and d2
  d2 <- delta/s
  
  d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  #d1 & d2 are both close to 0.4, but d1 is closer. Great news
  
  #Determine v1 and v2
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  #v1 tends to be smaller than v2, this is good news too
  
  df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)
  
  df_clean <- rbind(df_clean, df_clean_current)
  
  names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")
  
}

mean(df_clean$Y)  
mean(df_clean$S)  

mean(df_clean$d1)  
mean(df_clean$d2) 

mean(df_clean$v1)  
mean(df_clean$v2)  

mean(df_clean$v1_1)  
mean(df_clean$v2_1)

mean(df_clean$v1_2)  
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  

(mean(df_clean$v2))/(var(df_clean$d2)) 

var(df_clean$d1)
var(df_clean$d2)


###########################################################
#Move on to the third case where sigma_block=0.5 & sigma_inter=0.2 & sigma_e=0.3

df_clean <- data.frame()

for (j in 1:10000){ 
  m=2
  n=500
  N=2*m*n
  arms = 2
  var_e=0.3
  var_inter = 0.2
  var_block = 0.5
  var_total = var_e + var_inter + var_block
  
  grandmean <- 0
  alpha_g1 <- (1/2)*0.4
  alpha_g2 <- (1/2)*0.4
  
  
  rm_beta <- rnorm((m-1), 0, sqrt((var_block)))
  
  last_beta<- -1*sum(rm_beta)
  
  beta<- c(rm_beta, last_beta)
  
  rm_interaction <- rnorm( (m-1) , 0, sqrt(var_inter))
  
  interaction <- c(rm_interaction, -1*rm_interaction,  rm_interaction, -1*rm_interaction)  
  
  noise <- rnorm( N , 0, sqrt(var_e))
  
  b1_up= N/(arms*m)
  b2_lw= b1_up + 1
  b2_up= b1_up + n
  b3_lw= b2_up + 1
  b3_up= b2_up + n 
  b4_lw= b3_up + 1
  b4_up= b3_up + n
  
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
  
  y_trt<- c(Y_trt1, Y_trt2)
  y_cnt<- c(Y_cnt1, Y_cnt2)
  
  y_block1<- c(Y_trt1, Y_cnt1)
  y_block2<- c(Y_trt2, Y_cnt2)
  
  var_trt<- var(y_trt)
  var_cnt<- var(y_cnt)
  
  grandmeanpop<- mean(Y)
  
  group<- c(rep(1, N/2), rep(2, N/2))
  
  block<- rep(seq(1, m), (N/m))
  
  #calculate delta
  delta <- mean(y_trt) - mean(y_cnt)
  
  # ss1 <- ((n*m - 1)*var_trt+(n*m-1)*var_cnt)/(N-2)  #SS_Error/(N-2)
  # 
  # ss2<- (n*arms*((mean(y_block1) - grandmeanpop)^2  + (mean(y_block2) - grandmeanpop)^2))/(N-2) #SS_Block/(N-2)
  # 
  # ss3<- (n*(  (mean(Y_trt1) - mean(y_trt) - mean(y_block1) + grandmeanpop)^2 
  #             + (mean(Y_trt2)   - mean(y_trt) - mean(y_block2) + grandmeanpop)^2 
  #             + (mean(Y_cnt1) - mean(y_cnt) - mean(y_block1) + grandmeanpop)^2 
  #             + (mean(Y_cnt2) - mean(y_cnt) - mean(y_block2) + grandmeanpop)^2 ))/(N-2)  #SS_inter/(N-2)
  # 
  # s<- sqrt(ss1 + ss2 + ss3) # S = sqrt(SS/(N-2))
  
  ### Added KC 1/15
  # SS as sum of chi-squared distributions
  
  ssw <- rchisq(n = 1, df = (2 * m * n - 2 * m))
  ssb <- rchisq(n = 1, df = (m - 1))
  ssab <- rchisq(n = 1, df = (m - 1))
  
  s <- sqrt((ssw + ssb + ssab)/(N - 2))
  
  ### Added KC 1/15
  
  #Determine d1 and d2
  d2 <- delta/s
  
  d1<- (delta/s)*sqrt(((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))/((N - 2) *var_total)) 
  
  #Determine v1 and v2
  
  v2_term_1 <- (2/(m*n))
  
  v2_term_2 <- (d2^2)/(2 * (N - 2))
  
  v2 <- v2_term_1  + v2_term_2
  
  v1_term_1 <- (2/(m * n)) * ((var_e + n * var_inter)/var_total)
  v1_term_2_num <- ((N - 2 * m) * var_e^2 + (m - 1) * (var_e + 2 * n * var_block)^2 + (m - 1) * (var_e + n * var_inter)^2) * d1^2
  v1_term_2_denom <- 2 * ((N - 2 * m) * var_e + (m - 1) * (var_e + 2 * n * var_block) + (m - 1) * (var_e + n * var_inter))^2
  
  v1_term_2 <- v1_term_2_num/v1_term_2_denom
  
  v1 <- v1_term_1 + v1_term_2
  
  df_clean_current<- c(delta,s , d1, d2,v1, v2, v1_term_1 , v1_term_2 , v2_term_1, v2_term_2)
  
  df_clean <- rbind(df_clean, df_clean_current)
  
  names(df_clean) <- c("Y", "S","d1", "d2", "v1", "v2", "v1_1", "v1_2", "v2_1", "v2_2")
  
}


mean(df_clean$Y)  
mean(df_clean$S)  

mean(df_clean$d1)  
mean(df_clean$d2) 

mean(df_clean$v1)  
mean(df_clean$v2)  

mean(df_clean$v1_1)  
mean(df_clean$v2_1)

mean(df_clean$v1_2)  
mean(df_clean$v2_2)

(mean(df_clean$v1))/(var(df_clean$d1))  

(mean(df_clean$v2))/(var(df_clean$d2)) 

var(df_clean$d1)
var(df_clean$d2)


#############################################################
