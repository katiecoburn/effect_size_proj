library(MASS)
library(tidyverse)
library(mvtnorm)

grandmean=0

alpha= 0.5

#Case with m=2
i=2
j=2
k=10
var_error= 0.5
var_block=0.5
var_inter=0
var_total=1

beta= rnorm(j, 0, var_block)
interaction= rnorm(i*j, 0, var_inter)
error = rnorm (i*j*k, 0, var_error)

Y_trt_1= c()
Y_trt_2= c()
Y_cnt_1=c()
Y_cnt_2 = c()

#note that i should go to 40
for(i in 1:10) {
trt_1 = grandmean + alpha[1] + beta[1] + interaction[1] + error[i]
Y_trt_1 = rbind(Y_trt_1, trt_1)
trt_2 = grandmean + alpha[1] + beta[2] + interaction[2] + error[i]
Y_trt_2 = rbind(Y_trt_2, trt_2)
cnt_1 = grandmean + alpha[2] + beta[1] + interaction[3] + error[i]
Y_cnt_1 = rbind(Y_cnt_1, cnt_1)
cnt_2 = grandmean + alpha[2] + beta[2] + interaction[4] + error[i]
Y_cnt_2 = rbind(Y_cnt_2, cnt_2)
}


data<- c(Y_trt_1, Y_trt_2, Y_cnt_1, Y_cnt_2)  
  
data2<- as.data.frame(data)

data2$group<- c(rep("T", 20), rep("C", 20))

data2$block<- c( rep(1, 10), rep(2, 10), rep(1, 10), rep(2, 10))

names(data2)<- c("outcome", "Group", "Block")

model1<- aov(outcome ~ Group + Block + Group:Block, data=data2)

summary(model1)

model2<- lm(outcome ~ Group + Block + Group:Block, data=data2)
summary(model2) 


# #####################################################
# ######Covariance Matrix second option###############
# set.seed=1234
# 
# sigma<- matrix( c(0.5, 0.6, 0.9, 0.9, 
#                   0.6, 0.5, 0.9, 0.9,
#                   0.9, 0.9, 0.5, 0.6,
#                   0.9, 0.9, 0.6, 0.5), ncol=4)
# 
# 
# y<-rmvnorm(n=10, mean = c(1.1, 0.5, 0.2, 0.6), sigma=sigma)
# 
# data<- as.data.frame(y)
# 
# data2<- c(data$V1, data$V2, data$V3, data$V4)
# 
# data2<- as.data.frame(data2)
# 
# data2$group<- c(rep("T", 20), rep("C", 20))
# 
# data2$block<- c( rep(1, 10), rep(2, 10), rep(1, 10), rep(2, 10))
# 
# names(data2)<- c("outcome", "Group", "Block")
# 
# model1<- aov(outcome ~ Group + Block + Group:Block, data=data2)
# 
# summary(model1)
# 
# model2<- lm(outcome ~ Group + Block + Group:Block, data=data2)
# summary(model2)

