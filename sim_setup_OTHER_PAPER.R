library(MASS)
library(tidyverse)

# Create empty data frames for gains and ANCOVA
df_clean_rho_cumulative_gains <- data.frame()
df_clean_rho_cumulative_ancova <- data.frame()

# Vector of n values
n_vector <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

for (j in 1:length(n_vector)){    
    print(j)
    df_clean_rho_gains <- data.frame()
    df_clean_rho_ancova <- data.frame()
    
    for (i in 1:100000){
        #-------------------------------------------------------------
        rho <- 0.4
        theta <- 0.8
        q <- 1
        iterations <- 100000
        
    
        #X & Y Treatment
        n_trt <- n_vector[j]
        mu_trt <-  c(0,0.9)
        sigma <- matrix(c(1,rho,rho,1), ncol=2)
        # note that first column is x_trt, second is y_trt
        
        ##--#--#---#--- Next line is random sampling:
        A <- mvrnorm(n=n_trt, mu_trt, sigma)
        
        #Store gain in treatment group
        gain_trt<-(A[,2]- A[,1])
        within_trt_data <- cbind(A, gain_trt) 
        colnames(within_trt_data) <- c("x_trt", "y_trt", "gain_trt")
        Ancova_sims <- A %>% as_tibble() %>%  mutate(treatment_indicator = rep(1, n_trt))
        colnames(Ancova_sims) <- c("x", "y", "treatment_indicator")
        
        #-------------------------------------------------------------
        #X & Y Control
        n_cnt <- n_vector[j]
        mu_cnt <- c(0,0.1)
        sigma <- matrix(c(1,rho,rho,1), ncol=2)
        
        ##--#--#---#--- Next line is random sampling:
        B <- mvrnorm(n=n_cnt, mu_cnt, sigma)
        
        #Store gain in control group
        gain_cnt <- (B[,2]- B[,1])
        within_cnt_data <- cbind(B, gain_cnt)
        colnames(within_cnt_data) <- c("x_cnt", "y_cnt", "gain_cnt")
        Bncova_sims <- B %>% as_tibble() %>%  mutate(treatment_indicator = rep(0, n_cnt))
        colnames(Bncova_sims) <- c("x", "y", "treatment_indicator")
        
        #-------------------------------------------------------------
        #Store Treatment and Control variables in one dataframe
        trtANDcnt_df<- cbind(within_trt_data, within_cnt_data) %>% as.data.frame() 
        #Store variables for ANCOVA dataframe
        ancova_df <- rbind(Ancova_sims,Bncova_sims) %>% mutate(treatment_indicator = as.character(treatment_indicator),
                                           x = scale(x, center = FALSE, scale = FALSE),
                                           y = scale(y, center = FALSE, scale = FALSE))
        
        #--------------------------------------------------------------
        # Running the ANCOVA model
        ancova_mod <- lm(ancova_df$y ~ as.factor(ancova_df$treatment_indicator) + ancova_df$x)
        anc_mod_summary <- summary(ancova_mod)
        # Extract SD adjusted
        sd_adjusted <- anc_mod_summary$sigma
        # Extract the difference in adjusted means
        diffMeansAdjusted <- anc_mod_summary$coefficients[2,1]
        df_complete_ancova <- cbind(ancova_df, sd_adjusted) 
        colnames(df_complete_ancova) <- c("x", "y", "treatment_indicator", "sd_adjusted")
        #-------------------------------------------------------------
        
        # long dataframe for GAINS
        df_gains <- trtANDcnt_df %>% 
            mutate(
                   ###################################################
                   diff_means_post = mean(y_trt-y_cnt),
                   sd_post = sqrt((var(y_trt)+var(y_cnt))/2),
                   # d
                   t_1 = (diff_means_post/sd_post),
                   unbiased_t1 = (1-(3/(4*(n_trt+n_cnt-2)-1)))*(diff_means_post/sd_post),
                   v_1 = ((1/((n_trt*n_cnt)/(n_trt+n_cnt)))+(t_1^2/(2*(n_trt+n_cnt-2)))),
                   unbiased_v1 = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*
                       (((1/((n_trt*n_cnt)/(n_trt+n_cnt)))+(t_1^2/(2*(n_trt+n_cnt-2))))),
                   diff_gains = gain_trt - gain_cnt,
                   mean_diff_gains = mean(diff_gains),
                   var_gain_trt = var(gain_trt),
                   var_gain_cnt = var(gain_cnt),
                   sd_pooled = sqrt(((n_trt - 1)*var_gain_trt+(n_cnt-1)*var_gain_cnt)/(n_trt+n_cnt-2)),
                   cov_trt = cov(x_trt, y_trt),
                   cov_cnt = cov(x_cnt, y_cnt),
                   LB_t1 = t_1-1.96*sqrt(v_1),
                   UB_t1 = t_1+1.96*sqrt(v_1),
                   in_CI_t1 = ifelse((theta >= LB_t1 & theta <= UB_t1), 1, 0),
                   
                   combo_LB_t1 = unbiased_t1-1.96*sqrt(v_1),
                   combo_UB_t1 = unbiased_t1+1.96*sqrt(v_1),
                   combo_in_CI_t1 = ifelse((theta >= combo_LB_t1 & theta <= combo_UB_t1), 1, 0),
                   
                   LB_unbiased_t1 = unbiased_t1-1.96*sqrt(unbiased_v1),
                   UB_unbiased_t1 = unbiased_t1+1.96*sqrt(unbiased_v1),
                   in_CI_unbiased_t1 = ifelse((theta >= LB_unbiased_t1 & theta <= UB_unbiased_t1), 1, 0),
        
                   # note we have the assumption of a common covariance matrix
                   corr_pooled_trt = cov_trt/(sqrt(var(x_trt))*sqrt(var(y_trt))),
                   corr_pooled_cnt = cov_cnt/(sqrt(var(x_cnt))*sqrt(var(y_cnt))),
                   corr_pooled = (corr_pooled_trt+corr_pooled_cnt)/2,
                   part1_in_estT5 = mean_diff_gains/sd_pooled,
                   part2_in_estT5 = sqrt(2*(1-corr_pooled)),
                   # dG1
                   t_2 = (mean_diff_gains/sd_post),
                   unbiased_t2 = (1-(3/(4*(n_trt+n_cnt-2)-1)))*(mean_diff_gains/sd_post),
                   v_2 =  (((2*(1-rho))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_2^2)/(2*(n_trt+n_cnt-2)))),
                   unbiased_v2 = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*
                       ((((2*(1-rho))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_2^2)/(2*(n_trt+n_cnt-2))))),
                   LB_t2 = t_2-1.96*sqrt(v_2),
                   UB_t2 = t_2+1.96*sqrt(v_2),
                   in_CI_t2 = ifelse((theta >= LB_t2 & theta <= UB_t2), 1, 0),
                   
                   combo_LB_t2 = unbiased_t2-1.96*sqrt(v_2),
                   combo_UB_t2 = unbiased_t2+1.96*sqrt(v_2),
                   combo_in_CI_t2 = ifelse((theta >= combo_LB_t2 & theta <= combo_UB_t2), 1, 0),
                   
                   LB_unbiased_t2 = unbiased_t2-1.96*sqrt(unbiased_v2),
                   UB_unbiased_t2 = unbiased_t2+1.96*sqrt(unbiased_v2),
                   in_CI_unbiased_t2 = ifelse((theta >= LB_unbiased_t2 & theta <= UB_unbiased_t2), 1, 0),
                   
                   # dG2
                   t_5_trueRho = (part1_in_estT5*sqrt(2*(1-rho))),
                   unbiased_t5_trueRho = (1-(3/(4*(n_trt+n_cnt-2)-1)))*(part1_in_estT5*sqrt(2*(1-rho))),
                   v_5_trueRho = (((2*(1-rho))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_5_trueRho^2)/(2*(n_trt+n_cnt-2)))),
                   unbiased_v5_trueRho = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*
                       ((((2*(1-rho))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_5_trueRho^2)/(2*(n_trt+n_cnt-2))))),
                   LB_t5_trueRho = t_5_trueRho-1.96*sqrt(v_5_trueRho),
                   UB_t5_trueRho = t_5_trueRho+1.96*sqrt(v_5_trueRho),
                   in_CI_t5_trueRho = ifelse((theta >= LB_t5_trueRho & theta <= UB_t5_trueRho), 1, 0),
                   
                   combo_LB_t5_trueRho = unbiased_t5_trueRho-1.96*sqrt(v_5_trueRho),
                   combo_UB_t5_trueRho = unbiased_t5_trueRho+1.96*sqrt(v_5_trueRho),
                   combo_in_CI_t5_trueRho = ifelse((theta >= combo_LB_t5_trueRho & theta <= combo_UB_t5_trueRho), 1, 0),
                   
                   LB_unbiased_t5_trueRho = unbiased_t5_trueRho-1.96*sqrt(unbiased_v5_trueRho),
                   UB_unbiased_t5_trueRho = unbiased_t5_trueRho+1.96*sqrt(unbiased_v5_trueRho),
                   in_CI_unbiased_t5_trueRho = ifelse((theta >= LB_unbiased_t5_trueRho & theta <= UB_unbiased_t5_trueRho), 1, 0),
                   
                   # dG3
                   t_5 = (part1_in_estT5*part2_in_estT5),
                   unbiased_t5 = (1-(3/(4*(n_trt+n_cnt-2)-1)))*(part1_in_estT5*part2_in_estT5),
                   v_5 = ((2*(1-corr_pooled))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_5^2)/(2*(n_trt+n_cnt-2)))+
                                   (((t_5^2)*(1-corr_pooled^2)^2)/(2*(n_trt+n_cnt)*(1-corr_pooled))),
                   unbiased_v5 = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*(((2*(1-corr_pooled))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_5^2)/(2*(n_trt+n_cnt-2)))+
                       (((t_5^2)*(1-corr_pooled^2)^2)/(2*(n_trt+n_cnt)*(1-corr_pooled)))),
                   LB_t5 = t_5-1.96*sqrt(v_5),
                   UB_t5 = t_5+1.96*sqrt(v_5),
                   in_CI_t5 = ifelse((theta >= LB_t5 & theta <= UB_t5), 1, 0),
                   
                   combo_LB_t5 = unbiased_t5-1.96*sqrt(v_5),
                   combo_UB_t5 = unbiased_t5+1.96*sqrt(v_5),
                   combo_in_CI_t5 = ifelse((theta >= combo_LB_t5 & theta <= combo_UB_t5), 1, 0),
                   
                   LB_unbiased_t5 = unbiased_t5-1.96*sqrt(unbiased_v5),
                   UB_unbiased_t5 = unbiased_t5+1.96*sqrt(unbiased_v5),
                   in_CI_unbiased_t5 = ifelse((theta >= LB_unbiased_t5 & theta <= UB_unbiased_t5), 1, 0),
                   
                   #### Pooled across trt and ctrl across pre and post
                   sPOOLED = sqrt((var(y_trt)+var(y_cnt)+var(x_trt)+var(x_cnt))/4),
                   # dP
                   t_10 = diff_means_post/sPOOLED,
                   unbiased_t10 = (1-(3/(4*(n_trt+n_cnt-2)-1)))*diff_means_post/sPOOLED,
                   v_10 = (1/((n_trt*n_cnt)/(n_trt+n_cnt)))+((1+rho^2)*(t_10^2)/(4*(n_trt+n_cnt-2))),
                   unbiased_v10 = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*
                       ((1/((n_trt*n_cnt)/(n_trt+n_cnt)))+((1+rho^2)*(t_10^2)/(4*(n_trt+n_cnt-2)))),
                              
                   LB_t10 = t_10-1.96*sqrt(v_10),
                   UB_t10 = t_10+1.96*sqrt(v_10),
                   in_CI_t10 = ifelse((theta >= LB_t10 & theta <= UB_t10), 1, 0),
                   
                   combo_LB_t10 = unbiased_t10-1.96*sqrt(v_10),
                   combo_UB_t10 = unbiased_t10+1.96*sqrt(v_10),
                   combo_in_CI_t10 = ifelse((theta >= combo_LB_t10 & theta <= combo_UB_t10), 1, 0),
                   
                   LB_unbiased_t10 = unbiased_t10-1.96*sqrt(unbiased_v10),
                   UB_unbiased_t10 = unbiased_t10+1.96*sqrt(unbiased_v10),
                   in_CI_unbiased_t10 = ifelse((theta >= LB_unbiased_t10 & theta <= UB_unbiased_t10), 1, 0),
                   
                   # dPG
                   t_11 = mean_diff_gains/sPOOLED,
                   unbiased_t11 = (1-(3/(4*(n_trt+n_cnt-2)-1)))*(mean_diff_gains/sPOOLED),
                   v_11 = ((2*(1-rho))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((1+rho^2)*(t_11^2)/(4*(n_trt+n_cnt-2))),
                   unbiased_v11 = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*
                       (((2*(1-rho))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((1+rho^2)*(t_11^2)/(4*(n_trt+n_cnt-2)))),
                           
                   LB_t11 = t_11-1.96*sqrt(v_11),
                   UB_t11 = t_11+1.96*sqrt(v_11),
                   in_CI_t11 = ifelse((theta >= LB_t11 & theta <= UB_t11), 1, 0),
                   
                   combo_LB_t11 = unbiased_t11-1.96*sqrt(v_11),
                   combo_UB_t11 = unbiased_t11+1.96*sqrt(v_11),
                   combo_in_CI_t11 = ifelse((theta >= combo_LB_t11 & theta <= combo_UB_t11), 1, 0),
                   
                   LB_unbiased_t11 = unbiased_t11-1.96*sqrt(unbiased_v11),
                   UB_unbiased_t11 = unbiased_t11+1.96*sqrt(unbiased_v11),
                   in_CI_unbiased_t11 = ifelse((theta >= LB_unbiased_t11 & theta <= UB_unbiased_t11), 1, 0)
            )
        
        # long dataframe for ANCOVA
        df_adjusted <- df_complete_ancova %>% 
            group_by(treatment_indicator) %>% 
            mutate(corr_sample = cor(x,y),
                   var_post_group = var(y),
                   var_pre_group = var(x)) %>% 
            ungroup() %>% 
            mutate(
                sd_post_pooled = sqrt((var_post_group[1]+var_post_group[length(var_post_group)])/2),
                ###############
                # T3 
                t_3 = diffMeansAdjusted/sd_post_pooled,
                unbiased_t3 = (1-(3/(4*(n_trt+n_cnt-2-q)-1)))*diffMeansAdjusted/sd_post_pooled,
                v_3 = ((1-rho^2)/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_3^2)/(2*(n_trt+n_cnt-2))),
                unbiased_v3 = (1-(3/(4*(n_trt+n_cnt-2-q)-1)))^2*
                    (((1-rho^2)/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_3^2)/(2*(n_trt+n_cnt-2)))),
                LB_t3 = t_3-1.96*sqrt(v_3),
                UB_t3 = t_3+1.96*sqrt(v_3),
                in_CI_t3 = ifelse((theta >= LB_t3 & theta <= UB_t3), 1, 0),
                
                combo_LB_t3 = unbiased_t3-1.96*sqrt(v_3),
                combo_UB_t3 = unbiased_t3+1.96*sqrt(v_3),
                combo_in_CI_t3 = ifelse((theta >= combo_LB_t3 & theta <= combo_UB_t3), 1, 0),
                
                LB_unbiased_t3 = unbiased_t3-1.96*sqrt(unbiased_v3),
                UB_unbiased_t3 = unbiased_t3+1.96*sqrt(unbiased_v3),
                in_CI_unbiased_t3 = ifelse((theta >= LB_unbiased_t3 & theta <= UB_unbiased_t3), 1, 0),
                
                
                ###############
                # T9 with estimated rho
                part1_t9 = diffMeansAdjusted/sd_adjusted[1],
                corr_pooled = (corr_sample[1]+corr_sample[length(corr_sample)])/2,
                part2_t9 = sqrt(1-corr_pooled^2),
                t_9 = part1_t9*part2_t9,
                unbiased_t9 = (1-(3/(4*(n_trt+n_cnt-2-q)-1)))*part1_t9*part2_t9,
                v_9 = ((1-corr_pooled^2)/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_9^2)/(2*(n_trt+n_cnt-2-q)))+
                    ((t_9^2)*(corr_pooled^2)*(1-corr_pooled^2))/(n_trt+n_cnt),
                unbiased_v9 = (1-(3/(4*(n_trt+n_cnt-2-q)-1)))^2*
                    (((1-corr_pooled^2)/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_9^2)/(2*(n_trt+n_cnt-2-q)))+
                         ((t_9^2)*(corr_pooled^2)*(1-corr_pooled^2))/(n_trt+n_cnt)),
                LB_t9 = t_9-1.96*sqrt(v_9),
                UB_t9 = t_9+1.96*sqrt(v_9),
                in_CI_t9 = ifelse((theta >= LB_t9 & theta <= UB_t9), 1, 0),
                
                combo_LB_t9 = unbiased_t9-1.96*sqrt(v_9),
                combo_UB_t9 = unbiased_t9+1.96*sqrt(v_9),
                combo_in_CI_t9 = ifelse((theta >= combo_LB_t9 & theta <= combo_UB_t9), 1, 0),
                
                LB_unbiased_t9 = unbiased_t9-1.96*sqrt(unbiased_v9),
                UB_unbiased_t9 = unbiased_t9+1.96*sqrt(unbiased_v9),
                in_CI_unbiased_t9 = ifelse((theta >= LB_unbiased_t9 & theta <= UB_unbiased_t9), 1, 0),
                
                ##############
                # T9 with TRUE rho
                part1_t9_trueRho = diffMeansAdjusted/sd_adjusted[1],
                part2_t9_trueRho = sqrt(1-rho^2),
                t_9_trueRho = part1_t9_trueRho*part2_t9_trueRho,
                unbiased_t9_trueRho = (1-(3/(4*(n_trt+n_cnt-2-q)-1)))*part1_t9_trueRho*part2_t9_trueRho,
                v_9_trueRho = ((1-rho^2)/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_9_trueRho^2)/(2*(n_trt+n_cnt-2-q))),
                unbiased_v9_trueRho = (1-(3/(4*(n_trt+n_cnt-2-q)-1)))^2*
                    (((1-rho^2)/((n_trt*n_cnt)/(n_trt+n_cnt)))+((t_9_trueRho^2)/(2*(n_trt+n_cnt-2-q)))),
                LB_t9_trueRho = t_9_trueRho-1.96*sqrt(v_9_trueRho),
                UB_t9_trueRho = t_9_trueRho+1.96*sqrt(v_9_trueRho),
                in_CI_t9_trueRho = ifelse((theta >= LB_t9_trueRho & theta <= UB_t9_trueRho), 1, 0),
                
                combo_LB_t9_trueRho = unbiased_t9_trueRho-1.96*sqrt(v_9_trueRho),
                combo_UB_t9_trueRho = unbiased_t9_trueRho+1.96*sqrt(v_9_trueRho),
                combo_in_CI_t9_trueRho = ifelse((theta >= combo_LB_t9_trueRho & theta <= combo_UB_t9_trueRho), 1, 0),
                
                LB_unbiased_t9_trueRho = unbiased_t9_trueRho-1.96*sqrt(unbiased_v9_trueRho),
                UB_unbiased_t9_trueRho = unbiased_t9_trueRho+1.96*sqrt(unbiased_v9_trueRho),
                in_CI_unbiased_t9_trueRho = ifelse((theta >= LB_unbiased_t9_trueRho & theta <= UB_unbiased_t9_trueRho), 1, 0),
                
                # dPA
                sPOOLED_anc = sqrt((var_post_group[1]+var_post_group[length(var_post_group)]+
                                     var_pre_group[1]+var_pre_group[length(var_pre_group)])/4),
                
                t_12 = diffMeansAdjusted/sPOOLED_anc,
                unbiased_t12 = (1-(3/(4*(n_trt+n_cnt-2)-1)))*diffMeansAdjusted/sPOOLED_anc,
                v_12 = (((1-rho^2))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((1+rho^2)*(t_12^2)/(4*(n_trt+n_cnt-2))),
                unbiased_v12 = (1-(3/(4*(n_trt+n_cnt-2)-1)))^2*
                    ((((1-rho^2))/((n_trt*n_cnt)/(n_trt+n_cnt)))+((1+rho^2)*(t_12^2)/(4*(n_trt+n_cnt-2)))),
                LB_t12 = t_12-1.96*sqrt(v_12),
                UB_t12 = t_12+1.96*sqrt(v_12),
                in_CI_t12 = ifelse((theta >= LB_t12 & theta <= UB_t12), 1, 0),
                
                combo_LB_t12 = unbiased_t12-1.96*sqrt(v_12),
                combo_UB_t12 = unbiased_t12+1.96*sqrt(v_12),
                combo_in_CI_t12 = ifelse((theta >= combo_LB_t12 & theta <= combo_UB_t12), 1, 0),
                
                LB_unbiased_t12 = unbiased_t12-1.96*sqrt(unbiased_v12),
                UB_unbiased_t12 = unbiased_t12+1.96*sqrt(unbiased_v12),
                in_CI_unbiased_t12 = ifelse((theta >= LB_unbiased_t12 & theta <= UB_unbiased_t12), 1, 0)
            )
        
        #-----------------------------------------------------------------------------------------------
        # remove x_trt, y_trt, x_cnt, y_cnt, and diff_gains and clean dataframe ANCOVA
        
        df_clean_current_ancova <- df_adjusted %>% select(-c(x, y, treatment_indicator, 
                                             corr_sample)) %>% 
            distinct()
        #-----------------------------------------------------------------------------------------------
        
        # remove x_trt, y_trt, x_cnt, y_cnt, and diff_gains and clean dataframe GAINS
        
        df_clean_current_gains <- df_gains %>% select(-c(x_trt, y_trt, x_cnt, y_cnt, diff_gains, gain_trt, gain_cnt)) %>% 
            distinct()
        #-----------------------------------------------------------------------------------------------
        # final data frame GAINS
        df_clean_rho_gains <- rbind(df_clean_rho_gains, df_clean_current_gains)
        
        # final data frame ANCOVA
        df_clean_rho_ancova <- rbind(df_clean_rho_ancova, df_clean_current_ancova)
    }    

    ##########################
    # GAINS
    df_rho_gains <- df_clean_rho_gains %>% 
        mutate(n_trt = n_trt,
               n_cnt = n_cnt,
               #########
               # t5 properties
               mean_t5 = mean(t_5),
               var_t5 = var(t_5),
               bias_t5 = mean(t_5) - theta,
               MCse_bias_t5 = sqrt(var_t5/iterations),
               rel_bias_t5 = bias_t5/theta,
               mse_t5 = bias_t5^2+var_t5,
               mse_t5_true = (t_5-mean_t5)^2,
               mse_t5true_avg = (1/iterations)*sum(mse_t5_true),
               mean_v5 = mean(v_5),
               var_v5 = var(v_5),
               se_t5 = sd(t_5)/(sqrt(iterations)),
               se_empVar_t5 = sqrt(var_t5)/sqrt(2*(iterations-1)),
               ratio_var_t5 = mean_v5/var_t5,
               MC_SE_model_se_t5 = sqrt(var_v5/(4*iterations*sqrt(mean_v5))),
               cov_prob_t5 = sum(in_CI_t5)/iterations,
               # t5 unbiased properties
               mean_unbiased_t5 = mean(unbiased_t5),
               var_unbiased_t5 = var(unbiased_t5),
               bias_unbiased_t5 = mean_unbiased_t5 - theta,
               MCse_bias_unbiased_t5 = sqrt(var_unbiased_t5/iterations),
               rel_bias_unbiased_t5 = bias_unbiased_t5/theta,
               mse_unbiased_t5 = bias_unbiased_t5^2+var_unbiased_t5,
               mse_unbiased_t5_true = (unbiased_t5-mean_unbiased_t5)^2,
               mse_unbiased_t5true_avg = (1/iterations)*sum(mse_unbiased_t5_true),
               mean_unbiased_v5 = mean(unbiased_v5),
               var_unbiased_v5 = var(unbiased_v5),
               se_unbiased_t5 = sd(unbiased_t5)/(sqrt(iterations)),
               se_empVar_unbiased_t5 = sqrt(var_unbiased_t5)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t5 = sqrt(var_unbiased_v5/(4*iterations*sqrt(mean_unbiased_v5))),
               cov_prob_unbiased_t5 = sum(in_CI_unbiased_t5)/iterations,
               combo_cov_prob_unbiased_t5 = sum(combo_in_CI_t5)/iterations,
               
               
               #########
               # t5 with true Rho properties
               mean_t5_trueRho = mean(t_5_trueRho),
               mean_unbiased_t5_trueRho = mean(unbiased_t5_trueRho),
               
               var_t5_trueRho = var(t_5_trueRho),
               bias_t5_trueRho = mean(t_5_trueRho) - theta,
               MCse_bias_t5_trueRho = sqrt(var_t5_trueRho/iterations),
               
               rel_bias_t5_trueRho = bias_t5_trueRho/theta,
               mse_t5_trueRho = bias_t5_trueRho^2+var_t5_trueRho,
               mse_t5_true_trueRho = (t_5_trueRho-mean_t5_trueRho)^2,
               mse_t5true_avg_trueRho = (1/iterations)*sum(mse_t5_true_trueRho),
               mean_v5_trueRho = mean(v_5_trueRho),
               var_v5_trueRho = var(v_5_trueRho),
               se_t5_trueRho = sd(t_5_trueRho)/(sqrt(iterations)),
               se_empVar_t5_trueRho = sqrt(var_t5_trueRho)/sqrt(2*(iterations-1)),
               MC_SE_model_se_t5_trueRho = sqrt(var_v5_trueRho/(4*iterations*sqrt(mean_v5_trueRho))),
               cov_prob_t5_trueRho = sum(in_CI_t5_trueRho)/iterations,
               # t5 unbiased with true Rho properties
               mean_unbiased_t5_trueRho = mean(unbiased_t5_trueRho),
               var_unbiased_t5_trueRho = var(unbiased_t5_trueRho),
               bias_unbiased_t5_trueRho = mean_unbiased_t5_trueRho - theta,
               MCse_bias_unbiased_t5_trueRho = sqrt(var_unbiased_t5_trueRho/iterations),
               rel_bias_unbiased_t5_trueRho = bias_unbiased_t5_trueRho/theta,
               mse_unbiased_t5_trueRho = bias_unbiased_t5_trueRho^2+var_unbiased_t5_trueRho,
               mse_unbiased_t5_true_trueRho = (unbiased_t5_trueRho-mean_unbiased_t5_trueRho)^2,
               mse_unbiased_t5true_avg_trueRho = (1/iterations)*sum(mse_unbiased_t5_true_trueRho),
               mean_unbiased_v5_trueRho = mean(unbiased_v5_trueRho),
               var_unbiased_v5_trueRho = var(unbiased_v5_trueRho),
               se_unbiased_t5_trueRho = sd(unbiased_t5_trueRho)/(sqrt(iterations)),
               se_empVar_unbiased_t5_trueRho = sqrt(var_unbiased_t5_trueRho)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t5_trueRho = sqrt(var_unbiased_v5_trueRho/(4*iterations*sqrt(mean_unbiased_v5_trueRho))),
               cov_prob_unbiased_t5_trueRho = sum(in_CI_unbiased_t5_trueRho)/iterations,
               combo_cov_prob_unbiased_t5_trueRho = sum(combo_in_CI_t5_trueRho)/iterations,
               
               
               #########
               # t1 properties
               mean_t1 = mean(t_1),
               mean_unbiased_t1 = mean(unbiased_t1),
               
               var_t1 = var(t_1),
               bias_t1 = mean(t_1) - theta,
               MCse_bias_t1 = sqrt(var_t1/iterations),
               
               rel_bias_t1 = bias_t1/theta,
               mse_t1 = bias_t1^2+var_t1,
               mse_t1_true = (t_1-mean_t1)^2,
               mse_t1true_avg = (1/iterations)*sum(mse_t1_true),
               mean_v1 = mean(v_1),
               var_v1 = var(v_1),
               se_t1 = sd(t_1)/(sqrt(iterations)),
               se_empVar_t1 = sqrt(var_t1)/sqrt(2*(iterations-1)),
               ratio_var_t1 = mean_v1/var_t1,
               MC_SE_model_se_t1 = sqrt(var_v1/(4*iterations*sqrt(mean_v1))),
               cov_prob_t1 = sum(in_CI_t1)/iterations,
               # t1 unbiased properties
               mean_unbiased_t1 = mean(unbiased_t1),
               var_unbiased_t1 = var(unbiased_t1),
               bias_unbiased_t1 = mean_unbiased_t1 - theta,
               MCse_bias_unbiased_t1 = sqrt(var_unbiased_t1/iterations),
               rel_bias_unbiased_t1 = bias_unbiased_t1/theta,
               mse_unbiased_t1 = bias_unbiased_t1^2+var_unbiased_t1,
               mse_unbiased_t1_true = (unbiased_t1-mean_unbiased_t1)^2,
               mse_unbiased_t1true_avg = (1/iterations)*sum(mse_unbiased_t1_true),
               mean_unbiased_v1 = mean(unbiased_v1),
               var_unbiased_v1 = var(unbiased_v1),
               se_unbiased_t1 = sd(unbiased_t1)/(sqrt(iterations)),
               se_empVar_unbiased_t1 = sqrt(var_unbiased_t1)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t1 = sqrt(var_unbiased_v1/(4*iterations*sqrt(mean_unbiased_v1))),
               cov_prob_unbiased_t1 = sum(in_CI_unbiased_t1)/iterations,
               combo_cov_prob_unbiased_t1 = sum(combo_in_CI_t1)/iterations,
               
               
               #########
               # t2 properties
               mean_t2 = mean(t_2),
               mean_unbiased_t2 = mean(unbiased_t2),
               
               var_t2 = var(t_2),
               bias_t2 = mean(t_2) - theta,
               MCse_bias_t2 = sqrt(var_t2/iterations),
               
               rel_bias_t2 = bias_t2/theta,
               mse_t2 = bias_t2^2+var_t2,
               mse_t2_true = (t_2-mean_t2)^2,
               mse_t2true_avg = (1/iterations)*sum(mse_t2_true),
               mean_v2 = mean(v_2),
               var_v2 = var(v_2),
               se_t2 = sd(t_2)/(sqrt(iterations)),
               se_empVar_t2 = sqrt(var_t2)/sqrt(2*(iterations-1)),
               ratio_var_t2 = mean_v2/var_t2,
               MC_SE_model_se_t2 = sqrt(var_v2/(4*iterations*sqrt(mean_v2))),
               cov_prob_t2 = sum(in_CI_t2)/iterations,
               # t2 unbiased properties
               mean_unbiased_t2 = mean(unbiased_t2),
               var_unbiased_t2 = var(unbiased_t2),
               bias_unbiased_t2 = mean_unbiased_t2 - theta,
               MCse_bias_unbiased_t2 = sqrt(var_unbiased_t2/iterations),
               rel_bias_unbiased_t2 = bias_unbiased_t2/theta,
               mse_unbiased_t2 = bias_unbiased_t2^2+var_unbiased_t2,
               mse_unbiased_t2_true = (unbiased_t2-mean_unbiased_t2)^2,
               mse_unbiased_t2true_avg = (1/iterations)*sum(mse_unbiased_t2_true),
               mean_unbiased_v2 = mean(unbiased_v2),
               var_unbiased_v2 = var(unbiased_v2),
               se_unbiased_t2 = sd(unbiased_t2)/(sqrt(iterations)),
               se_empVar_unbiased_t2 = sqrt(var_unbiased_t2)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t2 = sqrt(var_unbiased_v2/(4*iterations*sqrt(mean_unbiased_v2))),
               cov_prob_unbiased_t2 = sum(in_CI_unbiased_t2)/iterations,
               combo_cov_prob_unbiased_t2 = sum(combo_in_CI_t2)/iterations,
               
               #########
               # t10 properties
               mean_t10 = mean(t_10),
               var_t10 = var(t_10),
               bias_t10 = mean(t_10) - theta,
               MCse_bias_t10 = sqrt(var_t10/iterations),
               
               rel_bias_t10 = bias_t10/theta,
               mse_t10 = bias_t10^2+var_t10,
               mse_t10_true = (t_10-mean_t10)^2,
               mse_t10true_avg = (1/iterations)*sum(mse_t10_true),
               mean_v10 = mean(v_10),
               var_v10 = var(v_10),
               se_t10 = sd(t_10)/(sqrt(iterations)),
               se_empVar_t10 = sqrt(var_t10)/sqrt(2*(iterations-1)),
               ratio_var_t10 = mean_v10/var_t10,
               MC_SE_model_se_t10 = sqrt(var_v10/(4*iterations*sqrt(mean_v10))),
               cov_prob_t10 = sum(in_CI_t10)/iterations,
               # t10 unbiased properties
               mean_unbiased_t10 = mean(unbiased_t10),
               var_unbiased_t10 = var(unbiased_t10),
               bias_unbiased_t10 = mean_unbiased_t10 - theta,
               MCse_bias_unbiased_t10 = sqrt(var_unbiased_t10/iterations),
               rel_bias_unbiased_t10 = bias_unbiased_t10/theta,
               mse_unbiased_t10 = bias_unbiased_t10^2+var_unbiased_t10,
               mse_unbiased_t10_true = (unbiased_t10-mean_unbiased_t10)^2,
               mse_unbiased_t10true_avg = (1/iterations)*sum(mse_unbiased_t10_true),
               mean_unbiased_v10 = mean(unbiased_v10),
               var_unbiased_v10 = var(unbiased_v10),
               se_unbiased_t10 = sd(unbiased_t10)/(sqrt(iterations)),
               se_empVar_unbiased_t10 = sqrt(var_unbiased_t10)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t10 = sqrt(var_unbiased_v10/(4*iterations*sqrt(mean_unbiased_v10))),
               cov_prob_unbiased_t10 = sum(in_CI_unbiased_t10)/iterations,
               combo_cov_prob_unbiased_t10 = sum(combo_in_CI_t10)/iterations,
               
               #########
               # t11 properties
               mean_t11 = mean(t_11),
               var_t11 = var(t_11),
               bias_t11 = mean(t_11) - theta,
               MCse_bias_t11 = sqrt(var_t11/iterations),
               
               rel_bias_t11 = bias_t11/theta,
               mse_t11 = bias_t11^2+var_t11,
               mse_t11_true = (t_11-mean_t11)^2,
               mse_t11true_avg = (1/iterations)*sum(mse_t11_true),
               mean_v11 = mean(v_11),
               var_v11 = var(v_11),
               se_t11 = sd(t_11)/(sqrt(iterations)),
               se_empVar_t11 = sqrt(var_t11)/sqrt(2*(iterations-1)),
               ratio_var_t11 = mean_v11/var_t11,
               MC_SE_model_se_t11 = sqrt(var_v11/(4*iterations*sqrt(mean_v11))),
               cov_prob_t11 = sum(in_CI_t11)/iterations,
               # t5 unbiased properties
               mean_unbiased_t11 = mean(unbiased_t11),
               var_unbiased_t11 = var(unbiased_t11),
               bias_unbiased_t11 = mean_unbiased_t11 - theta,
               MCse_bias_unbiased_t11 = sqrt(var_unbiased_t11/iterations),
               rel_bias_unbiased_t11 = bias_unbiased_t11/theta,
               mse_unbiased_t11 = bias_unbiased_t11^2+var_unbiased_t11,
               mse_unbiased_t11_true = (unbiased_t11-mean_unbiased_t11)^2,
               mse_unbiased_t11true_avg = (1/iterations)*sum(mse_unbiased_t11_true),
               mean_unbiased_v11 = mean(unbiased_v11),
               var_unbiased_v11 = var(unbiased_v11),
               se_unbiased_t11 = sd(unbiased_t11)/(sqrt(iterations)),
               se_empVar_unbiased_t11 = sqrt(var_unbiased_t11)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t11 = sqrt(var_unbiased_v11/(4*iterations*sqrt(mean_unbiased_v11))),
               cov_prob_unbiased_t11 = sum(in_CI_unbiased_t11)/iterations,
               combo_cov_prob_unbiased_t11 = sum(combo_in_CI_t11)/iterations
        ) %>% 
        select(-c(diff_means_post, sd_post, mean_diff_gains, var_gain_trt, var_gain_cnt, sd_pooled, cov_trt, cov_cnt,
                  corr_pooled_trt, corr_pooled_cnt, corr_pooled, part1_in_estT5, part2_in_estT5, mse_t5_true, 
                  mse_t1_true, mse_t2_true, t_1, v_1, t_2, v_2, t_5, v_5, t_5_trueRho, v_5_trueRho, mse_t5_true_trueRho,
                  LB_t1, UB_t1, in_CI_t1, LB_t2, UB_t2, in_CI_t2, LB_t5, UB_t5, in_CI_t5, LB_t5_trueRho,
                  UB_t5_trueRho, in_CI_t5_trueRho, t_11, v_11, t_10, v_10, mse_t11_true, mse_t10_true,
                  LB_t11, UB_t11, LB_t10, UB_t10, sPOOLED, in_CI_t10, in_CI_t11, 
                  
                  combo_LB_t11, combo_UB_t11, combo_in_CI_t11,
                  combo_LB_t10, combo_UB_t10, combo_in_CI_t10,
                  combo_LB_t1, combo_UB_t1, combo_in_CI_t1,
                  combo_LB_t2, combo_UB_t2, combo_in_CI_t2,
                  combo_LB_t5, combo_UB_t5, combo_in_CI_t5,
                  combo_LB_t5_trueRho, combo_UB_t5_trueRho, combo_in_CI_t5_trueRho,
                  
                  LB_unbiased_t11, UB_unbiased_t11, in_CI_unbiased_t11,
                  LB_unbiased_t10, UB_unbiased_t10, in_CI_unbiased_t10,
                  LB_unbiased_t1, UB_unbiased_t1, in_CI_unbiased_t1,
                  LB_unbiased_t2, UB_unbiased_t2, in_CI_unbiased_t2,
                  LB_unbiased_t5, UB_unbiased_t5, in_CI_unbiased_t5,
                  LB_unbiased_t5_trueRho, UB_unbiased_t5_trueRho, in_CI_unbiased_t5_trueRho,
                  
                  mse_unbiased_t5_true, mse_unbiased_t1_true, mse_unbiased_t10_true, mse_unbiased_t11_true,
                  mse_unbiased_t2_true, mse_unbiased_t5_true_trueRho,
                  
                  unbiased_t11, unbiased_v11, unbiased_t5_trueRho, unbiased_v5_trueRho,
                  unbiased_t5, unbiased_v5, unbiased_t10, unbiased_v10, unbiased_t2, unbiased_v2,
                  unbiased_t1, unbiased_v1
                  )) %>% 
        distinct()
    
    df_clean_rho_cumulative_gains <- rbind(df_clean_rho_cumulative_gains, df_rho_gains)
    
    #############################
    # ANCOVA
    df_rho_ancova <- df_clean_rho_ancova %>% 
        mutate(n_trt = n_trt,
               n_cnt = n_cnt,
               #########
               # t3 properties
               mean_t3 = mean(t_3),
               var_t3 = var(t_3),
               bias_t3 = mean(t_3) - theta,
               MCse_bias_t3 = sqrt(var_t3/iterations),
               
               rel_bias_t3 = bias_t3/theta,
               mse_t3 = bias_t3^2+var_t3,
               mse_t3_true = (t_3-mean_t3)^2,
               mse_t3true_avg = (1/iterations)*sum(mse_t3_true),
               mean_v3 = mean(v_3),
               var_v3 = var(v_3),
               se_t3 = sd(t_3)/(sqrt(iterations)),
               se_empVar_t3 = sqrt(var_t3)/sqrt(2*(iterations-1)),
               ratio_var_t3 = mean_v3/var_t3,
               MC_SE_model_se_t3 = sqrt(var_v3/(4*iterations*sqrt(mean_v3))),
               cov_prob_t3 = sum(in_CI_t3)/iterations,
               # t3 unbiased properties
               mean_unbiased_t3 = mean(unbiased_t3),
               var_unbiased_t3 = var(unbiased_t3),
               bias_unbiased_t3 = mean_unbiased_t3 - theta,
               MCse_bias_unbiased_t3 = sqrt(var_unbiased_t3/iterations),
               rel_bias_unbiased_t3 = bias_unbiased_t3/theta,
               mse_unbiased_t3 = bias_unbiased_t3^2+var_unbiased_t3,
               mse_unbiased_t3_true = (unbiased_t3-mean_unbiased_t3)^2,
               mse_unbiased_t3true_avg = (1/iterations)*sum(mse_unbiased_t3_true),
               mean_unbiased_v3 = mean(unbiased_v3),
               var_unbiased_v3 = var(unbiased_v3),
               se_unbiased_t3 = sd(unbiased_t3)/(sqrt(iterations)),
               se_empVar_unbiased_t3 = sqrt(var_unbiased_t3)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t3 = sqrt(var_unbiased_v3/(4*iterations*sqrt(mean_unbiased_v3))),
               cov_prob_unbiased_t3 = sum(in_CI_unbiased_t3)/iterations,
               combo_cov_prob_unbiased_t3 = sum(combo_in_CI_t3)/iterations,
               
               
               ###############
               # t9 properties with ESTIMATED RHO
               mean_t9 = mean(t_9),
               var_t9 = var(t_9),
               bias_t9 = mean_t9 - theta,
               MCse_bias_t9 = sqrt(var_t9/iterations),
               
               rel_bias_t9 = bias_t9/theta,
               mse_t9 = bias_t9^2+var_t9,
               mse_t9_true = (t_9-mean_t9)^2,
               mse_t9true_avg = (1/iterations)*sum(mse_t9_true),
               mean_v9 = mean(v_9),
               var_v9 = var(v_9),
               se_t9 = sd(t_9)/(sqrt(iterations)),
               se_empVar_t9 = sqrt(var_t9)/sqrt(2*(iterations-1)),
               ratio_var_t9 = mean_v9/var_t9,
               MC_SE_model_se_t9 = sqrt(var_v9/(4*iterations*sqrt(mean_v9))),
               cov_prob_t9 = sum(in_CI_t9)/iterations,
               # t9 unbiased with estimated Rho properties
               mean_unbiased_t9 = mean(unbiased_t9),
               var_unbiased_t9 = var(unbiased_t9),
               bias_unbiased_t9 = mean_unbiased_t9 - theta,
               MCse_bias_unbiased_t9 = sqrt(var_unbiased_t9/iterations),
               rel_bias_unbiased_t9 = bias_unbiased_t9/theta,
               mse_unbiased_t9 = bias_unbiased_t9^2+var_unbiased_t9,
               mse_unbiased_t9_true = (unbiased_t9-mean_unbiased_t9)^2,
               mse_unbiased_t9true_avg = (1/iterations)*sum(mse_unbiased_t9_true),
               mean_unbiased_v9 = mean(unbiased_v9),
               var_unbiased_v9 = var(unbiased_v9),
               se_unbiased_t9 = sd(unbiased_t9)/(sqrt(iterations)),
               se_empVar_unbiased_t9 = sqrt(var_unbiased_t9)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t9 = sqrt(var_unbiased_v9/(4*iterations*sqrt(mean_unbiased_v9))),
               cov_prob_unbiased_t9 = sum(in_CI_unbiased_t9)/iterations,
               combo_cov_prob_unbiased_t9 = sum(combo_in_CI_t9)/iterations,
               
              
               
               ###############
               # t9 properties with TRUE rho
               mean_t9_trueRho = mean(t_9_trueRho),
               var_t9_trueRho = var(t_9_trueRho),
               bias_t9_trueRho = mean(t_9_trueRho) - theta,
               MCse_bias_t9_trueRho = sqrt(var_t9_trueRho/iterations),
               
               rel_bias_t9_trueRho = bias_t9_trueRho/theta,
               mse_t9_trueRho = bias_t9_trueRho^2+var_t9_trueRho,
               mse_t9_true_trueRho = (t_9_trueRho-mean_t9_trueRho)^2,
               mse_t9true_avg_trueRho = (1/iterations)*sum(mse_t9_true_trueRho),
               mean_v9_trueRho = mean(v_9_trueRho),
               var_v9_trueRho = var(v_9_trueRho),
               se_t9_trueRho = sd(t_9_trueRho)/(sqrt(iterations)),
               se_empVar_t9_trueRho = sqrt(var_t9_trueRho)/sqrt(2*(iterations-1)),
               ratio_var_t9_trueRho = mean_v9_trueRho/var_t9_trueRho,
               MC_SE_model_se_t9_trueRho = sqrt(var_v9_trueRho/(4*iterations*sqrt(mean_v9_trueRho))),
               cov_prob_t9_trueRho = sum(in_CI_t9_trueRho)/iterations,
               # t9 unbiased with true Rho properties
               mean_unbiased_t9_trueRho = mean(unbiased_t9_trueRho),
               var_unbiased_t9_trueRho = var(unbiased_t9_trueRho),
               bias_unbiased_t9_trueRho = mean_unbiased_t9_trueRho - theta,
               MCse_bias_unbiased_t9_trueRho = sqrt(var_unbiased_t9_trueRho/iterations),
               rel_bias_unbiased_t9_trueRho = bias_unbiased_t9_trueRho/theta,
               mse_unbiased_t9_trueRho = bias_unbiased_t9_trueRho^2+var_unbiased_t9_trueRho,
               mse_unbiased_t9_trueRho_true = (unbiased_t9_trueRho-mean_unbiased_t9_trueRho)^2,
               mse_unbiased_t9_trueRho_true_avg = (1/iterations)*sum(mse_unbiased_t9_trueRho_true),
               mean_unbiased_v9_trueRho = mean(unbiased_v9_trueRho),
               var_unbiased_v9_trueRho = var(unbiased_v9_trueRho),
               se_unbiased_t9_trueRho = sd(unbiased_t9_trueRho)/(sqrt(iterations)),
               se_empVar_unbiased_t9_trueRho = sqrt(var_unbiased_t9_trueRho)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t9_trueRho = sqrt(var_unbiased_v9_trueRho/(4*iterations*sqrt(mean_unbiased_v9_trueRho))),
               cov_prob_unbiased_t9_trueRho = sum(in_CI_unbiased_t9_trueRho)/iterations,
               combo_cov_prob_unbiased_t9_trueRho = sum(combo_in_CI_t9_trueRho)/iterations,
               
               #########
               # t12 properties
               mean_t12 = mean(t_12),
               var_t12 = var(t_12),
               bias_t12 = mean(t_12) - theta,
               MCse_bias_t12 = sqrt(var_t12/iterations),
               
               rel_bias_t12 = bias_t12/theta,
               mse_t12 = bias_t12^2+var_t12,
               mse_t12_true = (t_12-mean_t12)^2,
               mse_t12true_avg = (1/iterations)*sum(mse_t12_true),
               mean_v12 = mean(v_12),
               var_v12 = var(v_12),
               se_t12 = sd(t_12)/(sqrt(iterations)),
               se_empVar_t12 = sqrt(var_t12)/sqrt(2*(iterations-1)),
               ratio_var_t12 = mean_v12/var_t12,
               MC_SE_model_se_t12 = sqrt(var_v12/(4*iterations*sqrt(mean_v12))),
               cov_prob_t12 = sum(in_CI_t12)/iterations,
               # t12 unbiased properties
               mean_unbiased_t12 = mean(unbiased_t12),
               var_unbiased_t12 = var(unbiased_t12),
               bias_unbiased_t12= mean_unbiased_t12 - theta,
               MCse_bias_unbiased_t12 = sqrt(var_unbiased_t12/iterations),
               rel_bias_unbiased_t12 = bias_unbiased_t12/theta,
               mse_unbiased_t12 = bias_unbiased_t12^2+var_unbiased_t12,
               mse_unbiased_t12_true = (unbiased_t12-mean_unbiased_t12)^2,
               mse_unbiased_t12true_avg = (1/iterations)*sum(mse_unbiased_t12_true),
               mean_unbiased_v12 = mean(unbiased_v12),
               var_unbiased_v12 = var(unbiased_v12),
               se_unbiased_t12 = sd(unbiased_t12)/(sqrt(iterations)),
               se_empVar_unbiased_t12 = sqrt(var_unbiased_t12)/sqrt(2*(iterations-1)),
               MC_SE_model_se_unbiased_t12 = sqrt(var_unbiased_v12/(4*iterations*sqrt(mean_unbiased_v12))),
               cov_prob_unbiased_t12 = sum(in_CI_unbiased_t12)/iterations,
               combo_cov_prob_unbiased_t12 = sum(combo_in_CI_t12)/iterations
            ) %>% 
        select(-c(var_post_group, sd_adjusted, sd_post_pooled, t_3, v_3, part1_t9, corr_pooled, part2_t9, t_9, v_9,
                  part1_t9_trueRho, part2_t9_trueRho, t_9_trueRho, v_9_trueRho, LB_t3, UB_t3, 
                  in_CI_t3, LB_t9, UB_t9, in_CI_t9, LB_t9_trueRho, UB_t9_trueRho, in_CI_t9_trueRho,
                  mse_t3_true, mse_t9_true, mse_t9_trueRho, mse_t9_true_trueRho, t_12, v_12, 
                  LB_t12, UB_t12, mse_t12_true, var_pre_group, sPOOLED_anc, in_CI_t12,
                  
                  
                  combo_LB_t12, combo_UB_t12, combo_in_CI_t12,
                  combo_LB_t9, combo_UB_t9, combo_in_CI_t9,
                  combo_LB_t9_trueRho, combo_UB_t9_trueRho, combo_in_CI_t9_trueRho,
                  combo_LB_t3, combo_UB_t3, combo_in_CI_t3,
                  
                  LB_unbiased_t12, UB_unbiased_t12, in_CI_unbiased_t12,
                  LB_unbiased_t9, UB_unbiased_t9, in_CI_unbiased_t9,
                  LB_unbiased_t9_trueRho, UB_unbiased_t9_trueRho, in_CI_unbiased_t9_trueRho,
                  LB_unbiased_t3, UB_unbiased_t3, in_CI_unbiased_t3,
                  
                  mse_unbiased_t9_true, mse_unbiased_t12_true, mse_unbiased_t9_trueRho_true,
                  mse_unbiased_t3_true,
                  
                  unbiased_t12, unbiased_v12, unbiased_t9_trueRho, unbiased_v9_trueRho,
                  unbiased_t9, unbiased_v9, unbiased_t3, unbiased_v3
        )) %>% 
        distinct()
    
    df_clean_rho_cumulative_ancova <- rbind(df_clean_rho_cumulative_ancova, df_rho_ancova)

}
df_clean_rho_cumulative_gains
print("--------------------------------------------------------------------------------")
df_clean_rho_cumulative_ancova




saveRDS(df_clean_rho_cumulative_gains, "df_rho_4_delta_8_gains.rds")   
saveRDS(df_clean_rho_cumulative_ancova, "df_rho_4_delta_8_ancova.rds")







