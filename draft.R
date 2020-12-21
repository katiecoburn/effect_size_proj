# devtools::install_github("arcaldwell49/Superpower", force = TRUE)
library(Superpower)
library(tidyverse)
# want es to be 0.4
# 
#   melt(as.data.frame(mvrnorm(n = n, mu = mu, Sigma = as.matrix(sigmatrix), 
#empirical = TRUE)))$value
#})
design_result <- ANOVA_design(design = "2b*2w",
                              n = 10, 
                              mu = c(0.8, 0.8, 0.2, 06), 
                              sd = 0.5, 
                              r = c(0.8, 0.5, 0.3, 0.9, 0.2, 0.1),
                              labelnames = c("intervention", "treatment", "control", "block", "one", "two"),
                              plot = TRUE)

test <- ANOVA_exact(design_result, alpha_level = 0.05)$dataframe %>% 
  as_tibble() %>% 
  select(intervention, block, y)
# test %>%
#   group_by(intervention) %>% 
#   summarise(var(y))
aov(test$y ~ test$intervention * test$block) %>% 
  summary()
