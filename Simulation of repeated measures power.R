# The goal of this file is to simulate trajectories to illustrate how
# much power there would be to find differences between treatment groups
# with repeated observations and varying levels of missing data.

library(dplyr)
library(lme4)
library(sqldf)
library(emmeans)
library(lmerTest)
library(simr)

icc_fun = function(model){
  #### compute ICC
  var.components = as.data.frame(VarCorr(model))$vcov
  ICC = var.components[1]/sum(var.components)
  
  #### find out average cluster size
  id.name = names(coef(model))
  clusters = nrow(matrix(unlist((coef(model)[id.name]))))
  n = length(residuals(model))
  average.cluster.size = n/clusters
  
  #### compute design effects
  design.effect = 1+ICC*(average.cluster.size-1)
  
  #### return stuff
  list(icc=ICC, design.effect=design.effect)
  
}


results = data.frame(matrix(ncol=2,nrow=0))
colnames(results) = c('p','ICC')

for (row in seq(1:1000)){

  print(row)
  
  # Parameters
  n <- 164  # Number of participants
  n_time <- 3  # Number of time points
  icc <- 0.5  # Desired ICC
  effect_size <- 0.5  # Full effect size (Cohen's d at time 3)
  missing = 0.20 # rate of missingness
  
  # Fixed effects for time, treatment, and their interaction
  beta0 <- 0  # Intercept (baseline value)
  beta1 <- 0  # No group difference at time 1
  beta2_time2 <- effect_size / 2  # Group difference at time 2 is half the effect size
  beta2_time3 <- effect_size  # Group difference at time 3 (full effect size)
  
  # Simulate random intercept variance for subjects
  subject_var <- icc * 1  # Between-subject variance
  residual_var <- (1 - icc) * 1  # Within-subject residual variance
  
  # Random effects for subjects
  subject_intercepts <- rnorm(n, mean = 0, sd = sqrt(subject_var))
  
  # Create a data frame
  df <- data.frame(
    subject = rep(1:n, each = n_time),
    time = rep(1:n_time, times = n),
    treatment = rep(sample(c(0, 1), n, replace = TRUE), each = n_time)
  )
  
  # Define group differences at each time point
  df$group_effect <- ifelse(df$time == 1, beta1,  # No group difference at time 1
                            ifelse(df$time == 2, beta2_time2, beta2_time3))  # Group difference at times 2 and 3
  
  # Generate outcome variable with interaction
  df$y <- beta0 + 
    df$group_effect * df$treatment +  # Group effect by time
    0.2 * df$time +  # Main effect of time
    subject_intercepts[df$subject] +  # Random intercepts for subjects
    rnorm(n * n_time, mean = 0, sd = sqrt(residual_var))  # Random noise
  
  #remove missing data
  df$miss2 = rbinom(n*n_time, 1, missing*n_time/2)
  df$miss3 = rbinom(n*n_time, 1, missing*n_time/2)
  
  df = df %>%
    mutate(y = ifelse(miss2==1 & time==2,NA,y)) %>%
    mutate(y = ifelse(miss3==1 & time==3,NA,y))
  
  df$time = factor(df$time)
  df$treatment = factor(df$treatment)
  
  # Fit a mixed-effects model
  model <- lmer(y ~ treatment * time + (1 | subject), data = df)
  
  summary(model)
 
  #extract p-value for end-of-treatment
  p = pairs(emmeans(model, "treatment", by="time")) %>% data.frame() %>% 
    filter(time==3) %>%
    dplyr::select(p.value) %>% as.numeric
  
  ICC = icc_fun(lmer(y ~ (1 | subject), data=df))[[1]]
  
  results[row,'p'] = p
  results[row,'ICC'] = ICC

}

#this returns the percentage of p-values < .05
ecdf(results$p)(0.05)


