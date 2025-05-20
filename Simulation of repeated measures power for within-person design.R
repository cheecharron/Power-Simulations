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
  n <- 30  # Number of participants
  n_time <- 4  # Number of time points
  icc <- 0.33  # Desired ICC
  effect_size <- 0.43  # Full effect size (Cohen's d at time 3)
  missing = 0.0 # rate of missingness
  
  # Fixed effects for time (i.e., within-person effect only)
  beta0 <- rnorm(1, mean=0, sd=1)  # Intercept (baseline value)
  beta1 <- rnorm(1, mean = effect_size, sd = 1)
  
  # Simulate random intercept variance for subjects
  subject_var <- icc * 1  # Between-subject variance
  residual_var <- (1 - icc) * 1  # Within-subject residual variance
  
  # Random effects for subjects
  subject_intercepts <- rnorm(n, mean = 0, sd = sqrt(subject_var))
  
  # Create a data frame
  df <- data.frame(
    subject = rep(1:n, each = n_time),
    time = rep(1:n_time, times = n)
  )
  

  # Generate outcome variable with interaction
  df$y <- beta0 + 
    beta1 * df$time +  # Main effect of time
    subject_intercepts[df$subject] +  # Random intercepts for subjects
    rnorm(n * n_time, mean = 0, sd = sqrt(residual_var))  # Random noise
  
  #remove missing data
  df$miss2 = rbinom(n*n_time, 1, missing*n_time/2)
  df$miss3 = rbinom(n*n_time, 1, missing*n_time/2)
  
  df = df %>%
    mutate(y = ifelse(miss2==1 & time==2,NA,y)) %>%
    mutate(y = ifelse(miss3==1 & time==3,NA,y))
  
  # Fit a mixed-effects model
  model <- lmer(y ~  time + (1 | subject), data = df)
  
  summary(model)
 
  #extract p-value for end-of-treatment
  p = summary(model)$coefficients[2,5] %>% as.numeric
  
  ICC = icc_fun(lmer(y ~ (1 | subject), data=df))[[1]]
  
  results[row,'p'] = p
  results[row,'ICC'] = ICC

}

#this returns the percentage of p-values < .05
ecdf(results$p)(0.05)


