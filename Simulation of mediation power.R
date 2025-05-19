#simulate power for mediation effects

library(sqldf)
library(dplyr)
library(MASS)
library(rstatix)
library(QuantPsyc)
library(ggplot2)

n = 49

iterations=50000

#this provides power on the direct, unmediated effect

df = data.frame(matrix(ncol=3,nrow=0))
colnames(df) = c('iteration','t','p')

for (i in 1:iterations){
  print(i)
  df[i,1] = i
  
  a_outcome = rnorm(n, 0, 1)
  b_outcome = rnorm(n, .6, 1)
  
  ttest = t.test(a_outcome, b_outcome, var.equal = TRUE)
  df[i,2] = ttest$statistic
  df[i,3] = ttest$p.value
}

nrow(df[df$t<0 & df$p<.05,])/iterations

#############################################################
# this provides power on the mediated effect accounting for #
# 'med' proportion of the direct, unmediated effect         #
#############################################################

# Define parameters

simulate_mediation <- function(N = 98, d_total = 0.57, med_perc = 0.5) {
  
  # 1. Treatment group (coded 0/1)
  treat <- rep(c(-.995, .995), each = N / 2)
  
  # 2. Total effect = d_total
  # Convert to correlation between treatment and outcome
  total_effect <- d_total/sqrt(2)
  
  # 3. Indirect effect = proportion of total effect
  ind_effect <- total_effect * med_perc
  
  # 4. Setting a to same association as total effect and solve for b such that a*b = indirect effect
  a <- total_effect
  b <- ind_effect / a
  
  # 5. Simulate mediator: M = a * treat + error
  M <- a * treat + rnorm(N, mean = 0, sd = 1)
  
  # 6. Direct effect = total - indirect
  direct_effect <- total_effect - ind_effect
  
  # 7. Simulate outcome: Y = b * M + direct_effect * treat + error
  Y <- b * M + direct_effect * treat + rnorm(N, mean = 0, sd = 1)
  
  # Return as data frame
  data.frame(treat = treat, M = M, Y = Y)
}


df = data.frame(matrix(ncol=8,nrow=0))
colnames(df) = c('iteration','c_cohensd','c_dir','a','b','c-prime','ind_p','ind_perc')
i=1

while (i <= 10000){
#  i=1
  
  samp = simulate_mediation(med_perc=.5)
  colnames(samp) = c('treatment','mediator','outcome')
  
  #standardize variables
  samp = samp %>%
    mutate_at(vars(c('treatment','mediator','outcome')), scale)
 
  cohensd = (mean(samp[samp$treatment<0,'outcome'])-mean(samp[samp$treatment>0,'outcome']))/sd(samp$outcome)
  
  if (cohensd < -.62 | cohensd > -.52) next 
  print(i)
  
  df[i,1] = i
  df[i,2] = cohensd
    
  mod = lm(outcome~treatment, data=samp)
  coef = summary(mod)$coefficients
  c = coef[2,1]
  
  df[i,3] = c
  
  mod1 = glm(mediator~treatment, data=samp)
  coef1 = summary(mod1)$coefficients 
  a = coef1[2,1]
  a_se = coef1[2,2]
  
  df[i,4] = a

  mod2 = glm(outcome~treatment+mediator, data=samp)
  coef2 = summary(mod2)$coefficients
  b = coef2[3,1]
  b_se = coef2[3,2]
  cp = coef2[2,1]
  cp_se = coef2[2,2]
  
  df[i,5] = b
  df[i,6] = cp
  
  df[i,7] = dt((a*b)/sqrt(a_se^2+b_se^2), df=n*2-2)
  
  df[i,8] = (a*b)/(a*b+cp)
  
  i=i+1
}

df$ind_perc2 = with(df, as.numeric(sprintf('%.2f',ind_perc)))
df$sig = ifelse(df$ind_p<.05,1,0)
df_agg = sqldf('select ind_perc2, avg(sig) as prop_sig
                from df
                group by ind_perc2')

summary(df2)

ggplot(df_agg) +
  geom_point(aes(x = ind_perc2, y = prop_sig))


######### modeling ############

model_log <- glm(sig ~ ind_perc, data = df, family = binomial)

df$predprob = predict(model_log, df, "response")

ggplot(df) +
  geom_point(aes(x = ind_perc2, y = predprob), size=.3) +
  labs(x = "Proportion of Total Effect Explained by Indirect Effect", y = "Power") +
  coord_cartesian(xlim = c(.5, 1)) +
  scale_y_continuous(labels = scales::percent, breaks = seq(.2,1,.2)) +
  theme_bw()
  
setwd('C:/Users/VHADURPaulD/OneDrive - Department of Veterans Affairs/Documents/R')
ggsave('power plot.jpeg', width=5, height=3, units='in')

