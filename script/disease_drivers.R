#Investigate factors influencing BCD infection probability (conventional PCR 0/1 diagnosis) 
  #and infection intensity (dPCR # of copies detected) in snow crab using Bayesian multivariate models  

#Note: using hierarchical modeling approach from Fedewa et al 2025 - but including 
  #NBS data in models. Not testing for sex effects because only females targeted in
  #2017-2019 NBS collections (no idea why.....)

#Haven't yet run infection load models- limited sample size and majority of low 
  #grade infections need some thought as to best approach for modeling drivers-
  #more rationale for splitting into two bins, high vrs low intensity infection, 
  #and use binomial models 

#load packages
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(marginaleffects)
library(emmeans)
library(MARSS)
library(corrplot)
library(factoextra)
library(patchwork)
library(modelr)
library(broom.mixed)
library(pROC)
library(ggthemes)
library(tidybayes)
library(RColorBrewer)
library(knitr)
library(loo)
library(sjPlot)
library(prettyglm)
source("./script/stan_utils.R")


# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#load color palettes
my_colors <- RColorBrewer::brewer.pal(8, "GnBu")

########################################
#Functions for model diagnostic checks

# Fisher Pearson skew 
skew <- function(y){
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}


#Plot PPC test statistic plots
posterior_stat_plot <- function(obs, model, samples=1000, statistic="mean"){
  fig <- ppc_stat(y = obs, 
                  yrep = posterior_predict(model, ndraws = samples),  #matrix of draws from posterior distribution
                  stat = statistic)
  
  return(fig)
}

########################################
#Data Manipulation
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(pcr_result %in% c(1, 0),
         size > 0, 
         sex %in% c(1,2),
         index_site != 2, #3 snow crab samples collected at a tanner crab index site
         no_positive_partitions >= 0) %>% 
  select(pcr_result, size, general_location, no_positive_partitions, sex, index_site, year, gis_station, julian, 
         mid_latitude, bottom_depth,gear_temperature, cpue) %>%
  rename(pcr = pcr_result, 
         station = gis_station,
         latitude = mid_latitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station)) -> opilio.dat 

#Check for missing data for PCA
opilio.dat %>%
  select(size, sex, year, index, station, julian, latitude, depth, temperature) %>%
  filter(!complete.cases(.)) #Looks good  

#Assess collinearity b/w covariates 
opilio.dat %>%
  group_by(year, station) %>%
  summarise(julian = mean(julian),
                   depth = mean(depth),
                   latitude = mean(latitude),
                   temperature = mean(temperature),
                   cpue = mean(cpue)) -> corr.opilio

cor(corr.opilio[,3:7]) #latitude and julian day highly correlated 
corrplot(cor(corr.opilio[,3:7]), method = 'number') 

# plot covariates
corr.opilio %>%
  rename("Depth (m)" = depth,
         "Latitude (N)" = latitude,
         "Bottom Temperature (C)" = temperature,
         "Snow Crab Density (CPUE)" = cpue) %>%
  pivot_longer(4:7, names_to = "variable", values_to = "data") %>% 
  ggplot(aes(julian, data)) +
  geom_point(aes(color = as.factor(year))) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = T, alpha = 0.2, 
              color = "black", lwd = 0.3) +
 theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(x= "Day of Year") +
  ggtitle("Snow Crab Covariates") +
  theme(plot.title = element_text(hjust = 0.5))

####################################################
#Probability of infection models: response is 0/1 diagnosis from conventional PCR
  #includes 2014-2023 data, both EBS and NBS

#Model 1: base model with size, julian day and random year/index intercept 
  #note that we're pooling NBS and EBS data here to look at drivers
opilio1_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + (1 | year/index)) 

## Show default priors
get_prior(opilio1_formula, opilio.dat, family = bernoulli(link = "logit"))

## fit binomial model 1
opilio1 <- brm(opilio1_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, 
               warmup = 1500, iter = 6000,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio1, file = "./output/opilio1.rds")
opilio1 <- readRDS("./output/opilio1.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio1$fit)
neff_lowest(opilio1$fit)
rhat_highest(opilio1$fit)
summary(opilio1)
bayes_R2(opilio1) #.23

#Diagnostic Plots
plot(opilio1, ask = FALSE)
plot(conditional_smooths(opilio1), ask = FALSE)
mcmc_plot(opilio1, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio1)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio1, pars = c("b_Intercept", "bs_ssize_1", "bs_spc1_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio1)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
y_obs <- opilio.dat$pcr #Observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio1) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio1, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1 

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio1)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio1) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well (though should preferably be >0.80)

#Final diagnostic: Correct classification rate
#calculate the predicted probabilities of 0/1 in the original data from the fitted model
  Pred <- predict(opilio1, type = "response")
  Pred <- if_else(Pred[,1] > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(opilio.dat, pcr)) #`pull` results in a vector
#correct classification rate
  sum(diag(ConfusionMatrix))/sum(ConfusionMatrix) #Model correctly classifies 77% of observations
  ConfusionMatrix

###################################################
#Model 2: add temperature to base model

opilio2_formula <- bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(temperature, k = 4) + (1 | year/index))  

opilio2 <- brm(opilio2_formula,
               data = opilio.dat,
               family = bernoulli(link = "logit"),
               cores = 4, chains = 4, 
               warmup = 1500, iter = 6000,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio2, file = "./output/opilio2.rds")
opilio2 <- readRDS("./output/opilio2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio2$fit)
neff_lowest(opilio2$fit)
rhat_highest(opilio2$fit)
summary(opilio2)
bayes_R2(opilio2) #0.23

#Diagnostic Plots
plot(opilio2, ask = FALSE)
plot(conditional_smooths(opilio2), ask = FALSE)
mcmc_plot(opilio2, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio2)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio2, pars = c("b_Intercept", "bs_ssize_1", "bs_sjulian_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio2)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
y_obs <- opilio.dat$pcr #Observed values

color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio2) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio2, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1 

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio2)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio2) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well
mean(auc)

#Final diagnostic: Correct classification rate (similar to above)
#calculate the predicted probabilities of 0/1 in the original data from the fitted model
Pred <- predict(opilio2, type = "response")
Pred <- if_else(Pred[,1] > 0.5, 1, 0)
ConfusionMatrix <- table(Pred, pull(opilio.dat, pcr)) #`pull` results in a vector
#correct classification rate
sum(diag(ConfusionMatrix))/sum(ConfusionMatrix) #Model correctly classifies 78% of observations
ConfusionMatrix

# model comparison
loo(opilio1, opilio2) #Temperature doesn't improve predictive capacity, but we'll 
  #keep in for now 

###############################################################################################
#Model 3: Add CPUE to base model 2

opilio3_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(temperature, k = 4) + 
                         s(cpue, k=4) + (1 | year/index))  

opilio3 <- brm(opilio3_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, 
               warmup = 1500, iter = 6000,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(opilio3, file = "./output/opilio3.rds")
opilio3 <- readRDS("./output/opilio3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio3$fit)
neff_lowest(opilio3$fit)
rhat_highest(opilio3$fit)
summary(opilio3)
bayes_R2(opilio3) #0.24

#Diagnostic Plots
plot(opilio3, ask = FALSE)
plot(conditional_smooths(opilio3), ask = FALSE)
mcmc_plot(opilio3, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio3)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio3, pars = c("b_Intercept", "bs_ssize_1", "bs_scpue_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio3)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio3) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio3, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1  

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio3)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 

# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio3) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well
mean(auc)

# model comparison
loo(opilio1, opilio2, opilio3) #CPUE improves predictive capacity of base model, 
  #though still fairly small improvement

######################################################
# Model 4: Add depth to model 3 

opilio4_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(temperature, k = 4) + 
                         s(cpue, k=4) + s(depth, k = 4) + (1 | year/index))                      

opilio4 <- brm(opilio4_formula,
               data = opilio.dat,
               family =bernoulli(link = "logit"),
               cores = 4, chains = 4, 
               warmup = 1500, iter = 6000,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save Output
saveRDS(opilio4, file = "./output/opilio4.rds")
opilio4 <- readRDS("./output/opilio4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(opilio4$fit)
neff_lowest(opilio4$fit)
rhat_highest(opilio4$fit)
summary(opilio4)
bayes_R2(opilio4) #0.24

#Diagnostic Plots
plot(opilio4, ask = FALSE)
plot(conditional_smooths(opilio4), ask = FALSE)
mcmc_plot(opilio4, type = "areas", prob = 0.95)
mcmc_rhat(rhat(opilio4)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(opilio4, pars = c("b_Intercept", "bs_ssize_1", "b_sex2"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(opilio4)) #Effective sample size: All ratios > 0.1

#Posterior Predictive check to assess predictive performance: mean and skewness 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, opilio4) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,opilio4, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1  

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(opilio4)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy fairly good 


# Compute AUC for predicting prevalence with the model
preds <- posterior_epred(opilio4) #posterior draws
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Model discriminates fairly well
mean(auc)

#Model comparison
loo(opilio1, opilio2, opilio3, opilio4)  #Opilio4 marginally better, although again- 
  #models very similar 

######################################################
#Extract and plot conditional effects of each predictor from best model, opilio4
#conditioning on the mean for all other predictors, yr/site effects ignored 

#Julian Day
## 95% CI
ce1s_1 <- conditional_effects(opilio4, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opilio4, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opilio4, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = julian, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Day of Year", y = "Probability of infection") +
  theme_bw() -> dayplot

#Depth
## 95% CI
ce1s_1 <- conditional_effects(opilio4, effect = "depth", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opilio4, effect = "depth", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opilio4, effect = "depth", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$depth
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$depth[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$depth[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$depth[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$depth[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = depth, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Depth (m)", y = "") +
  theme_bw() -> depthplot

#CPUE
## 95% CI
ce1s_1 <- conditional_effects(opilio4, effect = "cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opilio4, effect = "cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opilio4, effect = "cpue", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$cpue[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = cpue, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Snow crab CPUE", y = "Probability of infection") +
  theme_bw() -> cpueplot

#temperature
## 95% CI
ce1s_1 <- conditional_effects(opilio4, effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(opilio4, effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(opilio4, effect = "temperature", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temperature
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temperature[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temperature[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temperature[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temperature[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = opilio.dat, aes(x = temperature, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Bottom Temperature", y = "Probability of infection") +
  theme_bw() -> tempplot

#Combine plots 
(dayplot + depthplot) / (cpueplot + tempplot) + (sizeplot + plot_spacer()) +
  plot_annotation(title = "Factors associated with BCD Occurrence", 
                  theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("./figures/drivers_occurrence.png", height=9)


################################################################################
#Drivers of infection intensity models: response is # of positive partitions from dPCR
  #Note that this is 2018+ data and infected individuals only (much smaller dataset), 
  #and we're filtering out runs that didn't amplify DNA (<2%)

#Data Manipulation
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(no_positive_partitions > 0, #infected crab only
         nssu_pcr == 1, #excluding runs that didn't amplify dna
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  select(pcr_result, size, no_positive_partitions, general_location, sex, index_site, year, gis_station, julian, 
         mid_latitude, bottom_depth,gear_temperature, cpue) %>%
  rename(pcr = pcr_result, 
         station = gis_station,
         latitude = mid_latitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station)) -> dpcr.dat 

#look at histos/density curves of # of positive partitions- assign breaks via R based on modal 
  #analysis? 
dpcr.dat %>%
  ggplot()+
  geom_histogram(aes(no_positive_partitions))
#Wow, looks like we've got a lot of very low grade infections- not even sure 
#sample sizes will be large enough for modeling 

#To do: bin samples by Hamish categories and run models, response= prob of infection intensity
  #bin number of positive partitions by 1-10 partitions (light infection),
  #11-100 (moderate infection), 101-1000 (moderate-heavy infection) and 1000+ (heavy infection)
  #is there a less arbiratry way to bin? Don't think we have enough data/peaks for modal analysis
  #to identify breaks in the data, might be better to just use binomial approach below

#Given very right and left skewed data, bin samples in two bins, light 
  #infection (1-4000?) vrs heavy (4000+?), response=prob of high intensity infection

#Not sure how to structure these models- first is logistic, second is 
  #binomial (e.g 0% to 100%, two categories only) ?


#Let's first test a preliminary logistic regression model to look at how infection
  #load influences the probability of an infection
  #but the problem here is that we're using standard and dPCR products, which we know
  #have different sensitivities-so probably need to create a 0/1 infection status
  #column using positive partitions >0/=1 and not use standard PCR data (see outlier issue below)

model <- glm(pcr ~ no_positive_partitions, data = opilio.dat, family = binomial(link = "logit"))
summary(model) 
plot(model) #outliers are most commonly due to large dPCR load but standard PCR=0 diagnosis
plot(opilio.dat$pcr ~ opilio.dat$no_positive_partitions)
#right duh, we can't fit a binomial model to this 

pretty_relativities(feature_to_plot= 'no_positive_partitions',
                    model_object = model,
                    relativity_label = 'Probability of infection',
                    upper_percentile_to_cut = 0.1)

#Finding observations that are outliers
hist(model$residuals)
absres <- abs(model$residuals) #Absolute value of the residual.  You can do any transformation of the residuals you find appropriate.
q <- quantile(absres, probs = .95) #The 95th quantile of this distribution
opilio.dat[which(absres>q),]