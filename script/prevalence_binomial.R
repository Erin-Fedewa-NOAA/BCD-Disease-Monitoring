## Annual BCD prevalence estimated via a binomial model
  #This approach allows us to "correct" for Julian day and crab size effect, because
  #we know prevalence increases with Julian day and in smaller crab, and our sampling
  #design doesn't account for this
  
#Author: Erin Fedewa

library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

# load PCR data
dat = read.csv("./data/pcr_haul_master.csv")

#########################################################################
#Eastern Bering Sea model
dat %>%
  filter(general_location == "EBS",
         size > 0) %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date
  filter(pcr_result %in% c(1, 0), 
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_latitude, 
         bottom_depth, gear_temperature, cpue) %>%
  rename(pcr = pcr_result,
         station = gis_station,
         latitude = mid_latitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station)) -> ebs.opilio

## Fit a binomial prevalence model for EBS only using 0/1 PCR status 
ebs_opilio_formula <-  bf(pcr ~ year  + s(size, k = 4) + s(julian, k=4)) 

ebs_snow_pcr = brm(ebs_opilio_formula,
                   data = ebs.opilio,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, 
                   warmup = 1500, iter = 6000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(ebs_snow_pcr, file = "./output/ebs_snow_pcr.rds")
ebs_snow_pcr = readRDS("./output/ebs_snow_pcr.rds")

pp_check(ebs_snow_pcr, type = "dens_overlay", ndraws = 100)
bayes_R2(ebs_snow_pcr) #0.14
plot(conditional_smooths(ebs_snow_pcr), ask = FALSE)

#extract year conditional effect 
ce_snow_pcr = brms::conditional_effects(ebs_snow_pcr, effect = "year")
df_snow = ce_snow_pcr$year

#plot estimated annual prevalence in EBS
ggplot(df_snow) +
    geom_point(aes(year, estimate__ * 100), color = "#2171b5", size = 3, position = dodge) +
    geom_errorbar(aes(year, ymin = lower__ * 100, ymax = upper__ * 100),
                  width = 0.3, size = 0.5, position = dodge, color = "#2171b5") +
    ylab("Bitter Crab Disease Prevalence (%)") + xlab("") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank()) +
    theme(legend.title = element_blank()) +
    theme(axis.text.x = element_text(size=11)) +
    ggtitle("Eastern Bering Sea Snow Crab")+ 
    theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,100) -> ebs_plot
#This story is quite a bit different than observed annual prevalence estimates!
  #points to the importance of accounting for crab size during sampling efforts

#########################################################
#Northern Bering Sea Model
  #Note: I tried to fit a model with a year*location interaction but model didn't converge, 
  #likely due to unbalance sampling design and the NBS sampled more infrequently?

dat %>%
  filter(general_location == "NBS",
         size > 0) %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date
  filter(pcr_result %in% c(1, 0)) %>% 
  select(pcr_result, size, sex, index_site, year, gis_station, julian, mid_latitude, 
         bottom_depth, gear_temperature, cpue) %>%
  rename(pcr = pcr_result,
         station = gis_station,
         latitude = mid_latitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station)) -> nbs.opilio

## Fit a binomial prevalence model for NBS only using 0/1 PCR status 
nbs_opilio_formula <-  bf(pcr ~ year  + s(size, k = 4) + s(julian, k=4)) 

nbs_snow_pcr = brm(nbs_opilio_formula,
                   data = nbs.opilio,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, 
                   warmup = 1500, iter = 6000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(nbs_snow_pcr, file = "./output/nbs_snow_pcr.rds")
nbs_snow_pcr = readRDS("./output/nbs_snow_pcr.rds")

pp_check(nbs_snow_pcr, type = "dens_overlay", ndraws = 100)
bayes_R2(nbs_snow_pcr) #0.17
plot(conditional_smooths(nbs_snow_pcr), ask = FALSE)
#size effect much less apparent here b/c mostly small, immature crab in NBS

#extract year conditional effect 
ce_snow_nbs_pcr = brms::conditional_effects(nbs_snow_pcr, effect = "year")
df_snow_nbs = ce_snow_nbs_pcr$year

#plot estimated annual prevalence in NBS
ggplot(df_snow_nbs) +
  geom_point(aes(year, estimate__ * 100), color = "#238b45", size = 3, position = dodge) +
  geom_errorbar(aes(year, ymin = lower__ * 100, ymax = upper__ * 100),
                width = 0.3, size = 0.5, position = dodge, color = "#238b45") +
  ylab("Bitter Crab Disease Prevalence (%)") + xlab("") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=11)) +
  ggtitle("Northern Bering Sea Snow Crab")+ 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,100) -> nbs_plot
#Hmmm not sure why our point estimate for 2017 has so much associated uncertainty- 
  #probably something to troubleshoot in future iterations

#quick visual look until then
nbs.opilio %>% 
  filter(year == 2017) %>%
  group_by(station) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) %>%
  ggplot(aes(station, Prevalance)) +
  geom_col() 
#okay, definitely something to look into. All stations sampled in 2017 have 
  #prevalence > 50% but model is predicting 25% prevalence with lots of uncertainty

####################################################
#combine plots and save
ebs_plot + nbs_plot 
ggsave("./figures/annual_prev.png")
