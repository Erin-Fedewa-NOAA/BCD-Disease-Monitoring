# notes ----
#Baseline estimates of Hematodinium infection prevalence in EBS & NBS monitoring sites
  
# Author: Erin Fedewa

# load ----
library(tidyverse)
library(patchwork)
library(ggthemes)
library(lubridate)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

###########################################
#Data wrangling 
dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(pcr_result %in% c(1, 0), 
         index_site != 2) %>% #three snow crab samples collected at a tanner crab index site
  select(species_name, pcr_result, size, sex, general_location, index_site, year, gis_station, 
         julian, mid_longitude, bottom_depth, gear_temperature, cpue) %>%
  rename(pcr = pcr_result,
         station = gis_station,
         longitude = mid_longitude,
         depth = bottom_depth,
         temperature = gear_temperature,
         index = index_site) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         index = as.factor(index),
         station = as.factor(station),
         depth = as.numeric(depth)) -> dat2

#Size range of crab sampled by species
dat2 %>%
  group_by(species_name, year, general_location) %>%
  summarise(Min_size = min(size, na.rm=T),
            Max_size = max(size, na.rm=T),
            Avg_size = mean(size, na.rm=T))

#Calculate prevalence by year
dat2 %>%
  group_by(species_name, year, general_location) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) 

#Calculate prevalence and sample sizes within index sites
dat2 %>%
  group_by(species_name, index, year) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) -> prev

#Calculate sample size and binomial CI's
  a <- .05 #significance level
dat2 %>%
  group_by(species_name, index, year) %>%
  summarise(n = n()) %>%
  full_join(prev) %>%
  mutate(CI_upper = Prevalance + (qnorm(1-a/2))*sqrt((1/n)*Prevalance*(100-Prevalance)),
         CI_lower = Prevalance + (-qnorm(1-a/2))*sqrt((1/n)*Prevalance*(100-Prevalance))) -> prev_n
  
#Prevalence Plot by species/site/year
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[2:7]

prev_n %>%
  mutate(species = case_when(species_name == "Chionoecetes opilio" ~ "Snow crab")) %>%
  ggplot(aes(x=index, y=Prevalance)) +
  geom_bar(aes(fill=index, group=index),stat='identity', show.legend = FALSE) +
  scale_fill_manual(values = my_colors) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.2, color="grey") + 
  #geom_text(aes(label=paste0("n=",n), group=index), vjust = -2.2, size=2) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x="Monitoring Site", y="Prevalence (%)") +
  facet_wrap(~year) +
  theme_bw() 
ggsave("./figures/opilio_obs_prev.png")

#Prevalence plot by site/size version 2 
snow_colors <- RColorBrewer::brewer.pal(8, "Blues")[c(3,5,6,8)]

prev_n %>%
  mutate(species = case_when(species_name == "Chionoecetes opilio" ~ "Snow crab")) %>%
  ggplot(aes(x=year, y=Prevalance, fill=index)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75) + 
  scale_fill_manual(values = snow_colors) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.2, color="dark grey", position = position_dodge(width = .9)) + 
  #geom_text(aes(year, Prevalence, label=paste0("n=",n)), size=2,
            #position = position_dodge(width = .9), vjust=-5) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x="", y="Prevalence (%)", fill="Monitoring Site", title="Snow Crab") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(panel.grid.major = element_line())
ggsave("./figures/opilio_obs_prev2.png",dpi=300)






 
  
  











