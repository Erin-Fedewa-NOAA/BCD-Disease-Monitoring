#Correlation between visual BCD EBS-wide prevalence timeseries and index site prevalence 
   
#To do:
#calculate annual visual prevalence via ESP script - add NBS
#use model estimates for EBS and NBS annual prevelance 
#test for syncrony in timeseries with correlation and overlay of two timeseries on same plot 
#test for relationship with linear model 

#subset visual dataset for same size/sampling criteria for blood draws
#filter visual dataset for index site stations only for more refined comparision

# Author: Erin Fedewa
# last updated: 2023/1/3

# load ----
library(tidyverse)
library(lubridate)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")


#EBS haul data 
crab_ebs <- read.csv("./Data/crabhaul_opilio.csv")
#EBS strata data 
strata_ebs <- read_csv("./Data/crabstrata_opilio.csv")

#NBS haul data 
crab_nbs <- read.csv("./Data/crabhaul_opilio_nbs.csv")
#NBS strata data 
strata_nbs <- read_csv("./Data/crabstrata_opilio_nbs.csv")

########################################
#Calculate BCD visual prevalence from EBS survey timeseries 
  #To do this, we'll filter survey data for sampling criteria for index site 
  #sampling (i.e. immature females, <= 94mm males)

#compute cpue 
crab_ebs %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(HAUL_TYPE == 3,
         ((SEX == 2 & CLUTCH_SIZE == 0) |
            (SEX == 1 & WIDTH <= 94))) %>% 
  mutate(bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, bcs) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue = ncrab / AREA_SWEPT) -> cpue 
#add zero catch stations to bcs=T/F datasets
cpue %>% 
  filter(bcs == FALSE) %>%
  right_join(crab_ebs %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  replace_na(list(bcs = FALSE)) %>%
  rbind(cpue %>% 
          filter(bcs == TRUE) %>%
          right_join(crab_ebs %>% 
                       mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
                       filter(HAUL_TYPE ==3) %>%
                       distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
          replace_na(list(cpue = 0)) %>%
          replace_na(list(ncrab = 0)) %>%
          replace_na(list(bcs = TRUE))) -> catch
#Now our BCS categories should each contain zero catch stations so we can group by 
#bcs category to calculate abundances 

#BCS prevelance
catch %>%
  right_join(strata_ebs %>%
               rename(GIS_STATION=STATION_ID, YEAR=SURVEY_YEAR)) %>%
  group_by(YEAR, STRATUM, bcs) %>%
  summarise(total_area = mean(TOTAL_AREA),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(YEAR, bcs) %>%
  #sum across strata
  summarise(Total_abun=sum(abundance_mil)) %>% 
  filter(!is.na(bcs)) %>%
  #calculate prevalence 
  group_by(YEAR) %>%
  summarise(Perc_Prevalance = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+
                                    (Total_abun[bcs==TRUE])))*100) %>%
  filter(YEAR > 1988) -> prev_ebs

#plot
prev_ebs %>%
  ggplot(aes(x = YEAR, y = Perc_Prevalance)) +
  geom_point(aes(), size=3) +
  geom_line(aes(), size=1) +
  labs(y = "Disease Prevalence (%)", x = "") +
  theme_bw() +
  geom_hline(aes(yintercept = mean(Perc_Prevalance, na.rm=TRUE)), linetype = 5)+
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) 

########################################
#Calculate BCD visual prevalence from NBS survey timeseries 
#To do this, we'll filter survey data for sampling criteria for index site 
#sampling (i.e. immature females, <= 94mm males)

#compute cpue 
crab_ebs %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(HAUL_TYPE == 3,
         ((SEX == 2 & CLUTCH_SIZE == 0) |
            (SEX == 1 & WIDTH <= 94))) %>% 
  mutate(bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, bcs) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue = ncrab / AREA_SWEPT) -> cpue 
#add zero catch stations to bcs=T/F datasets
cpue %>% 
  filter(bcs == FALSE) %>%
  right_join(crab_ebs %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  replace_na(list(bcs = FALSE)) %>%
  rbind(cpue %>% 
          filter(bcs == TRUE) %>%
          right_join(crab_ebs %>% 
                       mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
                       filter(HAUL_TYPE ==3) %>%
                       distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
          replace_na(list(cpue = 0)) %>%
          replace_na(list(ncrab = 0)) %>%
          replace_na(list(bcs = TRUE))) -> catch
#Now our BCS categories should each contain zero catch stations so we can group by 
#bcs category to calculate abundances 

#BCS prevelance
catch %>%
  right_join(strata_ebs %>%
               rename(GIS_STATION=STATION_ID, YEAR=SURVEY_YEAR)) %>%
  group_by(YEAR, STRATUM, bcs) %>%
  summarise(total_area = mean(TOTAL_AREA),
            mean_cpue = mean(cpue),
            abundance_mil = mean(total_area) * mean_cpue / 1000000) %>%
  group_by(YEAR, bcs) %>%
  #sum across strata
  summarise(Total_abun=sum(abundance_mil)) %>% 
  filter(!is.na(bcs)) %>%
  #calculate prevalence 
  group_by(YEAR) %>%
  summarise(Perc_Prevalance = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+
                                                        (Total_abun[bcs==TRUE])))*100) %>%
  filter(YEAR > 1988) -> prev_ebs

#plot
prev_ebs %>%
  ggplot(aes(x = YEAR, y = Perc_Prevalance)) +
  geom_point(aes(), size=3) +
  geom_line(aes(), size=1) +
  labs(y = "Disease Prevalence (%)", x = "") +
  theme_bw() +
  geom_hline(aes(yintercept = mean(Perc_Prevalance, na.rm=TRUE)), linetype = 5)+
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) 


