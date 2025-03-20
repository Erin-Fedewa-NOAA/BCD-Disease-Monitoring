# notes ----
#GOAL: Create "pcr_combined_master.csv" dataset for modeling/analyses
#1) Add a column for maturity using clutch codes/chela heights
     #Distribution based cutlines via Richar and Foy 2022
#2) Add in haul data (bottom temp, depth, lat/long, day sampled)
#3) Calculate snow crab CPUE at each station
   
# Author: Erin Fedewa

# load ----
library(tidyverse)

afsc_master <- read.csv("./data/PCR_2014_2017.csv")
vims_master <- read.csv("./data/PCR_2018_2023.csv")

#############################
#Combine VIMS and AFSC Patho datasets
vims_master %>%
  mutate(PCR_result = as.integer(PCR_result),
         Index_Site = as.integer(Index_Site)) %>%
  bind_rows(afsc_master) %>%
  filter(Species_Name == "Chionoecetes opilio",
         Prioritization == 1) -> pcr_master

#Determine male maturity via distribution-based cutline method/clutch codes
pcr_master %>%
  mutate(maturity = case_when((Sex == 2 & Clutch > 0) ~ 1, #mature female (EBS & NBS)
                                         (Sex == 2 & Clutch == 0) ~ 0, #immature female (EBS & NBS)
                                         #Define conditions for EBS immature/mature male using cutlines 
                                         (Sex == 1 & General_Location == "EBS" & 
                                            log(Chela) < -2.20640 + 1.13523 * log(Size)) | #immature male EBS
                                           (Sex == 1 & Size < 50 & General_Location == "EBS") ~ 0,
                                           (Sex == 1 & General_Location == "EBS" &
                                            log(Chela) >= -2.20640 + 1.13523 * log(Size)) ~ 1, #mature male EBS
                                         #NBS male cutlines
                                         (Sex == 1 & General_Location == "NBS" & 
                                            log(Chela) < -1.916947 + 1.070620 * log(Size))| 
                                           (Sex == 1 & Size < 40 & General_Location == "NBS") ~ 0, #immature male NBS
                                         (Sex == 1 & General_Location == "NBS" &
                                            log(Chela) >= -1.916947 + 1.070620 * log(Size)) ~ 1)) -> pcr_mat

#############################
#Append EBS & NBS haul data 
ebs_haul <- read.csv("./data/crabhaul_opilio.csv")
nbs_haul <- read.csv("./data/crabhaul_opilio_nbs.csv")

#Combine files 
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR %in% c(2014:2023)) %>%
  select(VESSEL, CRUISE, START_DATE, HAUL, MID_LATITUDE, MID_LONGITUDE,GIS_STATION,
         BOTTOM_DEPTH,GEAR_TEMPERATURE) %>%
  distinct() -> haul

#Join haul and PCR datasets 
pcr_mat %>% 
  left_join(haul, by = c("CRUISE", "VESSEL", "HAUL")) %>%
  rename_with(tolower) -> mat_haul

#################################
#Add in CPUE data for each station 
  
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR %in% c(2014:2023)) %>%
  rename_with(tolower) %>%
  group_by(cruise, gis_station, area_swept) %>% 
  summarise(cpue = sum(sampling_factor, na.rm = T) / mean(area_swept)) %>%
  right_join(mat_haul, by = c("cruise", "gis_station")) -> cpue

##################################################
#Write new master csv                               
write_csv(cpue, file="./data/pcr_haul_master.csv")
  
  

