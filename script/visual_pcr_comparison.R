#Assess the correlation between visual BCD EBS-wide prevalence timeseries and 
  #index site PCR prevalence 
   
# Author: Erin Fedewa

# load ----
library(tidyverse)
library(lubridate)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#load model-based prevalence estimates (see prevalence_binomial.R)
prev_mod <- read.csv("./output/prevalence_estimates_model.csv")

#Create look up tables with index site stations
index <- read.csv("./data/index_station_lookup.csv")

index %>% 
  filter(design_year == "2014_2023",
         region == "EBS") %>% 
  pull(gis_station) -> index_ebs

index %>% 
  filter(design_year == "2014_2023",
         region == "NBS") %>% 
  pull(gis_station) -> index_nbs

#EBS haul data 
crab_ebs <- read.csv("./Data/crabhaul_opilio.csv")
#EBS strata data 
strata_ebs <- read_csv("./Data/crabstrata_opilio.csv")

#NBS haul data 
crab_nbs <- read.csv("./Data/crabhaul_opilio_nbs.csv")
#NBS strata data 
strata_nbs <- read_csv("./Data/crabstrata_opilio_nbs.csv")

########################################
#Calculate EBS-wide BCD visual prevalence  
  #To do this, we'll filter survey data for sampling criteria for index site 
  #sampling (i.e. immature females & <= 94mm males)

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

#BCS prevalence
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
  summarise(Perc_Prevalence = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+
                                    (Total_abun[bcs==TRUE])))*100) %>%
  filter(YEAR > 1988) %>%
  complete(YEAR = full_seq(YEAR, period = 1), fill = list(Value = NA)) %>%
  mutate(region = "EBS") -> prev_ebs

#plot
prev_ebs %>%
  ggplot(aes(x = YEAR, y = Perc_Prevalence)) +
  geom_point(aes(), size=3) +
  geom_line(aes(), size=1) +
  labs(y = "Disease Prevalence (%)", x = "") +
  theme_bw() +
  geom_hline(aes(yintercept = mean(Perc_Prevalence, na.rm=TRUE)), linetype = 5)+
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12)) 

########################################
#Calculate BCD visual prevalence from NBS survey timeseries 
#It appears that 2017-2019 protocols specified to collect immature females >30mm
  #only. We'll filter by size but not worry about females only 

#compute cpue 
crab_nbs %>% 
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(((SEX == 2 & CLUTCH_SIZE == 0 & WIDTH > 30) |
            (SEX == 1 & WIDTH > 30))) %>% 
  mutate(bcs = ifelse(DISEASE_CODE != 2 | is.na(DISEASE_CODE), F, T)) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT, bcs) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue = ncrab / AREA_SWEPT) -> cpue_nbs
#add zero catch stations to bcs=T/F datasets
cpue_nbs %>% 
  filter(bcs == FALSE) %>%
  right_join(crab_nbs %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
              distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  replace_na(list(bcs = FALSE)) %>%
  rbind(cpue_nbs %>% 
          filter(bcs == TRUE) %>%
          right_join(crab_nbs %>% 
                       mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
                       distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
          replace_na(list(cpue = 0)) %>%
          replace_na(list(ncrab = 0)) %>%
          replace_na(list(bcs = TRUE))) -> catch_nbs
#Now our BCS categories should each contain zero catch stations so we can group by 
#bcs category to calculate abundances 

#BCS prevelance
catch_nbs %>%
  right_join(strata_nbs %>%
               rename(YEAR=SURVEY_YEAR)) %>%
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
  summarise(Perc_Prevalence = (Total_abun[bcs==TRUE]/((Total_abun[bcs==FALSE])+
                                                        (Total_abun[bcs==TRUE])))*100) %>%
  filter(YEAR > 1988) %>%
  complete(YEAR = full_seq(YEAR, period = 1), fill = list(Value = NA)) %>%
  mutate(region = "NBS") -> prev_nbs

#plot
prev_nbs %>%
  ggplot(aes(x = YEAR, y = Perc_Prevalence)) +
  geom_point(aes(), size=3) +
  geom_line(aes(), size=1) +
  labs(y = "Disease Prevalence (%)", x = "") +
  theme_bw() +
  geom_hline(aes(yintercept = mean(Perc_Prevalence, na.rm=TRUE)), linetype = 5)+
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12))

#Combine timeseries and save as output
prev_ebs %>%
  full_join(prev_nbs) %>%
  rename_with(tolower) -> prev_visual
  write.csv(prev_visual, "./output/prevalence_estimates_visual.csv")

#Compare visual prevalence to PCR prevalence (model based estimates from conventional PCR)
ggplot() +
  geom_point(data=prev_visual, aes(year, perc_prevalence, group=region, color=region)) +
  geom_line(data=prev_visual, aes(year, perc_prevalence, group=region, color=region)) +
  geom_point(data=prev_mod, aes(year, prevalence_perc, group=region, color=region)) +
  geom_line(data=prev_mod, aes(year, prevalence_perc, group=region, color=region)) +
  theme_bw() + 
  labs(y="Disease Prevalence (%)", x="")
#hard to compare on same scale b/c visual prevalence is magnitudes lower 

#stacked plots
prev_mod %>%
  select(year, region, prevalence_perc) %>%
  rename(perc_prevalence = prevalence_perc) %>%
  mutate(method = "Molecular Survelliance") %>%
  rbind(prev_visual %>%
          mutate(method = "Visual Survelliance: EBS-wide")) %>%
  ggplot(aes(year, perc_prevalence, group=region, color=region)) +
  geom_point() +
  geom_line() +
  theme_bw() + 
  labs(y="Disease Prevalence (%)", x="") +
  facet_wrap(~method, scales="free_y") 

#And now let's look at correlation - combining both regions here
prev_mod %>%
  select(year, region, prevalence_perc) %>%
  rename(perc_prevalence = prevalence_perc) %>%
  mutate(method = "Molecular Survelliance") %>%
  rbind(prev_visual %>%
          mutate(method = "Visual Survelliance: EBS-wide")) %>%
  pivot_wider(names_from = method, values_from = perc_prevalence) -> comb_prev

r <- round(cor(comb_prev$`Molecular Survelliance`, comb_prev$`Visual Survelliance: EBS-wide`, 
               use="pairwise.complete.obs"), 2)
p <- cor.test(comb_prev$`Molecular Survelliance`, comb_prev$`Visual Survelliance: EBS-wide`)$p.value
 
ggplot(comb_prev, aes(y=`Molecular Survelliance`, x=`Visual Survelliance: EBS-wide`)) + 
  #geom_point() +  
   geom_text(aes(label=year)) +
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=2.5, y=4.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=2.5, y=8, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() +
  labs(x="Prevalence estimated from visual survelliance", 
       y="Prevalence estimated from molecular survelliance")
#For a first pass this is actually much better than I'd expected! 


#####################################################################
#And now we'll repeat the same process as above, but calculate visual prevalence
  #within index site stations only using our lookup table to see if this improves
  #fit- though we'd hope index sites are representative of EBS/NBS population as a whole

#stopping here, but to do: use same approach as above but filter for EBS and then 
  #NBS haul files for stations in index site look up. Calculate prevalence from 
  #these stations only and then assess relationship with model-based PCR prevalences

gis_station %in% index_ebs
