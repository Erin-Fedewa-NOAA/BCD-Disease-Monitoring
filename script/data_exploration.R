#Data exploration: 
  #Maturity/crab size/sample size 
  #Index site maps by year
  #Quick glance at covariates of interest for modeling drivers

# Author: Erin Fedewa

# load ----
library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(akgfmaps)
library(ggridges)
library(patchwork)

#load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#####################################################
# data wrangling

dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(pcr_result %in% c(1, 0)) %>%
  rename(pcr = pcr_result,
                station = gis_station,
                latitude = mid_latitude,
                longitude = mid_longitude,
                depth = bottom_depth,
                temperature = gear_temperature,
                index = index_site) %>%
  mutate(year = as.factor(year),
                sex = as.factor(sex),
                index = as.factor(index),
                station = as.factor(station)) -> opilio.dat 

###############################################
#Spatial maps

## SET COORDINATE REFERENCE SYSTEMS (CRS) --------------------------------------
in.crs <- "+proj=longlat +datum=NAD83" #CRS is in lat/lon
map.crs <- "EPSG:3338" # final crs for mapping/plotting: Alaska Albers

## LOAD SHELLFISH ASSESSMENT PROGRAM GEODATABASE -------------------------------
survey_gdb <- "./data/SAP_layers" 
survey_strata <- terra::vect(survey_gdb, layer = "EBS.NBS_surveyarea")
#EBS/NBS Boundary line
boundary <- st_read(layer = "EBS_NBS_divide", survey_gdb)

## LOAD ALASKA REGION LAYERS (FROM AKGFMAPS R package) -----------------------------------
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ebs_survey_areas <- ebs_layers$survey.area
ebs_survey_areas$survey_name <- c("Eastern Bering Sea", "Northern Bering Sea")

#Transform crab data into spatial data frame
opilio.dat  %>% 
  group_by(year, general_location, latitude, longitude) %>%
  summarise(n_crab=n()) %>% 
  # Convert lat/long to an sf object
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326)) %>%
  #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
  st_transform(crs = st_crs(3338)) -> map_dat

#Add crab data to map
ggplot() +
  geom_sf(data = ebs_layers$survey.grid, fill=NA, color=alpha("grey80"))+
  geom_sf(data = ebs_survey_areas, fill = NA) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
  #add crab data
  geom_sf(data=map_dat, aes(size = n_crab, color=general_location), alpha = .6) +
  geom_sf(data= boundary, linewidth = 1, color = "grey40") +
  scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                     breaks = ebs_layers$lon.breaks) +
  scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                     breaks = ebs_layers$lat.breaks) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  facet_wrap(~year) +
  scale_color_manual(values = c("#034e7b", "#238b45")) +
  theme(plot.margin = margin(0,-5,0,-5)) +
  theme(axis.text=element_text(size=8)) +
  theme(axis.text.x=element_blank()) 
ggsave("./figures/map_effort.png")




















###############################################

#examine <70mm opilio cpue distribution                                    
ggplot(opilio.dat, aes(snow70under_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#4th root
ggplot(opilio.dat, aes(snow70under_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#log
ggplot(opilio.dat, aes(log(snow70under_cpue))) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#examine protocol target size opilio cpue distribution to decide on transf.
ggplot(opilio.dat, aes(snowimm_cpue)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#4th root
ggplot(opilio.dat, aes(snowimm_cpue^0.25)) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#log
ggplot(opilio.dat, aes(log(snowimm_cpue))) +
  geom_histogram(bins = 30, fill = "grey", color = "black")

#Transform both CPUE indices with 4th root
opilio.dat %>%
  mutate(fourth.root.cpue70 = snow70under_cpue^0.25,
         fouth.root.cpueimm = snowimm_cpue^0.25) -> opilio.dat 

nrow(opilio.dat) # 1520 samples!

#Sample sizes
opilio.dat %>%
  group_by(year, index, station) %>%
  count() %>%
  print(n=100) #Only one crab sampled at  some stations..though better than tanner

# BCS+/- by index/year
opilio.dat %>%
  select(year, index, pcr) %>%
  group_by(year, index) %>%
  summarise(PCR_0 = sum(pcr == 0),
            PCR_1 = sum(pcr == 1)) %>%
  pivot_longer(cols = c(-year, -index)) -> plot2

#% +/- stacked barplot
ggplot(plot2, aes(fill=name, y=value, x=year)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ index)

#Sample sizes by maturity
dat %>%
  #filter(maturity != "NA") %>%
  group_by(year,maturity) %>%
  count() %>%
  ggplot() +
  geom_bar(aes(x=as.factor(maturity), y= n), stat='identity') +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Maturity Status", y = "Sample size") #lots of missing maturity info in 14/15

#Size range sampled across years
opilio.dat %>% 
  summarize(avg_cw = mean(size, na.rm=T), 
            max_cw = max(size, na.rm=T), 
            min_cw = min(size, na.rm=T))

#Size composition sampled by index site/yr
opilio.dat %>%
  mutate(Sex = recode_factor(sex, '1' = "M", '2' = "F")) %>%
  group_by(year, index) %>%
  ggplot() +
  geom_histogram(aes(x=size, fill=Sex), position = "stack", bins=50) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
  facet_wrap(~year) +
  theme_bw() +
  labs(x= "Snow crab carapace width (mm)", y = "Count")
ggsave("./figs/opilio_size.png", width=6.75)

 
#Size-frequency distribution of uninfected vrs infected
opilio.dat %>%
  ggplot(aes(size, fill=as.factor(pcr), color=as.factor(pcr))) +  
  geom_histogram(position="identity",alpha=0.5) +
  theme_bw()

#Percent prevalence by size bin
opilio.dat %>% 
  mutate(size_bin = cut(size, breaks=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120))) %>%
  group_by(size_bin) %>%
  summarise(Prevalance = (sum(pcr)/n())*100) %>%
  filter(size_bin != "NA") -> size
  ggplot(size, aes(as.factor(size_bin), Prevalance)) +
  geom_col() 
  
#This is something to keep in mind when interpreting changes in prev. across
#site and year- at site 6 prevalence was very low, but is likely due to 
#most samples being taken from mature males (despite protocol specifying imm).
#Without pulling them, it's difficult to interpret changes in prevalence as 
#a true measure of the disease prev, vrs. just a factor of the subset of population 
#being sampled. 

#############################################
#Covariate data exploration 

#Plot range of observed data by year/index site 
opilio.dat %>%
  group_by(year, index, station) %>%
  summarise(temperature = mean(temperature),
            depth = mean(depth),
            julian = mean(julian),
            cpue = mean(snow70under_cpue)) -> plot 

#Temperature
ggplot(plot, aes(temperature)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#CPUE
ggplot(plot, aes(cpue)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#Depth
ggplot(plot, aes(depth)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index)

#Julian Day 
ggplot(plot, aes(julian)) +
  geom_histogram(bins = 12, fill = "dark grey", color = "black") +
  facet_grid(year ~ index) # big differences among index areas

#Julian date vrs temperature 
ggplot(plot, aes(julian, temperature)) +
  geom_point() # Stronger correlation than tanner temp vrs date 

#Plot explanatory variables as predictors of proportion BCS+ by year/station 
opilio.dat %>%
  group_by(year, station) %>%
  summarise(size = mean(size), 
            julian = mean(julian),
            temperature = mean(temperature),
            CPUE70 = mean(fourth.root.cpue70),
            CPUEimm = mean(fouth.root.cpueimm),
            proportion.positive = sum(pcr == 1) / n()) -> plot3

# Julian day vrs  %positive plot
ggplot(plot3, aes(julian, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #Seasonality effect

#Mean size-at-station vrs %positive plot 
ggplot(plot3, aes(size, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #size effect

#Temp-at-station vrs %positive plot 
ggplot(plot3, aes(temperature, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam") #temperature effect

# <70mm CPUE-at-station vrs %positive plot 
ggplot(plot3, aes(CPUE70, proportion.positive)) +
  geom_point() + 
  #facet_wrap(~year) +
  geom_smooth(method = "gam")

# Immature CPUE-at-station vrs %positive plot 
ggplot(plot3, aes(CPUEimm, proportion.positive)) +
  geom_point() + 
  facet_wrap(~year) +
  geom_smooth(method = "gam") #Very similar to above plot 


                              





