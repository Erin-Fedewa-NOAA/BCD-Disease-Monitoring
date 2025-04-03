#Compare testing accuracy of PCR, dPCR and visual diagnosis
#Sensitivity and specificity analysis between PCR and visual diagnosis with 
#additional years

#This script includes lots of questions for Hamish, and filtering criteria 
  #for all further conventional and dPCR analyses 

# Author: Erin Fedewa

# load ----
library(tidyverse)
library(lubridate)
library(viridis)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#Load color palette
my_colors <- c("#8fd7d7","#01665e", "#c99b38", "#eddca5")
my_colors2 <- c("#3594cc", "#8cc5e3","#c99b38", "#eddca5")
my_colors3 <- RColorBrewer::brewer.pal(8, "Blues")[c(3,5,6,8)]

############################################
#Sample sizes
dat %>%
  filter(pcr_result %in% c(1, 0), 
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  summarize(total = n(),
            pcr_pos = sum(pcr_result==1),
            vis_pos = sum(visual_positive==1, na.rm=TRUE))

#Calculate sensitivity and specificity
dat %>%
  filter(index_site != 2,
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) %>% #true negative
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#15% probability of diagnosing visually + when PCR is +
  #this is up from 7% reported in Fedewa et al for 2015-2017 dataset
#99% probability of diagnosing visually neg when hemato is absent 
  #25 crab visually positive but no hemato- either incorrect visual diagnosis, or
  #symptoms presenting for a different bug? 


#Plot Classification Accuracy  
dat %>%
  filter(index_site != 2,
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) %>% #true negative
  summarise(TP_n = sum(sen_spec=="TP", na.rm=TRUE)/n()*100,
            FP_n = sum(sen_spec=="FP", na.rm=TRUE)/n()*100,
            FN_n = sum(sen_spec=="FN", na.rm=TRUE)/n()*100,
            TN_n = sum(sen_spec=="TN", na.rm=TRUE)/n()*100) %>%
  rename("True Positive" = "TP_n", "False Positive" = "FP_n",
         "True Negative" = "TN_n", "False Negative" = "FN_n") %>%
  pivot_longer(c(1:4), names_to="Classification", values_to = "Percent") %>%
  mutate(Diagnosis = case_when(Classification=="False Negative" | Classification=="True Positive" ~ "Infected", 
                               Classification=="False Positive" | Classification=="True Negative" ~ "Not Infected")) %>%
  ggplot(aes(Diagnosis, Percent, fill=Classification)) +
    geom_bar(position = "stack", stat="identity") +
    theme_bw() +
    scale_fill_manual(values = my_colors2) 
ggsave("./figures/sens_spec.png")

########################################################
#data exploration of dPCR data 

#number of samples with runs that failed to amplify DNA
dat %>%
  filter(year >= 2018,
         nssu_pcr == 0) %>%
  summarize(total = n()) #n=30, <2%
#we should be removing these samples for both conventional and dPCR analyses  
  #moving forward, very small sample size 

#number of samples with faint positives (often questionable and subjective)
dat %>%
  filter(year >= 2018,
         faint_pos == 1) %>%
  summarize(total = n()) #n=134, ~8.8%

#sample sizes of conventional PCR data
  #NOTE: dropping nssu fails from here on out
dat %>%
  filter(pcr_result %in% c(1, 0), 
         nssu_pcr == 1,
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  group_by(year, general_location) %>%
  summarize(total = n(),
            pos = sum(pcr_result==1),
            neg = sum(pcr_result==0),
            prevalence = (pos/total)*100) %>%
  mutate(method = "Standard PCR") -> pcr

#sample sizes of conventional PCR data filtered for faint positives
  #note that we're keeping in faint positive samples that were diagnosed positive with dPCR-
  #need to check with Hamish on this 
dat %>%
  filter(pcr_result %in% c(1, 0), 
         nssu_pcr == 1, 
         faint_pos == 0 | faint_pos == 1 & no_positive_partitions > 0, 
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  group_by(year, general_location) %>%
  summarize(total = n(),
            pos = sum(pcr_result==1),
            neg = sum(pcr_result==0),
            prevalence = (pos/total)*100) %>%
  mutate(method = "Standard PCR minus faint positives") -> pcr_filter

#sample sizes of all digital PCR data
dat %>%
  filter(year >= 2018, 
         nssu_pcr == 1,
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  group_by(year, general_location) %>%
  summarize(total = n(),
            pos = length(no_positive_partitions[no_positive_partitions > 0]),
            neg = length(no_positive_partitions[no_positive_partitions == 0]),
            prevalence = (pos/total)*100) %>%
  mutate(method = "Digital PCR") -> dig_pcr

#combine datasets and plot prevalence between methods 
pcr %>%
  rbind(dig_pcr) %>%
  rbind(pcr_filter) %>%
  filter(year > 2017) %>%
  ggplot(aes(fill=method, y=prevalence, x=year)) + 
  geom_bar(position="dodge", stat="identity") + 
  facet_wrap(~general_location) + 
  theme_bw() + 
  scale_fill_manual(values = my_colors) +
  labs(y="Prevalence (%)", x="")
ggsave("./figures/prev_by_method.png")
#so seems that standard pcr prevalance estimates are always higher- better
  #sensitivity for low grade infections due to dilution needed for dPCR?

#plot sample sizes of # of positive
pcr %>%
  rbind(dig_pcr) %>%
  rbind(pcr_filter) %>%
  filter(year > 2017) %>%
  ggplot(aes(fill=method, y=pos, x=year)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values = my_colors) +
  labs(y="Total # of samples diagnosed as positive", x="")
ggsave("./figures/no_positive_by_method.png")

# of samples where standard PCR=1 and dPCR=0
dat %>%
  filter(year >= 2018, 
         nssu_pcr == 1,
         index_site != 2,
         pcr_result == 1, 
         no_positive_partitions == 0) %>% 
  summarize(total = n()) %>%
  mutate(method = "Standard PCR=1") -> pcr1

#of samples where standard PCR=1 and dPCR=0, but removing faint positives 
  #from standard PCR data
dat %>%
  filter(year >= 2018, 
         nssu_pcr == 1,
         index_site != 2,
         pcr_result == 1, 
         faint_pos == 0,
         no_positive_partitions == 0) %>% 
  summarize(total = n()) %>%
  mutate(method = "Standard PCR=1 minus faint positives") -> pcr2

#and plot both
pcr1 %>%
  rbind(pcr2) %>%
  ggplot(aes(x=method, y=total)) +
  geom_bar(position="dodge", stat="identity") + 
  theme_bw() + 
  labs(y="Total # of samples", x="", 
       title="Total # of samples diagnosed as infected via standard \n PCR but uninfected via dPCR") +
  theme(plot.title = element_text(size = 9)) +
  theme(axis.text.x = element_text(size=7.5)) 
ggsave("./figures/standardPCR_vrs_dPCR.png")
#so it seems like removing faint positives greatly reduces the discrepancy- support
  #that these faint positives are probably not true infections?
 
# of samples where standard PCR=0 and dPCR=1
dat %>%
  filter(year >= 2018, 
         nssu_pcr == 1,
         index_site != 2,
         pcr_result == 0, 
         no_positive_partitions > 0) %>% 
  summarize(total = n())
#but we've got almost the same number of samples in the opposite direction!


#Follow up questions for Hamish:
#1 This is a first pass, but how do we determine limit of detection of dPCR, 
  #and compare sensitivities of standard PCR vrs dPCR?
#2 When developing dPCR assay, shouldn't we be using different dilutions to quantify 
  #the limits of detection and limits of quantification- because this is important 
  #to our understanding of disease onset, progression and occurrence 
  #i.e. wouldn't it make the most sense to use 1:10 dilutions on standard PCR products
  #(can you even dilute these?) as well as dPCR?
  #because we can't really compare diagnosis/sensitives from the two methods without 
  #3) The concern of reduced sensitivity of dPCR (i.e. 95 samples infected via standard PCR
  #but uninfected via dPCR is a concern, though removing subjective faint positives helps)
  #on the other end of the spectrum though, there are 75 samples diagnosed positive via dPCR, 
  #but not standard PCR. Hmmmm....
#4) all faint positives are diagnosed as 1 for standard PCR, but should we be dropping all 
  #"pcr_result=1 & faint_positive=1" samples? - 134 samples as faint positive, but 28 of these 
  #have positive partitions > 1, so should we leave these in as confirmation for infection?
#5) We should be removing nssu fails from entire dataset, correct? 
            