#Compare testing accuracy of PCR, dPCR and visual diagnosis
#Sensitivity and specificity analysis between PCR and visual diagnosis with 
#additional years

#With 2018+ data:
#also- how to determine limit of detection of dPCR, and compare sensitivities of 
  #cPCR vrs dPCR? 
#things to look at: compare diagnosis/sample size between cPCR, case when nssu_pcr=0 (no PCR for
  #amplification, only n=30), case when nssu_pcr=0 AND faint_positive=F (n=134)
#Then bin number of positive partitions by 0 partitions (negative), 1-10 partitions (light infection),
#11-100 (moderate infection), 101-1000 (moderate-heavy infection) and 1000+ (heavy infection) - double
#check distribution and max of positive partitions. Is there a less arbitrary way to bin these?

# Author: Erin Fedewa
# last updated: 2023/1/3

# load ----
library(tidyverse)
library(lubridate)

# load PCR data 
dat <- read.csv("./data/pcr_haul_master.csv")

#Load color palette
my_colors <- c("#8fd7d7","#01665e", "#c99b38", "#eddca5")
my_colors2 <- c("#3594cc", "#8cc5e3","#c99b38", "#eddca5")

############################################
#Sample sizes
dat %>%
  filter(pcr_result %in% c(1, 0), 
         index_site != 2) %>% #3 snow crab samples collected at a tanner crab index site
  summarize(total = n(),
            pcr_pos = sum(pcr_result==1),
            vis_pos = sum(visual_positive==1))

#Calculate sensitivity and specificity
dat %>%
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
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

#7% probability of diagnosing visually + when PCR is +
  #99% probability of diagnosing visually neg when hemato is absent 

###############################################
#Does specificity/sensitivity differ by crab size/seasonality? 
#let's assess variation in sensitivity in subgroups of 
#the population (i.e. just using the median julian date/size to split data)

dat %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) -> dat2 #true negative
range(dat2$julian)
median(dat2$julian) #198
range(dat2$size, na.rm=T)
median(dat2$size, na.rm=T) #55mm

######Effect of Seasonality on sensitivity: 

#Crab sampled on day 198 or earlier 
dat2 %>%
  filter(julian <= 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
          FP_n = nrow(filter(.,sen_spec=="FP")),
          FN_n = nrow(filter(.,sen_spec=="FN")),
          TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Crab sampled post day 198 
dat2 %>%
  filter(julian > 198) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

######Effect of crab size on sensitivity: 

#Crab <=55mm carapace width
dat2 %>%
  filter(size <= 55) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Crab >55mm carapace width 
dat2 %>%
  filter(size > 55) %>%
  summarise(TP_n = nrow(filter(.,sen_spec=="TP")),
            FP_n = nrow(filter(.,sen_spec=="FP")),
            FN_n = nrow(filter(.,sen_spec=="FN")),
            TN_n = nrow(filter(.,sen_spec=="TN"))) %>%
  mutate(sens = TP_n/(TP_n+FN_n), #Sensitivity of visual diagnosis
         spec = TN_n/(FP_n+TN_n)) #Specificity of visual diagnosis

#Sample size of visually positive crab is so small that we really don't see much 
  #difference in sensitivity by date sampled or crab size 

######################################################
#Plot Classification Accuracy  

#All data
dat %>%
  filter(index_site %in% c(1:6),
         year %in% c(2015:2017),
         sex %in% c(1, 2),
         pcr_result %in% c(1, 0)) %>%
  mutate(sen_spec = case_when(pcr_result==1 & visual_positive==1 ~ "TP", #true positive
                              pcr_result==0 & visual_positive==1 ~ "FP", #false positive
                              pcr_result==1 & visual_positive==0 ~ "FN", #false negative
                              pcr_result==0 & visual_positive==0 ~ "TN")) %>% #true negative
  summarise(TP_n = sum(sen_spec=="TP")/n()*100,
            FP_n = sum(sen_spec=="FP")/n()*100,
            FN_n = sum(sen_spec=="FN")/n()*100,
            TN_n = sum(sen_spec=="TN")/n()*100) %>%
  mutate(Test = print("All Samples")) -> total

#Crab sampled on day 198 or earlier 
dat2 %>%
  filter(julian <= 198) %>%
  summarise(TP_n = sum(sen_spec=="TP")/n()*100,
            FP_n = sum(sen_spec=="FP")/n()*100,
            FN_n = sum(sen_spec=="FN")/n()*100,
            TN_n = sum(sen_spec=="TN")/n()*100) %>%
  mutate(Test = print("June-July Samples")) -> pre

#Crab sampled post day 198 
dat2 %>%
  filter(julian > 198) %>%
  summarise(TP_n = sum(sen_spec=="TP")/n()*100,
            FP_n = sum(sen_spec=="FP")/n()*100,
            FN_n = sum(sen_spec=="FN")/n()*100,
            TN_n = sum(sen_spec=="TN")/n()*100) %>%
  mutate(Test = print("July-Aug Samples")) -> post 
  
#Combine datasets and plot 
total %>%
  rbind(post) %>%
  rbind(pre) %>%
  rename("True Positive" = "TP_n", "False Positive" = "FP_n",
         "True Negative" = "TN_n", "False Negative" = "FN_n") %>%
  pivot_longer(c(1:4), names_to="Classification", values_to = "Percent") %>%
  mutate(Diagnosis = case_when(Classification=="False Negative" | Classification=="True Positive" ~ "Infected", 
                               Classification=="False Positive" | Classification=="True Negative" ~ "Not Infected")) %>%
  ggplot(aes(Diagnosis, Percent, fill=Classification)) +
    geom_bar(position = "stack", stat="identity") +
    theme_bw() +
    scale_fill_manual(values = my_colors2) +
    facet_wrap(~factor(Test, levels=c("All Samples", "June-July Samples", "July-Aug Samples")))
ggsave("./figs/sens_spec.png", width=7)

#Follow up questions for Hamish:
#1 When developing dPCR assay, shouldn't we be using different dilutions to quantify 
  #the limits of detection and limits of quantification- because this is important 
  #to our understanding of disease onset, progression and occurrence 
#wouldn't it make the most sense to use PCR on 1:10 dilutions as well as dPCR?
  #because we can't really compare diagnosis/sensitivies from the two methods 
  #without doing this. And then to point 1, should we run these same samples 



