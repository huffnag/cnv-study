
# Set path to data files and load data -----
library(tidyverse)
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'
regtable_path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/Regression Tables/Duplications'
#load(paste(path,'illumina.25january.RData',sep = '/')) # cleaned chip dataset
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset
load(paste(path,'pTSpHI.cag.annotated.26feb.RData',sep = '/')) # the LOF and Burden CNV data

# Merge together the LOF dataset with cleaned chip subjects and PNC CNB 
adataset <- merge(cogdata,subject.annotated,by='cag_id')

adataset <- adataset %>% 
  mutate(sex=as.factor(sex),
         race2=as.factor(race2),
         pHI_0=ifelse(pHI>0,1,0)) %>% 
  mutate(pHI_0=as.factor(pHI_0))

covars <- c('envSES','sex2','race2 [2]','race2 [3]')


names(adataset)[31:36] <- substr(names(adataset[31:36]),4,100)
cog.list <- names(adataset[c(28:29,31:36)])



# Table with results for logistic regression model for psychosis spectrum, ADHD, anxiety, depression. 
# (Take out ASD for now because it’s too complicated I think – I’ll discuss with other faculty about whether to include or not)
better.psych.info <- read.csv(paste(path,'pncdxcagid2.csv',sep = '/'))

cnv.psych <- merge(adataset,better.psych.info,by='cag_id')


cnv.psych <- cnv.psych %>% 
  mutate(Psychosis_Spectrum=as.factor(Psychosis_Spectrum),
         ADHD=as.factor(ADHD),
         AnyAnxietyDis=as.factor(AnyAnxietyDis),
         Depression=as.factor(Depression))



write_clip(rbind(
  tidy(glm(Psychosis_Spectrum~pHI+sex+race2+envSES,cnv.psych,family = 'binomial'))[2,],
  tidy(glm(ADHD~pHI+sex+race2+envSES,cnv.psych,family = 'binomial'))[2,],
  tidy(glm(AnyAnxietyDis~pHI+sex+race2+envSES,cnv.psych,family = 'binomial'))[2,],
  tidy(glm(Depression~pHI+sex+race2+envSES,cnv.psych,family = 'binomial'))[2,]))


