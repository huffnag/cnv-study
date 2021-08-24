# Load packages -----
library(tidyverse) # for manipulating data 
library(data.table)
library(clipr) 
library(ggpubr)
library(broom)

# Set path to data files and load data -----
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset

# need to load deletions and duplications together
load(paste(path,'dupes.annotated.RData',sep = '/')) # the LOF and Burden CNV data
names(subject.annotated)[2:length(subject.annotated)] <- paste("dupe.",names(subject.annotated)[2:length(subject.annotated)],sep = '')
duplications <- subject.annotated %>% 
  mutate(dupe.pTS_0=ifelse(dupe.pTS>0,1,0)) %>% 
  mutate(dupe.pTS_0=as.factor(dupe.pTS_0))

load(paste(path,'cag.annotated.3march.RData',sep = '/')) # the LOF and Burden CNV data
names(subject.annotated)[2:length(subject.annotated)] <- paste("delet.",names(subject.annotated)[2:length(subject.annotated)],sep = '')
deletions <- subject.annotated %>% 
  select(!delet.nmarkers) %>% 
  mutate(delet.pHI_0=ifelse(delet.pHI>0,1,0)) %>% 
  mutate(delet.pHI_0=as.factor(delet.pHI_0))

subject.annotated <- merge(deletions,duplications,by = 'cag_id')

# Merge together the LOF dataset with cleaned chip subjects and PNC CNB 
adataset <- merge(cogdata,subject.annotated,by='cag_id')

adataset <- adataset %>% 
  mutate(sex=as.factor(sex),
         race2=as.factor(race2),
         Trauma=as.numeric(traumaExposure))
covars <- c('envSES','sex2','race2 [2]','race2 [3]','Trauma')

# Iterate regression table outputs
model_func = function(response,input1,input2,covars=FALSE) {
  if(class(adataset[,input1])=='numeric'){
    
    if(covars==TRUE){
      form <- paste(response, "~", "scale(", input1,")","+","scale(", input2,")","+envSES+sex+race2+Trauma")
      cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
            tidy(conf.int=T,lm(as.formula(form),
                               data = adataset))[,c(1:3,5:7)],response)}
    
    else {
      form <- paste(response, "~", "scale(", input1,")","+","scale(", input2,")")
      cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
            tidy(conf.int=T,lm(as.formula(form),
                               data = adataset))[,c(1:3,5:7)],response)}}
  
  else {
    
    if(covars==TRUE){
      form <- paste(response, "~", input1,"+",input2,"+envSES+sex+race2+Trauma")
      cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
            tidy(conf.int=T,lm(as.formula(form),
                               data = adataset))[,c(1:3,5:7)],response)}
    
    else {
      form <- paste(response, "~", input1,"+",input2)
      cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
            tidy(conf.int=T,lm(as.formula(form),
                               data = adataset))[,c(1:3,5:7)],response)}}}

input.dels <- names(adataset[grep(pattern = 'del.',names(adataset))])

input.dupes <- names(adataset[grep(pattern = 'dup.',names(adataset))])

drop.list <- c('(Intercept)','envSES',"sex2","race22","race23",'Trauma')

numer_mods <- rbindlist(map2(.x = input.dels, .y=input.dupes,
                             ~model_func(input1 = .x,input2 = .y, response = "Overall_Accuracy",covars= TRUE)))

reg.result <- numer_mods %>% 
  filter(!term %in% drop.list) %>% 
  select(term,estimate,conf.high,conf.low,p.value,AIC) %>% 
  mutate(conf.low=round(conf.low,3),
         conf.high=round(conf.high,3)) %>% 
  mutate(CI=as.character(paste(conf.low,'->',conf.high))) %>% 
  select(!conf.low & !conf.high) %>% 
  select(term,estimate,CI,p.value,AIC) %>% 
  mutate(estimate=round(estimate,3))

write_clip(reg.result)

