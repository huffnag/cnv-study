# Supplemental code for CNV-STUDY

# Load packages
library(tidyverse)
library(data.table)
library(clipr) 
library(ggpubr)
library(broom)
library(wesanderson)

### STable 2 ----
# Set path to data files and load data
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
         race2=as.factor(race2))
covars <- c('envSES','sex2','race2 [2]','race2 [3]')

# Iterate regression table outputs
model_func = function(response,input1,input2,covars=FALSE) {
  if(class(adataset[,input1])=='numeric'){
    
    if(covars==TRUE){
      form <- paste(response, "~", "scale(", input1,")","+","scale(", input2,")","+envSES+sex+race2")
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
      form <- paste(response, "~", input1,"+",input2,"+envSES+sex+race2")
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

drop.list <- c('(Intercept)','envSES',"sex2","race22","race23")

numer_mods <- rbindlist(map2(.x = input.dels, .y=input.dupes, ~model_func(input1 = .x,input2 = .y, response = "Overall_Accuracy",covars= FALSE)))

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

### STable 4 ----
# Set path to data files and load data
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset
load(paste(path,'subs.w22q.RData',sep = '/'))

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
         Trauma=as.numeric(traumaExposure)) %>% 
  filter(!cag_id%in% subs.w22q)

covars <- c('envSES','sex2','race2 [2]','race2 [3]')

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

numer_mods2 <- rbindlist(map2(.x = input.dels, .y=input.dupes, ~model_func(input1 = .x,input2 = .y, response = "Overall_Accuracy",covars= TRUE)))

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

### SFIG. 3 ----
# Set path to data files and load data
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset

# need to load deletions and duplications together
load(paste(path,'dupes.annotated.RData',sep = '/')) # the LOF and Burden CNV data
names(subject.annotated)[2:length(subject.annotated)] <- paste("dupe.",names(subject.annotated)[2:length(subject.annotated)],sep = '')
duplications <- subject.annotated %>% 
  mutate(dupe.LOF_1=ifelse(dupe.LOF>1,1,0),
         dupe.LOF_0=ifelse(dupe.LOF>0,1,0),
         dupe.pTS_0=ifelse(dupe.pTS>0,1,0),
         dupe.LOF_gh=ifelse(dupe.LOF>1/.35,1,0)) %>% 
  mutate(dupe.LOF_1=as.factor(dupe.LOF_1),
         dupe.LOF_0=as.factor(dupe.LOF_0),
         dupe.pTS_0=as.factor(dupe.pTS_0),
         dupe.LOF_gh=as.factor(dupe.LOF_gh))

load(paste(path,'cag.annotated.3march.RData',sep = '/')) # the LOF and Burden CNV data
names(subject.annotated)[2:length(subject.annotated)] <- paste("delet.",names(subject.annotated)[2:length(subject.annotated)],sep = '')
deletions <- subject.annotated %>% 
  select(!delet.nmarkers) %>% 
  mutate(delet.LOF_1=ifelse(delet.LOF>1,1,0),
         delet.LOF_0=ifelse(delet.LOF>0,1,0),
         delet.pHI_0=ifelse(delet.pHI>0,1,0),
         delet.LOF_gh=ifelse(delet.LOF>1/.35,1,0)) %>% 
  mutate(delet.LOF_1=as.factor(delet.LOF_1),
         delet.LOF_0=as.factor(delet.LOF_0),
         delet.pHI_0=as.factor(delet.pHI_0),
         delet.LOF_gh=as.factor(delet.LOF_gh))

subject.annotated <- merge(deletions,duplications,by = 'cag_id')

# Merge together the LOF dataset with cleaned chip subjects and PNC CNB 
adataset <- merge(cogdata,subject.annotated,by='cag_id')

adataset <- adataset %>% 
  mutate(sex=as.factor(sex),
         race2=as.factor(race2),
         pHI_0=ifelse(delet.pHI>0,1,0)) %>% 
  mutate(pHI_0=as.factor(pHI_0),
         Trauma=as.numeric(traumaExposure))

covars <- c('envSES','sex2','race2 [2]','race2 [3]','Trauma')

names(adataset)[31:36] <- substr(names(adataset[31:36]),4,100)
cog.list <- names(adataset[c(28:29,31:36)])

ttest_fun = function(response,input){
  form <- paste(response, "~ scale(",input,")",'+ envSES+scale(dupe.Burden)+Trauma+race2+sex')
  cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
        tidy(conf.int=T,lm(as.formula(form),
                           data = adataset))[,c(1:3,5:7)],response)}

list <- lapply(cog.list,FUN = ttest_fun,input="delet.Burden")
list <- as.data.frame(do.call(rbind,list))

list <- list %>% 
  filter(!term %in% c('(Intercept)','sex2','race22','race23','envSES','Trauma')) %>% 
  arrange(estimate) %>% 
  mutate(var='Cognitive Outcome')

list$response <- gsub('_', " ", list$response)

list$response <- gsub('Cog', "Cognition", list$response)

list$response <- gsub('Comp', "Complex", list$response)

list$response <- gsub('Exec', "Executive", list$response)

full.list <- list

names(adataset)[c(21:24,41:45)] #<- substr(names(adataset[31:40]),4,100)
corr.list <- names(adataset)[21:24]
bifac.list <- names(adataset)[41:45]

# Model Psych ~ pHI
list <- lapply(corr.list,FUN = ttest_fun,input="delet.Burden")
list <- as.data.frame(do.call(rbind,list))
list <- list %>%
  arrange(estimate) %>% 
  filter(!term %in% c('(Intercept)','sex2','race22','race23','envSES','Trauma')) %>% 
  mutate(var='Correlative Traits')

list$response <- gsub('_ar', "", list$response)

list$response <- gsub('_', " ", list$response)

list$response <- gsub('MOOD', "Mood", list$response)

list$response <- gsub('EXTERNALIZING', "Externalizing", list$response)

list$response <- gsub('PSYCHOSIS', "Psychosis", list$response)

list$response <- gsub('FEAR', "Fear", list$response)

list$response <- gsub('CorrTraits', "", list$response)

full.list <- rbind(full.list,list)


list <- lapply(bifac.list,FUN = ttest_fun,input="delet.Burden")
list <- as.data.frame(do.call(rbind,list))
list <- list %>% 
  arrange(estimate) %>% 
  filter(!term %in% c('(Intercept)','sex2','race22','race23','envSES','Trauma')) %>% 
  mutate(var='Bi-factor')

list$response <- gsub('_ar', "", list$response)

list$response <- gsub('_', " ", list$response)

list$response <- gsub('Bifactor', "", list$response)

list$response <- factor(list$response, levels = list$response)

full.list <- rbind(full.list,list)

full.list <- full.list %>% 
  mutate(fdr=p.adjust(p.value,method = 'fdr')) %>%
  mutate(signif=ifelse(fdr<0.05,'yes','no'),
         var=factor(var,levels = c('Cognitive Outcome','Bi-factor','Correlative Traits')))


cpanel <- ggplot(data=full.list,aes(x=abs(estimate), y=response,alpha=factor(signif),color=term)) +
  scale_alpha_discrete(range=c(0.2, 1),name = "pFDR < 0.05", labels = c("False", "True")) +
  geom_point() +
  geom_vline(xintercept = 0,lty=2) +
  xlab("Standardized Effect (Absolute Value)") +
  ylab("") +
  theme_linedraw()+
  scale_color_discrete(name = "Term",labels=c('Total size deletions','Total size duplications'))+
  facet_wrap(~var,ncol = 1,shrink = T,drop = T,scales = 'free')+
  xlim(0,.25)+
  theme(legend.text = element_text(size = 12),
        text = element_text(size=12))

cpanel



### SFIG. 4 ----
# Set path to data files and load data
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
  mutate(pHI_0=as.factor(pHI_0),
         Trauma=as.numeric(traumaExposure))

covars <- c('envSES','sex2','race2 [2]','race2 [3]','Trauma')

names(adataset)[31:36] <- substr(names(adataset[31:36]),4,100)
cog.list <- names(adataset[c(28:29,31:36)])

adataset <- adataset %>% 
  mutate(fact.traumaExposure = as.factor(ifelse(traumaExposure==0,0,ifelse(traumaExposure>0 & traumaExposure<4,1,2))))

tidy(lm(PSYCHOSIS_CorrTraits_ar~pHI_0*fact.traumaExposure+sex+race2+envSES,adataset))

adataset <- adataset %>% 
  mutate(fact.traumaExposure = as.factor(ifelse(traumaExposure==0,0,ifelse(traumaExposure>0 & traumaExposure<4,1,2))))

adataset <- adataset %>% drop_na(fact.traumaExposure)

jitter <- position_jitterdodge(jitter.width = 0.4, jitter.height = 0.1)

p1 <- ggplot(adataset,aes(fact.traumaExposure,PSYCHOSIS_CorrTraits_ar,color=pHI_0))+
  geom_jitter(alpha=.25,position = jitter)+
  theme_light()+
  ylab('Psychosis CorrTraits')+
  xlab("Number of traumatic exposures")


### SFIG. 5 ----
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'

read_excel_allsheets <- function(filename, tibble = T) { # a function to read all 3 sheets from the PRS excel file
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x}

prs_data <- read_excel_allsheets(paste(path,'GO_PNC_PRS_MatchingIDs.xlsx',sep = '/')) # dataset with PRS scores for the PNC
Residual_PNC_EUR_IQ_PRS <- as.data.frame(prs_data[1])
Residual_PNC_EUR_SCZ_PRS <- as.data.frame(prs_data[2])
Linker <- as.data.frame(prs_data[3])
load(paste(path,'illumina.11feb.RData',sep = '/')) # cleaned chip dataset
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset
load(paste(path,'pTSpHI.cag.annotated.26feb.RData',sep = '/')) # the LOF and Burden CNV data
rm(prs_data)

# Merge together the PRS dataset with cleaned chip subjects and PNC CNB 
prs.in.cnv <- Residual_PNC_EUR_IQ_PRS %>%
  filter(Residual_PNC_EUR_IQ_PRS.MATCHING_CAG_CHIP_ID %in% cnv.data$ChipID) %>% 
  dplyr::rename(ChipID=Residual_PNC_EUR_IQ_PRS.MATCHING_CAG_CHIP_ID) %>% 
  select(3:18)

linked <- cnv.data %>%
  select(ChipID,cag_id) %>%
  filter(ChipID %in% prs.in.cnv$ChipID) %>% 
  distinct(ChipID,.keep_all = T)

newprs.in.cnv <- merge(prs.in.cnv,linked,by='ChipID')

full_datamerge <- merge(newprs.in.cnv,cogdata,by='cag_id') # Get cognition data and merge with prs data

full_datamerge <- merge(full_datamerge,subject.annotated,by='cag_id') # Get LOF and Burden data and merge with prs data

full_datamerge.euro <- full_datamerge %>%
  mutate(sex=as.factor(sex),
         race2=as.factor(race2))

scz.in.cnv <- Residual_PNC_EUR_SCZ_PRS %>%
  filter(Residual_PNC_EUR_SCZ_PRS.MATCHING_CAG_CHIP_ID %in% cnv.data$ChipID) %>%
  dplyr::rename(ChipID=Residual_PNC_EUR_SCZ_PRS.MATCHING_CAG_CHIP_ID) %>%
  select(3:18)

full_datamerge.euro <- merge(full_datamerge.euro,scz.in.cnv,by='ChipID') # Get LOF and Burden data and merge with prs data

pcs <- grep('.PC',names(full_datamerge.euro))

termlist <- names(full_datamerge.euro)[pcs] # all the PC's

adataset <- full_datamerge.euro %>% 
  mutate(Trauma=as.numeric(traumaExposure))

covars <- c('envSES','sex2','race2 [2]','race2 [3]')

names(adataset)[47:52] <- substr(names(adataset[47:52]),4,100)

cog.list <- names(adataset[c(44:45,47:52)])
names(adataset)[16] <- "PRS.IQ"
names(adataset)[88] <- "PRS.SCZ"
psy.list <- names(adataset[57:61])

adataset$tip_3tile <- ntile(adataset$PRS.IQ, 3)

head(adataset$tip_3tile)

x <- adataset$PRS.IQ

adataset$PRS.IQgroup <-
  case_when(x > mean(x)+sd(x) ~ "SD Above",
            x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "Average",
            x < mean(x)-sd(x) ~ "SD Below")

count(adataset, PRS.IQgroup)

adataset %>% 
  ggplot() +
  aes(x = pHI, y = Overall_Accuracy, group = PRS.IQgroup, color = PRS.IQgroup) +
  geom_point(color = "grey", alpha = .7) +
  geom_smooth(method = "lm")+
  labs(colour="PRS.IQ Score")+
  ylab('Overall Accuracy')

### SFIG. 6 ----
path <- "/Users/huffnaglen/"

brainlof <- read.csv(paste(path,'brainlof.csv',sep = '/'))

brainlof <- brainlof[,-(2:3)]

brain.deviance.cont <- apply(abs(brainlof[,c('GMV','WMV', 'sGMV')] -.5),1,sum)
brain.deviance <- brainlof[,c('GMV','WMV', 'sGMV')]

sup.thresh <- .9
supranormal <- brain.deviance > sup.thresh
supranormal.all <- apply(supranormal,1,any)

inf.thresh <- .1
infranormal <- brain.deviance < inf.thresh
infranormal.all <- apply(infranormal,1,any)


deviat.all <- infranormal.all | supranormal.all

lof.cat <- rep('low_risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$pHI>(1) | brainlof$pTS>(1)] <- 'high_risk'

toplot.raw <- table(lof.cat,deviat.all)
chisq.test(toplot.raw)

toplot <- table(lof.cat,deviat.all)
toplot[1,] <- toplot[1,]/sum(toplot[1,])
toplot[2,] <- toplot[2,]/sum(toplot[2,])

toplotdf <- data.frame(proportion=c(toplot[,1],toplot[,2]), CNV_Intolerance_Score=c(rownames(toplot), rownames(toplot)), 
                       growthcart=c('low_deviance', 'low_deviance','high_deviance','high_deviance'))

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=proportion, fill=growthcart)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

lab <- c('','n=919','n=60','')
toplotdf$lab <- lab
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))

p+
  geom_text(data=toplotdf,aes(x=CNV_Intolerance_Score,y=proportion,label=lab),vjust=-.4,hjust=1)+
  scale_fill_manual(values=c('#999999','#E69F00'))+ 
  ggtitle("CNV risk and brain deviance")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab('Proportion')+
  xlab('CNV Intolerance Score')+
  theme(text=element_text(size=16))

ggsave('~/Documents/pncCNV/figures/braindev_main.pdf')

### SFIG. 7 ----
###supp plots by phenotype

lof.cat <- rep('low_risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$pHI>(1) | brainlof$pTS>(1)] <- 'high_risk'

toplotdf <- data.frame(CNV_Intolerance_Score=lof.cat, GMV_quantile=brainlof$GMV,  WMV_quantile=brainlof$WMV, sGMV_quantile=brainlof$sGMV)

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=GMV_quantile)) +
  geom_violin()+geom_jitter() + theme_minimal()
p
ggsave('~/Documents/pncCNV/figures/braindev_suppGMV.pdf')

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=WMV_quantile)) +
  geom_violin()+geom_jitter() + theme_minimal()
p
ggsave('~/Documents/pncCNV/figures/braindev_suppWMV.pdf')

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=sGMV_quantile)) +
  geom_violin()+geom_jitter() + theme_minimal()
p
ggsave('~/Documents/pncCNV/figures/braindev_suppSGMV.pdf')

# Use brewer color palettes
##redo with LOF
lof.cat <- rep('low_risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$LOF>(1/.34)] <- 'high_risk'

#lof.cat <- rep('0', length(brainlof$LOF))
#lof.cat[brainlof$LOF>0] <- '0-3'
#lof.cat[brainlof$LOF>(1/.34)] <- '3+'

toplot.raw <- table(lof.cat,deviat.all)
chisq.test(toplot.raw)

toplot <- table(lof.cat,deviat.all)
toplot[1,] <- toplot[1,]/sum(toplot[1,])
toplot[2,] <- toplot[2,]/sum(toplot[2,])

toplotdf <- data.frame(proportion=c(toplot[,1],toplot[,2]), CNV_Intolerance_Score=c(rownames(toplot), rownames(toplot)), 
                       growthcart=c('low_deviance', 'low_deviance','high_deviance','high_deviance'))

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=proportion, fill=growthcart)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

lab <- c('','n=941','n=38','')
toplotdf$lab <- lab
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))

p+geom_text(data=toplotdf,aes(x=CNV_Intolerance_Score,y=proportion,label=lab),vjust=-1)+ scale_fill_manual(values=c('#999999','#E69F00')) +
  ggtitle("CNV risk and brain deviance \n (LOEUF annotation)")  + theme(plot.title = element_text(hjust = 0.5))
ggsave('~/Documents/pncCNV/figures/braindev_supp1.pdf')

##redo with LOF cont

lof.cat <- rep('C_low_risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$LOF>0] <- 'B_medium_risk'
lof.cat[brainlof$LOF>(1/.34)] <- 'A_high_risk'

toplot.raw <- table(lof.cat,deviat.all)
chisq.test(toplot.raw)
toplot <- table(lof.cat,deviat.all)
toplot[1,] <- toplot[1,]/sum(toplot[1,])
toplot[2,] <- toplot[2,]/sum(toplot[2,])
toplot[3,] <- toplot[3,]/sum(toplot[3,])

toplotdf <- data.frame(proportion=c(toplot[,1],toplot[,2]), CNV_Intolerance_Score=c(rownames(toplot), 
                                                                                    rownames(toplot)), growthcart=c('low_deviance', 'low_deviance','low_deviance', 
                                                                                                                    'high_deviance','high_deviance', 'high_deviance'))

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=proportion, fill=growthcart)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

lab <- c('','n=217','n=724','n=38','', '')
toplotdf$lab <- lab
# Use custom colors

p+geom_text(data=toplotdf,aes(x=CNV_Intolerance_Score,y=proportion,label=lab),vjust=-1)+ scale_fill_manual(values=c('#999999','#E69F00')) +
  ggtitle("CNV risk and brain deviance \n (LOEUF annotation with medium risk)")  + theme(plot.title = element_text(hjust = 0.5))
ggsave('~/Documents/pncCNV/figures/braindev_supp2.pdf')

### SFIG. 8 ----
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'

read_excel_allsheets <- function(filename, tibble = T) { # a function to read all 3 sheets from the PRS excel file
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x}

prs_data <- read_excel_allsheets(paste(path,'GO_PNC_PRS_MatchingIDs.xlsx',sep = '/')) # dataset with PRS scores for the PNC
Residual_PNC_EUR_IQ_PRS <- as.data.frame(prs_data[1])
Residual_PNC_EUR_SCZ_PRS <- as.data.frame(prs_data[2])
Linker <- as.data.frame(prs_data[3])
load(paste(path,'illumina.11feb.RData',sep = '/')) # cleaned chip dataset
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset
load(paste(path,'pTSpHI.cag.annotated.26feb.RData',sep = '/')) # the LOF and Burden CNV data
rm(prs_data)

# Merge together the PRS dataset with cleaned chip subjects and PNC CNB 
prs.in.cnv <- Residual_PNC_EUR_IQ_PRS %>%
  filter(Residual_PNC_EUR_IQ_PRS.MATCHING_CAG_CHIP_ID %in% cnv.data$ChipID) %>% 
  dplyr::rename(ChipID=Residual_PNC_EUR_IQ_PRS.MATCHING_CAG_CHIP_ID) %>% 
  select(3:18)

linked <- cnv.data %>%
  select(ChipID,cag_id) %>%
  filter(ChipID %in% prs.in.cnv$ChipID) %>% 
  distinct(ChipID,.keep_all = T)

newprs.in.cnv <- merge(prs.in.cnv,linked,by='ChipID')

full_datamerge <- merge(newprs.in.cnv,cogdata,by='cag_id') # Get cognition data and merge with prs data

full_datamerge <- merge(full_datamerge,subject.annotated,by='cag_id') # Get LOF and Burden data and merge with prs data

full_datamerge.euro <- full_datamerge %>%
  mutate(sex=as.factor(sex),
         race2=as.factor(race2))

scz.in.cnv <- Residual_PNC_EUR_SCZ_PRS %>%
  filter(Residual_PNC_EUR_SCZ_PRS.MATCHING_CAG_CHIP_ID %in% cnv.data$ChipID) %>% 
  dplyr::rename(ChipID=Residual_PNC_EUR_SCZ_PRS.MATCHING_CAG_CHIP_ID) %>% 
  select(3:18)

full_datamerge.euro <- merge(full_datamerge.euro,scz.in.cnv,by='ChipID') # Get LOF and Burden data and merge with prs data

pcs <- grep('.PC',names(full_datamerge.euro))

termlist <- names(full_datamerge.euro)[pcs] # all the PC's

adataset <- full_datamerge.euro %>% 
  mutate(Trauma=as.numeric(traumaExposure))

covars <- c('envSES','sex2','race2 [2]','race2 [3]')

names(adataset)[47:52] <- substr(names(adataset[47:52]),4,100)

cog.list <- names(adataset[c(44:45,47:52)])
names(adataset)[16] <- "PRS.IQ"
names(adataset)[88] <- "PRS.SCZ"
psy.list <- names(adataset[57:61])

#prs.iq with covar
ttest_fun = function(response,input) {
  form <- paste(response, "~ scale(",input,")+scale(Trauma)+scale(pHI)+scale(pTS)+scale(PRS.SCZ)+scale(envSES) + race2 + sex",'+',paste(termlist[1:10],collapse = '+'))
  cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
        tidy(conf.int=T,lm(as.formula(form),
                           data = adataset))[,c(1:3,5:7)],response)}

list <- lapply(cog.list,FUN = ttest_fun,input="PRS.IQ")
list <- as.data.frame(do.call(rbind,list))
list <- list %>% arrange(estimate)

list$response <- gsub('_', " ", list$response)

list$response <- gsub('Cog', "Cognition", list$response)

list$response <- gsub('Comp', "Complex", list$response)

list$response <- gsub('Exec', "Executive", list$response)

vars <- c('scale(envSES)',"scale(Trauma)","scale(PRS.SCZ)",'scale(PRS.IQ)','scale(pHI)','scale(pTS)','scale(Trauma)')

list <- list %>% 
  mutate(fdr=p.adjust(p.value,method = 'fdr')) %>% 
  mutate(signif=ifelse(fdr<0.05,'yes','no')) %>% 
  filter(term %in% vars)

listbind <- list %>% 
  mutate(va='Models of cognition')

iq <- ggplot(data=list,aes(x=abs(estimate), y=response,alpha=factor(signif),color=term)) +
  scale_alpha_discrete(range=c(0.35, 1),name = "pFDR < 0.05", labels = c("False", "True")) +
  geom_point() +
  geom_vline(xintercept = 0,lty=2) +
  xlab("Standardized Effect") +
  ylab("") +
  scale_color_discrete(name = "Term")+
  theme_linedraw()

#prs.scz with covar
list <- lapply(psy.list,FUN = ttest_fun,input="PRS.IQ")
list <- as.data.frame(do.call(rbind,list))
list <- list %>% arrange(estimate)

list <- list %>% 
  mutate(fdr=p.adjust(p.value,method = 'fdr')) %>% 
  mutate(signif=ifelse(fdr<0.05,'yes','no')) %>% 
  filter(term %in% vars)

list$response <- gsub('_ar', "", list$response)

list$response <- gsub('_', " ", list$response)

list$response <- gsub('MOOD', "Mood", list$response)

list$response <- gsub('EXTERNALIZING', "Externalizing", list$response)

list$response <- gsub('PSYCHOSIS', "Psychosis", list$response)

list$response <- gsub('FEAR', "Fear", list$response)

list$response <- factor(list$response, levels = list$response)

list$response <- factor(list$response, levels = list$response)

listbind2 <- list %>% 
  mutate(va='Models of psychopathology')

scz <- ggplot(data=list,aes(x=abs(estimate), y=response,alpha=factor(signif),color=term)) +
  scale_alpha_discrete(range=c(0.35, 1),name = "pFDR < 0.05", labels = c("False", "True")) +
  geom_point() +
  geom_vline(xintercept = 0,lty=2) +
  xlab("Standardized Effect (Absolute Value)") +
  ylab("") +
  scale_color_discrete(name = "Term")+
  theme_linedraw()+
  xlim(0,.3)

full <- rbind(listbind,listbind2)

ggplot(data=full,aes(x=estimate, y=response,alpha=factor(signif),color=term)) +
  scale_alpha_discrete(range=c(0.25, 1),name = "pFDR < 0.05", labels = c("False", "True")) +
  geom_jitter(height = .25) +
  geom_vline(xintercept = 0,lty=2) +
  xlab("Standardized Effect") +
  ylab("") +
  scale_color_discrete(name = "Term",labels=c('SES','pHI','PRS-G','PRS-SCZ','pTS','Trauma'))+
  theme_linedraw()+
  xlim(-.1,.4)+
  facet_wrap(~va,ncol = 1,shrink = T,drop = T,scales = 'free')+
  theme(text = element_text(size=12))
