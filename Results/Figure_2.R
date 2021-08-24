



# NOTE #
# Run Annotation_1 before this script to get 'pTSpHI.cag.annotated.26feb.RData'

library(broom)
library(wesanderson)

# Set path to data files and load data -----
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

ttest_fun = function(response,input){
  form <- paste(response, "~ scale(",input,")",'+ envSES+scale(pTS)+Trauma+race2+sex')
  cbind(glance(lm(as.formula(form), data = adataset))[c(1:4,6:12)],
        tidy(conf.int=T,lm(as.formula(form),
                           data = adataset))[,c(1:3,5:7)],response)}

list <- lapply(cog.list,FUN = ttest_fun,input="pHI")
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
list <- lapply(corr.list,FUN = ttest_fun,input="pHI")
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

list <- lapply(bifac.list,FUN = ttest_fun,input="pHI")
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
  theme_light()+
  scale_color_discrete(name = "Term", labels = c("pHI", "pTS"))+
  facet_wrap(~var,ncol = 1,shrink = T,drop = T,scales = 'free')+
  xlim(0,.25)+
  theme(legend.text = element_text(size = 12),
        text = element_text(size=12))+
  theme(strip.text = element_text(colour = 'white'))+
  theme(strip.background =element_rect(fill="black"))
  





# Load packages -----
library(tidyverse) # for manipulating data 
library(sjPlot) # for regression model tables
library(data.table) # for turning lists into dataframes
library(ggpubr)
library(gridExtra)
library(broom)
library(clipr)

# Set path to data files and load data -----
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'
regtable_path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/Regression Tables/Duplications'
#load(paste(path,'illumina.25january.RData',sep = '/')) # cleaned chip dataset
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset
load(paste(path,'cag.annotated.3march.RData',sep = '/')) # the LOF and Burden CNV data

# Merge together the LOF dataset with cleaned chip subjects and PNC CNB 
adataset <- merge(cogdata,subject.annotated,by='cag_id')

adataset <- adataset %>% 
  mutate(sex=as.factor(sex),
         race2=as.factor(race2),
         pHI_0=ifelse(pHI>0,1,0)) %>% 
  mutate(pHI_0=as.factor(pHI_0))


apanel <- ggplot(adataset,aes(x=log(pHI+1),y=Overall_Accuracy))+
  geom_hex(bins=50)+
  scale_fill_gradientn(name='# subjects',colors = wes_palette(n=5, name="Zissou1"),trans='log10',guide = 'legend')+
  ylim(-10,9)+
  ylab('Overall Accuracy')+
  xlab('log(pHI+1)')+
  geom_smooth(method = 'lm',linetype='dotted',se=F,color="black")+
  theme_light()+
  theme(legend.position = c(.5, .8),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.title  = element_text(size = 8))

bpanel <- ggplot(adataset,aes(x=log(pHI+1),y=Bifactor_Psychosis_ar))+
  geom_hex(bins=50)+
  scale_fill_gradientn(name='n subjects',colors = wes_palette(n=5, name="Zissou1"),trans='log10',guide = 'legend')+
  ylim(-9,9)+
  ylab('Bifactor Psychosis')+
  xlab('log(pHI+1)')+
  geom_smooth(method = 'lm',linetype='dotted',se=F,color="black")+
  theme_light()+
  theme(legend.position = 'none')

pdf(file = paste(path,'fig2.v2.pdf',sep = '/'),width = 12,height = 6)

grid.arrange(apanel, cpanel,bpanel, nrow=2,ncol=2,
             layout_matrix = rbind(c(1,2), c(3,2)),
             widths = c(1,1.75))



dev.off()
