library(gridExtra)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(tidyverse)

# Set path to data files and load data -----
path <- '/Users/aa2227/Documents/pncCNV/clean'
figpath <- paste(path,'figs',sep='/')
regtable_path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/Regression Tables/Duplications' ####FLAG
load(paste(path,'pnc_cnb_data.RData',sep = '/')) # the PNC cognitive neurobehavioral battery dataset
load(paste(path,'pTSpHI.cag.annotated.26feb.RData',sep = '/')) # the LOF and Burden CNV data
load(paste(path,'dup.cluster.RData',sep = '/'))

# Merge together the LOF dataset with cleaned chip subjects and PNC CNB
adataset <- merge(cogdata,subject.annotated,by='cag_id')

adataset <- adataset %>% 
  mutate(sex=as.factor(sex),
         race2=as.factor(race2),
         pHI_0=ifelse(pHI>0,1,0)) %>% 
  mutate(pHI_0=as.factor(pHI_0))

covars <- c('envSES','sex2','race2 [2]','race2 [3]')

hist.set <- adataset %>% 
  select(pHI,pTS)

histo <- ggplot(hist.set)+
  geom_histogram(data=hist.set %>% select(pHI) %>% mutate(group='pHI'),aes(pHI,fill=group),bins=250,alpha = 0.5)+
  geom_histogram(data=hist.set %>% select(pTS) %>% mutate(group='pTS'),aes(pTS,fill=group), bins=250,alpha = 0.5 )+
  scale_y_log10()+
  ylab('# subjects (log scale)')+
  xlab('pHI/pTS score')+
  scale_fill_manual(name="Score:", values=c("red","#3B9AB2"),labels=c("pHI","pTS"))+
  theme(text = element_text(size=16))+
  theme(legend.position = c(.5, .85),
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.background = element_blank())

dup.cluster <- ggplot(cluster.frame %>% filter(csize<=20),aes(csize,component.mean.pli))+ # mean for cluster
  geom_hex(bins=55)  +
  scale_fill_gradientn(colors =  wes_palette('Zissou1'),
                       guide = 'legend' ,breaks=c(1,10,100,1000,3000),trans='log')+
  xlab('Duplication frequency')+
  ylab('pTS mean')+
  xlim(0,20)+
  geom_label_repel(size=4,data = cluster.frame %>% filter(component.mean.pli>10), 
                   aes(label=paste('Chr',X.chr,': ',spot,sep = '')),nudge_x = 3.5)+
  scale_x_continuous(breaks = round(seq(min(cluster.frame$csize), 20, by = 2),1))+
  theme(text = element_text(size=16))+
  theme(legend.position = 'none')

load(paste(path,'del.cluster.RData',sep = '/'))

del.cluster <- ggplot(cluster.frame %>% filter(csize<=20),aes(csize,component.mean.pli))+ # mean for cluster
  geom_hex(bins=55)  +
  scale_fill_gradientn(colors =  wes_palette('Zissou1'),
                       guide = 'legend' ,breaks=c(1,10,100,1000,3000),trans='log')+
  xlab('Deletion frequency')+
  ylab('pHI mean')+
  xlim(0,20)+
  geom_label_repel(size=4,data = cluster.frame %>% filter(component.mean.pli>10), 
                   aes(label=paste('Chr',X.chr,': ',spot,sep = '')),nudge_x = 3.5)+
  scale_x_continuous(breaks = round(seq(min(cluster.frame$csize), 20, by = 2),1))+
  ylim(-.01,25)+
  theme(text = element_text(size=16))+
  theme(legend.position = c(.5, .95),
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.background = element_blank())

pdf(file = paste(figpath,'fig1.pdf',sep = '/'),width = 12,height = 6)

ggarrange(del.cluster, dup.cluster,histo, nrow=1,labels = c("A)", "B)",'C)'))

dev.off()
