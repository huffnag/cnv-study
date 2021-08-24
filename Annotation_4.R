# for each cnv, check for 50% overlap with another cnv 
# if it has overlap, then count those as the same cnv
#########################################
# Directory ----------------------------- 
path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'

# Packages ----------------------------- 
library(tidyverse)
library(IRanges)
library(data.table)
library(igraph)
library(sjPlot)
library(effects)


path <- '/Users/huffnaglen/PNC CNV Project copy/Analysis/RData'

load(paste(path,'dupe.overlaps.data.23june2021.RData',sep = '/'))

# cnv cluster membership is releative to chromosome now
# would now be good to have nother column with the sum of 
# the size of the chr relative cluster 

#---- Get cluster sizes per chromosomes
cluster.sizes <- final.bind %>% 
  group_by(X.chr,membership) %>% 
  mutate(csize = length(membership))

#---- pLI scores
load(paste(path,'ploi.RData',sep = '/'))
phi <- read.csv(paste(path,'phi.csv',sep = '/'))
rm(thing,cnv_dels.dups,cogdata)
# Calculate oeuf Scores for Individual CNV's -----
ploi <- ploi %>% 
  filter(chromosome!="X" & chromosome!="Y") %>% 
  mutate(chromosome=as.numeric(chromosome))

length(setdiff(ploi$gene,phi$gene))
# 1883 genes are in gnomad but not phi

with.phi <- merge(ploi,phi,by='gene')

no.phi <- ploi %>% 
  filter(gene %in% setdiff(ploi$gene,phi$gene)) %>% 
  mutate(pHI=as.numeric(0),pTS=as.numeric(0))

ploi <- rbind(with.phi,no.phi)

#---- Calculate pLI Scores for Individual CNV's and CNV Components
ploi <- ploi %>% 
  filter(chromosome!="X" & chromosome!="Y") %>% 
  mutate(chromosome=as.numeric(chromosome))

pli.data <- data.frame()

cluster.sizes[,ncol(cluster.sizes)+1] <- 0 # initialize an empty col for pLI

cluster.sizes <- ungroup(cluster.sizes)

for(chrome in 1:22){
  print(chrome)
  
  lengths <- cluster.sizes %>% # get deletion sizes  
    filter(X.chr==chrome) %>% 
    select(start.hg19.,end.hg19.) %>%
    mutate(dif = end.hg19. - start.hg19.)
  
  ploi.loop <- ploi %>% # filter values so that only genes which fall completely within a cnv are kept
    filter(chromosome==chrome) %>%
    filter(end_position <= max(lengths$end.hg19.)) %>%
    filter(start_position >= min(lengths$start.hg19.)) %>%
    select(start_position,end_position,pTS)
  
  dels.only.pli <- cluster.sizes %>% 
    filter(X.chr==chrome)
  
  for(each in 1:nrow(dels.only.pli)){ # for each row in the cnv data
    cnv.location <- c(dels.only.pli[each,4],dels.only.pli[each,5]) # set the start and stop base location of the current cnv
    condition <- ploi.loop[,1] >= cnv.location[1] & ploi.loop[,2] <= cnv.location[2] # test if a gene falls completely within the current cnv row (the deletion space)
    pli.temp <- as.data.frame(matrix(nrow = nrow(ploi.loop),ncol = 1)) # create a temporary genome sized data frame for each cnv loop around
    for(ij in 1:length(condition)){ # go through each row in the gene data
      if(condition[ij]){ # if condition is true (gene is within cnv)
        pli.temp[ij,1] <- ploi.loop[ij,3]}} # then save the pLI value by inserting into the temp df
    dels.only.pli[each,24] <- as.numeric(sum(pli.temp$V1,na.rm = T)) # after all 'true' genes saved get a final summed pLI score for the cnv
  }
  pli.data <- rbind(pli.data,dels.only.pli)
}

names(pli.data)[24] <- 'cnv.specific.pli'

# 
pli.data.grouped <- pli.data %>%
  group_by(membership, X.chr) %>%
  mutate(component.mean.pli = mean(cnv.specific.pli)) %>% 
  mutate(component.median.pli = median(cnv.specific.pli)) %>% 
  mutate(spot=paste(round(mean(start.hg19.)/1000000,1),'-',round(mean(end.hg19.)/1000000,1),sep = ''),
         fullspot=paste(round(mean(start.hg19.),0),'-',round(mean(end.hg19.),0),sep = ''))

cluster.frame <- pli.data.grouped %>% 
  group_by(X.chr,membership) %>% 
  distinct(membership,.keep_all = T)


dup.cluster <- ggplot(cluster.frame %>% filter(csize<=20),aes(csize,component.mean.pli))+ # mean for cluster
  geom_hex(bins=55)  +
  scale_fill_gradientn(colors =  wes_palette('Zissou1'),
                       guide = 'legend' ,breaks=c(1,10,100,1000,3000),trans='log')+
  xlab('cluster size')+
  ylab('cluster pTS mean')+
  xlim(0,20)+
  geom_label_repel(size=4,data = cluster.frame %>% filter(component.mean.pli>10), 
                   aes(label=paste('Chr',X.chr,': ',spot,sep = '')),nudge_x = 3.5)+
  scale_x_continuous(breaks = round(seq(min(cluster.frame$csize), 20, by = 2),1))+
  theme(text = element_text(size=16))+
  theme(legend.position = c(.75, .85),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),legend.background = element_blank())
dup.cluster


save(cluster.frame,file = paste(path,'/dup.cluster.RData',sep = ''))

# cnv table





