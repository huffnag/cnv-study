# Exclusion_1 script
# Takes data from CNVision merge format and applies chip quality control and markers criteria (detailed in SFigs. 1-2).

# Load packages
library(tidyverse)

# Load data
path <- "/Users/huffnaglen/PNC CNV Project copy/Analysis/RData"
load(paste(path,'chip.qc.cag.RData',sep = '/'))
load(paste(path,'full.cnvision_merge.RData',sep = '/'))
chip.markers <- read.csv(paste(path,'chipmarkers.csv',sep = '/'))

ncnvs <- as.data.frame(table(dtcnv$SampleID))

names(ncnvs)[1] <- 'chipids'

newmerge <- merge(chip.markers,ncnvs,by = 'chipids')

# Markers - CNV criteria
cleanmerge <- newmerge %>%
  filter(Freq<50 | nmarkers > 1000000) %>%
  filter(Freq<200)

cleanmerge$ChipID <- unlist(strsplit(cleanmerge$chipids,'.int.csv'))

# Remove GH/Seb chips
chip.qc.cag <- chip.qc.cag %>% 
  filter(ChipID %in% cleanmerge$ChipID)

# Exclude baf, lrr, wf
chip.qc.cag <- chip.qc.cag %>% 
  filter(lrrsd<0.35)

chip.qc.cag <- chip.qc.cag %>% 
  filter(bafsd<0.08)

chip.qc.cag <- chip.qc.cag %>% 
  filter(abs(wavef)<0.05)

list.a <- chip.qc.cag %>% # Save a simplified cnv quality dataset for later use.
  select(1,4:7)
save(list.a,file = paste(path,'list.a2.RData',sep = '/'))

# Look at cnv level w/ these good quality chips
cnv.set <- dtcnv %>% 
  filter(ChipID %in% chip.qc.cag$ChipID)

# Algorithm cutoffs
cnv.set$two.algs <- as.numeric(unlist(strsplit(cnv.set$X.Two.Algs,"%")))

# Keep cnvs w/ confidence score >= 30
cnv.set <- cnv.set %>% 
  filter(Max.Score>=30)

# Keep cnvs w/ both algorithm
cnv.set <- cnv.set %>%
  filter(X.Algos==2)

# Keep cnvs w/ 70% algorithm overlap
cnv.set <- cnv.set %>% 
  filter(two.algs>=70)

list.b <- chip.qc.cag %>% 
  filter(ChipID %in% cnv.set$ChipID)

count.set <- as.data.frame(table(cnv.set$ChipID))
names(count.set) <- c("ChipID","Freq")

null.set <- as.data.frame(setdiff(list.a$ChipID,count.set$Var1))
null.set[,2] <- 0
names(null.set) <- c("ChipID","Freq")

qspncv.set <- rbind(count.set,null.set)

save(cnv.set,file = paste(path,'cnv.set.11feb.RData',sep = '/'))
