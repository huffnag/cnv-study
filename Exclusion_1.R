# Load packages
library(tidyverse)

# Load data
path <- "/Users/huffnaglen/PNC CNV Project copy/Analysis/RData"
load(paste(path,'chip.qc.cag.RData',sep = '/')) # this is illumina chip dataset with all chip level qc
load(paste(path,'full.cnvision_merge.RData',sep = '/'))
chip.markers <- read.csv(paste(path,'chipmarkers.csv',sep = '/'))

ncnvs <- as.data.frame(table(dtcnv$SampleID))

names(ncnvs)[1] <- 'chipids'

newmerge <- merge(chip.markers,ncnvs,by = 'chipids')

cleanmerge <- newmerge %>%
  filter(Freq<50 | nmarkers > 1000000) %>%
  filter(Freq<200)

cleanmerge$ChipID <- unlist(strsplit(cleanmerge$chipids,'.int.csv'))

# 1. count raw
length(unique(chip.qc.cag$ChipID)) # 10622
length(unique(chip.qc.cag$cag_id)) # 8067

# Remove GH/Seb chips
chip.qc.cag <- chip.qc.cag %>% 
  filter(ChipID %in% cleanmerge$ChipID)

length(unique(chip.qc.cag$ChipID)) # 9483
length(unique(chip.qc.cag$cag_id)) # 7697

# 2. exclude baf, lrr, wf
chip.qc.cag <- chip.qc.cag %>% 
  filter(lrrsd<0.35)
length(unique(chip.qc.cag$ChipID)) # 9458
length(unique(chip.qc.cag$cag_id)) # 7679

chip.qc.cag <- chip.qc.cag %>% 
  filter(bafsd<0.08)
length(unique(chip.qc.cag$ChipID)) # 9458
length(unique(chip.qc.cag$cag_id)) # 7679

chip.qc.cag <- chip.qc.cag %>% 
  filter(abs(wavef)<0.05)
length(unique(chip.qc.cag$ChipID)) # 9388
length(unique(chip.qc.cag$cag_id)) # 7671

list.a <- chip.qc.cag %>% 
  select(1,4:7)

save(list.a,file = paste(path,'list.a2.RData',sep = '/'))

# 5. Look at cnv level w/ these chips
cnv.set <- dtcnv %>% 
  filter(ChipID %in% chip.qc.cag$ChipID)

# 6. Algorithm cutoffs
cnv.set$two.algs <- as.numeric(unlist(strsplit(cnv.set$X.Two.Algs,"%")))

# Keep cnvs w/ confidence score >= 30
cnv.set <- cnv.set %>% 
  filter(Max.Score>=30)
length(unique(cnv.set$ChipID)) # 8997
length(unique(cnv.set$cag_id)) # 7610

# Keep cnvs w/ both algorithm
cnv.set <- cnv.set %>%
  filter(X.Algos==2)
length(unique(cnv.set$ChipID)) # 8880
length(unique(cnv.set$cag_id)) # 7578

# Keep cnvs w/ 70% algorithm overlap
cnv.set <- cnv.set %>% 
  filter(two.algs>=70)
length(unique(cnv.set$ChipID)) # 8761
length(unique(cnv.set$cag_id)) # 7541

list.b <- chip.qc.cag %>% 
  filter(ChipID %in% cnv.set$ChipID)

count.set <- as.data.frame(table(cnv.set$ChipID))
names(count.set) <- c("ChipID","Freq")

null.set <- as.data.frame(setdiff(list.a$ChipID,count.set$Var1))
null.set[,2] <- 0
names(null.set) <- c("ChipID","Freq")

qspncv.set <- rbind(count.set,null.set)

save(cnv.set,file = paste(path,'cnv.set.11feb.RData',sep = '/'))
