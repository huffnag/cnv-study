# Exclusion_2
# Takes data from Exclusion_1 and applies remaining CNV cleaning criteria (detailed in SFig. 1).

# Set data path and load packages
path <- "/Users/huffnaglen/PNC CNV Project copy/Analysis/RData"
library(tidyverse)
library(IRanges)
library(igraph)
library(readxl)

# Read Data ----------------------------- 
load(paste(path,'cnv.set.11feb.RData',sep = '/'))

# 0. Get CNVs -----------------------------
cnv.data <- cnv.set

#### 1. Variant Size >= 50 kb -----------------------------
cnv.data <- cnv.data %>% 
  filter(sizedel >= 50000)

### 2. Variant Confidence Score >= 30 -----------------------------
cnv.data$ConfidenceScore <- cnv.data$Max.Score
cnv.data$ConfidenceScore <- as.numeric(cnv.data$ConfidenceScore)
cnv.data <- cnv.data %>%
  filter(ConfidenceScore >= 30)

### 3. Variant has < 50% overlap w/ Segmental Duplications -----------------------------
load(paste(path,'segmental.duplications.RData',sep = '/'))
# fix chromosome labels
data.segments$chrom.str <- substr(data.segments$chrom,1,5)
x <- strsplit(data.segments$chrom.str,'chr')
data.segments$splitter<-sapply(x,"[[",2)

y <- strsplit(data.segments$splitter,'_')
data.segments$splitter2<-sapply(y,"[[",1)

data.segments <- data.segments %>% filter(splitter2!='Un' & splitter2!= 'X' & splitter2!='Y')
table(data.segments$splitter2)
data.segments$splitter2 <- as.numeric(data.segments$splitter2)
table(data.segments$splitter2)

rm(x,y)

# sd overlap analysis
sd.bind <- data.frame()
for(i in 1:22){
  print(i)
  pnc.data <- cnv.data %>% 
    filter(X.chr==i) %>%
    select(start.hg19.,end.hg19.,Chr.Start.Stop.hg19.)
  startcnv <- pnc.data$start.hg19.
  endcnv <- pnc.data$end.hg19.
  pnc.range <- IRanges(start=startcnv, end=endcnv)
  
  sd.dat <- data.segments %>% 
    filter(splitter2==i) %>%
    select(chromStart,chromEnd)
  startsd <- sd.dat$chromStart
  endsd <- sd.dat$chromEnd
  sd.range <- IRanges(start=startsd, end=endsd)
  
  sd.overlaps <- findOverlaps(pnc.range, sd.range)
  
  sd.overlaps <- as.matrix(sd.overlaps)
  
  percentoverlaprec <- vector()
  
  for(f in 1:nrow(sd.overlaps)){ #
    index1 <- sd.overlaps[f,1] #
    index2 <- sd.overlaps[f,2] #
    
    start1<-startcnv[index1] #
    end1<-endcnv[index1] #
    
    start2<-startsd[index2] #
    end2 <-endsd[index2] #
    
    overlap <- min(end1,end2) - max(start1,start2)
    pnc.length <- end1 - start1
    percentoverlaprec <- c(percentoverlaprec,overlap/pnc.length)}
  
  sd.overlaps.trim <- matrix(sd.overlaps[which(percentoverlaprec >= .5),],ncol = 2)
  
  bad.rows.chr1 <- unique(sd.overlaps.trim[,1])
  
  pnc.data.bad <- pnc.data[c(bad.rows.chr1),]
  
  sd.bind <- rbind(pnc.data.bad,sd.bind)}
# end loop

badlist <- sd.bind$Chr.Start.Stop.hg19.

cnv.data <- cnv.data %>% 
  filter(!Chr.Start.Stop.hg19. %in% badlist) # removes cnv that overlap > 50% with SD's

### 4. Remove Centromere overlaps -----------------------------
# remove variants that span the centromere
load(paste(path,'hg19.RData',sep = '/'))

centro.bind <- data.frame() # setup data frame
full.overlap <- vector()
centro.overlap <- vector()

for(abc in 1:22){ # for each chromosome
  pnc.data. <- cnv.data %>% filter(X.chr==abc) %>%  # subset cnv data by chr.
    select(X.chr,start.hg19.,end.hg19.,Chr.Start.Stop.hg19.)
  
  cnv.start <- pnc.data.$start.hg19. # set variant start
  cnv.end <- pnc.data.$end.hg19. # set variant end
  
  cnv.range <- IRanges(start=cnv.start, end=cnv.end) # set ranges of variants
  
  hg19.1 <- hg19 %>% filter(chrom==abc) # subset centromere data by chr.
  
  centro.start <- hg19.1$centromerStart # set centromere start
  centro.end<- hg19.1$centromerEnd # set centromere end
  
  centro.range <- IRanges(start=centro.start, end=centro.end) # set range of centromere
  
  centro.overlaps <- findOverlaps(cnv.range,centro.range) # find overlaps between variants and the centromere
  
  centro.overlaps <- as.matrix(centro.overlaps) # as matrix
  
  print(centro.overlaps) # see progress/review
  
  centro.overlap <- vector() # setup vector
  
  for(fx in 1:nrow(centro.overlaps)){ # for each row of all overlaps
    if(nrow(centro.overlaps)<1){next} # if no overlaps, skip loop
    
    index1 <- centro.overlaps[fx,1] # the variant locations
    index2 <- centro.overlaps[fx,2] # the centromere
    
    start1<-cnv.start[index1] # variant start
    end1<-cnv.end[index1] # variant end
    
    start2<-centro.start[index2] # centromere start (always same)
    end2 <-centro.end[index2] # centromere end (always same)
    
    overlap <- min(end1,end2) - max(start1,start2) # the total overlap
    centro.length <- end2 - start2 # the centormere length
    centro.overlap <- c(centro.overlap,overlap/centro.length)} # the proprotion of centromere covered by variant
  
  full.overlap <- c(full.overlap,centro.overlap)
  
  
  centro.overlaps.trim <- matrix(centro.overlaps[which(centro.overlap == 1),],
                                 ncol = 2) # overlaps which cover the entire centromere
  
  bad.rows.chr <- unique(centro.overlaps.trim[,1]) # unique cnvs in data that cover a centromere
  
  pnc.data.bad <- pnc.data.[c(bad.rows.chr),] # the pnc data overlaps
  
  centro.bind <- rbind(pnc.data.bad,centro.bind)} # bind all centromere overlaps together
### end loop ###
###          ###

badlist <- centro.bind$Chr.Start.Stop.hg19.

cnv.data <- cnv.data %>% 
  filter(!Chr.Start.Stop.hg19. %in% badlist) %>%  # removes cnv that overlap > 50% with SD's
  arrange(X.chr)

### 5. Remove Telomere overlaps -----------------------------
end.telomere.1 <- 100*1000 # always 100k for each chromosme
telomere.hits <- data.frame()
for(each in 1:22){
  pnc.chr <- cnv.data %>% filter(X.chr==each) %>% 
    select(X.chr,start.hg19.,end.hg19.,Chr.Start.Stop.hg19.)
  
  hg.chr <- hg19 %>% filter(chrom==each)
  
  start.telomere2 <- hg.chr$length - (100*1000)
  
  pnc.chr.hits <- pnc.chr %>% 
    filter(start.hg19. < end.telomere.1 | end.hg19. > start.telomere2)
  
  telomere.hits <- rbind(telomere.hits,pnc.chr.hits)
  
  print(nrow(pnc.chr.hits))}

cnv.data <- cnv.data %>%
  filter(!Chr.Start.Stop.hg19. %in% telomere.hits$Chr.Start.Stop.hg19.)

### 6. Variant has < 50% overlap with major histocompatibility complex  -----------------------------
pnc.data.clean.hmc <- cnv.data %>% filter(X.chr==6)

pnc6.start <- pnc.data.clean.hmc$start.hg19.
pnc6.end <- pnc.data.clean.hmc$end.hg19.

pnc6.ranges <- IRanges(pnc6.start,pnc6.end)

hmc.start <- 28477797
hmc.end <- 33448354

hmc.range <- IRanges(hmc.start,hmc.end)

hmc.overlaps <- findOverlaps(pnc6.ranges,hmc.range)

hmc.overlaps <- as.matrix(hmc.overlaps)

hmc.overlap.percent <- vector()

for(f in 1:nrow(hmc.overlaps)){ #
  index1 <- hmc.overlaps[f,1] #
  index2 <- hmc.overlaps[f,2] #
  
  start1<-pnc6.start[index1] #
  end1<-pnc6.end[index1] #
  
  start2<-hmc.start[index2] #
  end2 <-hmc.end[index2] #
  
  overlap <- min(end1,end2) - max(start1,start2)
  pnc.length <- end1 - start1
  hmc.overlap.percent <- c(hmc.overlap.percent,overlap/pnc.length)}

hmc.overlaps.trim <- matrix(hmc.overlaps[which(hmc.overlap.percent >= .5),],ncol = 2)

bad.rows.chr1 <- unique(hmc.overlaps.trim[,1])

pnc.data.bad <- pnc.data.clean.hmc[c(bad.rows.chr1),]

badlist <- pnc.data.bad$Chr.Start.Stop.hg19.

cnv.data <- cnv.data %>% 
  filter(!Chr.Start.Stop.hg19. %in% badlist) # removes cnv that overlap > 50% with SD's

### Chip check
# Make a data frame of chips that are now gone but were of good quality from Exclusion_1
# They have zero cnv and zero burden
load(paste(path,'list.a2.RData',sep = '/'))

list.b <- list.a %>% 
  filter(!ChipID %in% cnv.data$ChipID) %>%
  mutate(Freq=0,Burden=0) %>% 
  select(ChipID,Burden,Freq)

newncnvset <- as.data.frame(table(cnv.data$ChipID))
names(newncnvset)[1] <- c('ChipID')

burdenset <- cnv.data %>% 
  group_by(ChipID) %>% 
  summarise(Burden=sum(sizedel)/1000000)

has.cnv <- merge(burdenset,newncnvset,by='ChipID')

chip.data <- rbind(list.b,has.cnv)

### 7. Chips with >= 30 cnv's excluded -----------------------------
### 8. Chip burden Size <= 8Mb -----------------------------
qual.chip.data <- chip.data %>% 
  filter(Burden<8)

cnv.data <- cnv.data %>% 
  filter(ChipID %in% qual.chip.data$ChipID) 

# Add the cnvless quality chips to cnv-level dataset with null/NA values for cnv variables
load(paste(path,'chip.qc.cag.RData',sep = '/'))
adding.cnv <- chip.qc.cag %>% 
  filter(ChipID %in% setdiff(qual.chip.data$ChipID,cnv.data$ChipID)) %>% 
  select(ChipID,cag_id)
adding.cnv[,3:19] <- NA 
adding.cnv[,17] <- adding.cnv$cag_id
adding.cnv[,2] <- NA
adding.cnv[,16] <- adding.cnv$ChipID
adding.cnv[,1] <- NA
names(adding.cnv) <- names(cnv.data)
cnv.data <- rbind(cnv.data,adding.cnv) # this cnv level data now has the 'empty' chips

save(cnv.data,file = paste(path,'cnv_pre_subjectlevel.RData'))

### 9. Pick lowest lrrsd for subs w/ more than 1 chip
chipppers <- cnv.data %>% 
  group_by(cag_id) %>% 
  summarise(number_distinct=length(unique(ChipID))) # number of distinct chips a subject has

good.subs <- chipppers %>% 
  filter(number_distinct<2) # each of these subjects has only 1 chip, their data is good
good.vect <- as.vector(good.subs$cag_id) # 5953 subjects have only 1 chip
perfect.subs <- cnv.data %>% 
  filter(cag_id %in% good.vect)

bad.subs <- chipppers %>% 
  filter(number_distinct>1) # each of these subjects has more than 1 chip
bad.vect <- as.vector(bad.subs$cag_id) # 502 subjects have more than 1 chip
nonperfect.subs <- cnv.data %>% 
  filter(cag_id %in% bad.vect)

load(paste(path,'chipqualitydf.RData',sep = '/'))

new.chipqualitydf <- chipqualitydf %>% 
  filter(chipids %in% nonperfect.subs$ChipID) %>% 
  select(chipids,LRRSD)

empty.df <- data.frame()
for(i in unique(nonperfect.subs$cag_id)){
  loopdf <- nonperfect.subs %>% 
    filter(cag_id==i) # get a subject, by cag id
  
  loop.quality <- new.chipqualitydf %>% 
    filter(chipids %in% loopdf$ChipID) # get chip quality data for that cag id
  
  qualitysmall2 <- loop.quality %>% 
    filter(LRRSD==min(LRRSD)) # get the data for the min lrrsd
  
  if(nrow(qualitysmall2)>1){
    print('glitch')
    print(qualitysmall2)
    qualitysmall2 <- qualitysmall2[1,]
    print(qualitysmall2)}
  
  smalltest2 <- loopdf %>% 
    filter(ChipID==qualitysmall2$chipids)
  
  empty.df <- rbind(empty.df,smalltest2)}

cnv.data <- rbind(empty.df,perfect.subs)

### 10. Pick random family member per family (see fam id script for methods)

# Read ancestry data and cnv data ----------------------------- 
fam.id <- read_xlsx(paste(path,'PNC_GO_IDs_dbGaP_Related_Ethnicity_Files_Report.xlsx',sep = '/'))

# Strsplit pnc chips -----
lof.chips <- strsplit(cnv.data$ChipID,'.int.csv')
cnv.data$chip <- sapply(lof.chips,'[[',1)

# remove na first?
fam.id <- fam.id %>% na.omit(ID1,ID2)

# get fam id's with chips together -----------------------------
# build a df with col for chip and col for famid
fam.df <- data.frame()
fam.id.list <- unique(fam.id$FID)
for(each in fam.id.list){
  this.fam <- fam.id %>% filter(FID==each)
  split1 <- unlist(strsplit(this.fam$ID1,','))
  split1 <- unique(split1)
  
  split2 <- unlist(strsplit(this.fam$ID2,','))
  split2 <- unique(split2)
  
  full.list <- as.data.frame(c(split1,split2))
  
  this.fam.df <- data.frame(full.list,rep(each,times=length(full.list)))
  names(this.fam.df) <- c('chip','fam.id')
  
  fam.df <- rbind(fam.df,this.fam.df)}

# see how many fam chips in our clean data
check.fams <- cnv.data %>% filter(chip %in% fam.df$chip)

# how do I give each check.fams subject their family id name as a variable?
length(unique(check.fams$cag_id))

fams.distinct <- distinct(fam.df, chip, .keep_all=T)

length(unique(fams.distinct$chip))

fam.ids.pncs <- merge(check.fams,fams.distinct,by='chip')
length(unique(fam.ids.pncs$fam.id))
# 459 fams with fam id's

# no family data (later for binding) -----
lof.nofam <- cnv.data %>% filter(!chip %in% fam.df$chip)

# full sample
all.fams <- vector()

set.seed(9454444)
for(subid in unique(fam.ids.pncs$fam.id)){
  fam.data <- fam.ids.pncs %>% filter(fam.id==subid)
  the.fam.member <- sample(fam.data$cag_id,size = 1)
  all.fams <- c(all.fams,the.fam.member)}

length(all.fams) # 459

length(unique(all.fams))

nrow(fam.ids.pncs) # 2940

length(unique(fam.ids.pncs$cag_id)) # 957

# the sample to keep is the vector allfams
# so %in% from fam.ids.pncs
fam.sampled <- fam.ids.pncs %>% filter(cag_id %in% all.fams)
length(unique(fam.sampled$chip))
length(unique(fam.sampled$cag_id))

fam.sampled <- fam.sampled %>%  # remove fam id for rbinding
  select(!fam.id)

# clean data no fam members
cnv.data <- as.data.frame(rbind(lof.nofam,fam.sampled))

save(cnv.data,file = paste(path,'illumina.11feb.RData',sep = '/'))
