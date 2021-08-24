#analyze quanttraits
library(tidyverse)

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

