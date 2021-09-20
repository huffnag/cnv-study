#Structural MRI normative model results
library(tidyverse)

path <- '/Users/aa2227/Documents/pncCNV/clean'
figpath <- paste(path,'figs',sep='/')

brainlof <- read.csv(paste(path,'brainlof_qc.csv',sep = '/'))

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

lof.cat <- rep('low risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$pHI>(1) | brainlof$pTS>(1)] <- 'high risk'

toplot.raw <- table(lof.cat,deviat.all)
chisq.test(toplot.raw)

toplot <- table(lof.cat,deviat.all)
toplot[1,] <- toplot[1,]/sum(toplot[1,])
toplot[2,] <- toplot[2,]/sum(toplot[2,])

toplotdf <- data.frame(proportion=c(toplot[,1],toplot[,2]), CNV_Intolerance_Score=c(rownames(toplot), rownames(toplot)), 
  growthcart=c('low deviance', 'low deviance','high deviance','high deviance'))
p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=proportion, fill=growthcart)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

lab <- c('','n=891','n=59','')
toplotdf$lab <- lab
# Use custom colors
p + scale_fill_manual(values=c("red","#3B9AB2"))

p+
  geom_text(data=toplotdf,aes(x=CNV_Intolerance_Score,y=proportion,label=lab),vjust=-.4,hjust=1)+
  scale_fill_manual(values=c("red","#3B9AB2"))+ 
  ggtitle("CNV risk and brain deviance")  + 
  theme(plot.title = element_text(hjust = 0.5))+
  ylab('Proportion')+
  xlab('CNV Intolerance Score')+
  theme(text=element_text(size=16))


ggsave(paste(figpath,'fig4.pdf',sep = '/'))

###Supplemental plots by phenotype

lof.cat <- rep('low risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$pHI>(1) | brainlof$pTS>(1)] <- 'high risk'

toplotdf <- data.frame(CNV_Intolerance_Score=lof.cat, GMV_quantile=brainlof$GMV,  WMV_quantile=brainlof$WMV, sGMV_quantile=brainlof$sGMV)


p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=GMV_quantile)) +
geom_violin()+geom_jitter() + theme_minimal() +  xlab("CNV Intolerance Score") + ylab("GMV centile")
  p
ggsave(paste(figpath,'suppGMV.pdf',sep = '/'))

  p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=WMV_quantile)) +
geom_violin()+geom_jitter() + theme_minimal() +  xlab("CNV Intolerance Score") + ylab("WMV centile")
  p
ggsave(paste(figpath,'suppWMV.pdf',sep = '/'))


p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=sGMV_quantile)) +
geom_violin()+geom_jitter() + theme_minimal() +  xlab("CNV Intolerance Score") + ylab("sGMV centile")
  p
ggsave(paste(figpath,'suppSubGMV.pdf',sep = '/'))

##redo with LOEUF


lof.cat <- rep('low risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$LOF>(1/.34)] <- 'high risk'

#lof.cat <- rep('0', length(brainlof$LOF))
#lof.cat[brainlof$LOF>0] <- '0-3'
#lof.cat[brainlof$LOF>(1/.34)] <- '3+'

toplot.raw <- table(lof.cat,deviat.all)
chisq.test(toplot.raw)

toplot <- table(lof.cat,deviat.all)
toplot[1,] <- toplot[1,]/sum(toplot[1,])
toplot[2,] <- toplot[2,]/sum(toplot[2,])

toplotdf <- data.frame(proportion=c(toplot[,1],toplot[,2]), CNV_Intolerance_Score=c(rownames(toplot), rownames(toplot)), 
  growthcart=c('low deviance', 'low deviance','high deviance','high deviance'))

p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=proportion, fill=growthcart)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

lab <- c('','n=912','n=38','')
toplotdf$lab <- lab
# Use custom colors
p + scale_fill_manual(values=c("red","#3B9AB2"))

p+geom_text(data=toplotdf,aes(x=CNV_Intolerance_Score,y=proportion,label=lab),vjust=-1)+ scale_fill_manual(values=c("red","#3B9AB2")) +
ggtitle("CNV risk and brain deviance \n (LOEUF annotation)")  + theme(plot.title = element_text(hjust = 0.5)) +  xlab("CNV Intolerance Score")

ggsave(paste(figpath,'supp_braindev_LOEUF.pdf',sep = '/'))


##redo with LOF cont

lof.cat <- rep('C) low risk', length(brainlof$pHI))
#lof.cat[brainlof$pHI>0 | brainlof$pTS>0 ] <- '0-1'
lof.cat[brainlof$LOF>0] <- 'B) medium risk'
lof.cat[brainlof$LOF>(1/.34)] <- 'A) high risk'

toplot.raw <- table(lof.cat,deviat.all)
chisq.test(toplot.raw)
toplot <- table(lof.cat,deviat.all)
toplot[1,] <- toplot[1,]/sum(toplot[1,])
toplot[2,] <- toplot[2,]/sum(toplot[2,])
toplot[3,] <- toplot[3,]/sum(toplot[3,])

toplotdf <- data.frame(proportion=c(toplot[,1],toplot[,2]), CNV_Intolerance_Score=c(rownames(toplot), 
  rownames(toplot)), growthcart=c('low deviance', 'low deviance','low deviance', 
  'high deviance','high deviance', 'high deviance'))



p <- ggplot(data=toplotdf, aes(x=CNV_Intolerance_Score, y=proportion, fill=growthcart)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

lab <- c('','n=210','n=702','n=38','', '')
toplotdf$lab <- lab
# Use custom colors

p+geom_text(data=toplotdf,aes(x=CNV_Intolerance_Score,y=proportion,label=lab),vjust=-1)+ scale_fill_manual(values=c("red","#3B9AB2")) +
ggtitle("CNV risk and brain deviance \n (LOEUF annotation with medium risk)")  + theme(plot.title = element_text(hjust = 0.5)) +  xlab("CNV Intolerance Score")
ggsave(paste(figpath,'supp_braindev_LOEUF_v2.pdf',sep = '/'))

