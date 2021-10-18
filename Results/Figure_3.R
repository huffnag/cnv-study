
# Load packages -----
library(tidyverse) # for manipulating data 
library(readxl) # for reading microsoft excel files
library(sjPlot) # for regression model tables
library(gridExtra)
library(broom)
path <- '/Users/aa2227/Documents/pncCNV/clean'
figpath <- paste(path,'figs',sep='/')


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

vars <- c('scale(envSES)','scale(PRS.IQ)','scale(pHI)','scale(Trauma)')

list <- list %>% 
  mutate(fdr=p.adjust(p.value,method = 'fdr')) %>% 
  mutate(signif=ifelse(fdr<0.05,'yes','no')) %>% 
  filter(term %in% vars)
  

listbind <- list %>% 
  mutate(va='Models of cognition')


listbind$response <- factor(listbind$response,levels = rev(c('Overall Accuracy','Executive Complex Cognition Accuracy','Memory Accuracy',
                                                             "Social Cognition Accuracy","Overall Speed","Memory Speed",
                                                             "Slow Speed","Fast Speed")))

iq <- ggplot(data=list,aes(x=abs(estimate), y=response,alpha=factor(signif),color=term)) +
  scale_alpha_discrete(range=c(0.35, 1),name = "pFDR < 0.05", labels = c("False", "True")) +
  geom_point() +
  geom_vline(xintercept = 0,lty=2) +
  xlab("Standardized Effect") +
  ylab("") +
  scale_color_discrete(name = "Term")+
  theme_linedraw()#+
#xlim(0,.3)

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

listbind2 <- list %>% 
  mutate(va='Models of psychopathology')

listbind2$response <- factor(listbind2$response,levels = rev(c("Bifactor Psychosis","Bifactor Overall Psychopathology",
                                                               "Bifactor Mood","Bifactor Fear","Bifactor Externalizing")))

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

tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)}

# figure 3
fig <- ggplot(data=full,aes(x=estimate, y=response,alpha=factor(signif),color=term)) +
  scale_alpha_discrete(range=c(0.25, 1),name = "pFDR < 0.05", labels = c("False", "True")) +
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high),position = position_jitter(height = 0.2)) +
  geom_vline(xintercept = 0,lty=2) +
  xlab("Standardized Effect") +
  ylab("") +
  scale_color_discrete(name = "Term", labels = c("SES", "pHI del.",'PGS-g','Trauma'))+
  theme_light()+
  facet_wrap(~va,ncol = 1,shrink = T,drop = T,scales = 'free')+
  theme(text = element_text(size=12))+
  theme(strip.text = element_text(colour = 'white'))+
  coord_cartesian(xlim = c(-.15,.4),clip = 'off')+
  theme(strip.background =element_rect(fill="black"))

pdf(file = paste(figpath,'fig3.pdf',sep = '/'),width = 10,height = 8)
tag_facet2(fig,tag_pool = LETTERS,open = '',x=-.25,vjust = .4)
dev.off()


  
