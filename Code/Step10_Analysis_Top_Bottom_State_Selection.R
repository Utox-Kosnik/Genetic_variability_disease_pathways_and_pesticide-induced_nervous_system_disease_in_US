############################################################################

### Identify top 10 and bottom 10 states
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())

library(data.table)
library(tidyverse)
library (plyr)
library(readxl)
library(grid)
library(gridExtra)
library(MASS)
library(scales)
library(ggrepel)
library(usa)

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)

setwd("")
load("Workspaces/Analysis_Workspace.RData")
#Recall: gene match describes the actual gene the SNP is in - genesymbol/ID describe the c-g-p-d link

############################################################################
#Function to prep data for a given year
year.select<- function(chem.dat, pest.dat, year.in){
  
  #Select pesticides applied in or before the year of study
  pest.list<- pest.dat[YEAR <= year.in, list(YEAR, COMPOUND, KG_APPLIED, State, CASRN)]
  pest.list<- pest.list[!duplicated(pest.list)]
  
  in.full.dat<- chem.dat[, .(.N), by=.(CASRN, ChemicalName, UMLS_Name, GeneID, Gene_Match, snpId)]
  usgs.year<- in.full.dat[CASRN %in% unique(pest.list$CASRN)]
  usgs.year<- merge(usgs.year, pest.list, by="CASRN", allow.cartesian=TRUE)
  usgs.year #Not combining years in this preparation
}

############################################################################
###Identify top 10 states (closest to the top right corner)
#and bottom 10 states (closest to bottom left corner)
#in plot of mean SNPs/pesticide vs disease occurrence

#Prepare input data
cpd.dat<- year.select(chem.snp.dis.dat, usgs.year.pest, 2018)
cpd.dat<- cpd.dat[, .(Var=length(unique(na.omit(snpId)))), by=.(State, UMLS_Name, CASRN, YEAR)] #SNPs per chem-year
cpd.dat<- cpd.dat[, .(Var=mean(Var)), by=.(location_name=State, Disease=UMLS_Name)] #SNP hits/pesticide-year

inc.mean.2018<- merge(cpd.dat, gbd.dat[measure_name %in% "Incidence" & year == 2018], allow.cartesian=TRUE)
prev.mean.2018<- merge(cpd.dat, gbd.dat[measure_name %in% "Prevalence" & year == 2018], allow.cartesian=TRUE)
mean.2018<- merge(cpd.dat, gbd.dat[year == 2018], allow.cartesian=TRUE)

############################################################################
#Function to plot top/bottom states

state.match<- as.data.table(cbind(state.name, state.abb))

trend.plot<- function(in.plot.dat){
  
  #Label states using abbreviations
  in.plot.dat$StateLabel<- in.plot.dat$location_name
  in.plot.dat[Label == 0]$StateLabel<- ""
  in.plot.dat$StateLabel<- state.match[match(in.plot.dat$StateLabel, state.match$state.name)]$state.abb
  
  in.plot.dat$Label<- as.factor(in.plot.dat$Label)
  
  #Alzheimer's disease
  us.alz<- in.plot.dat[Disease %in% "Alzheimer's Disease"]
  
  plot.alz.inc<- ggplot(data = us.alz, aes(x = Var, y = Incidence)) + geom_point(aes(fill=Label), pch=21) + theme_bw() + 
    ylab(paste("Disease incidence rate")) + xlab("SNP Hits/Pesticide") + geom_text_repel(aes(label=StateLabel), size=1.5) +
    facet_wrap(.~Disease, scales="free") + 
    geom_smooth(method = "rlm", se=FALSE) + 
    scale_fill_manual(values = c("2"="red", "1"="blue", "0"="grey"), guide=NULL)+
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  plot.alz.prev<- ggplot(data = us.alz, aes(x = Var, y = Prevalence)) + geom_point(aes(fill=Label), pch=21) + theme_bw() + 
    ylab(paste("Disease pravelence rate")) + xlab("SNP Hits/Pesticide") + geom_text_repel(aes(label=StateLabel), size=1.5) +
    facet_wrap(.~Disease, scales="free") + 
    geom_smooth(method = "rlm", se=FALSE) + 
    scale_fill_manual(values = c("2"="red", "1"="blue", "0"="grey"), guide=NULL)+
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  #Multiple sclerosis
  us.ms<- in.plot.dat[Disease %in% "Multiple Sclerosis"]
  
  plot.ms.inc<- ggplot(data = us.ms, aes(x = Var, y = Incidence)) + geom_point(aes(fill=Label), pch=21) + theme_bw() + 
    ylab(paste("Disease incidence rate")) + xlab("SNP Hits/Pesticide") + geom_text_repel(aes(label=StateLabel), size=1.5) +
    facet_wrap(.~Disease, scales="free") + 
    geom_smooth(method = "rlm", se=FALSE) + 
    scale_fill_manual(values = c("2"="red", "1"="blue", "0"="grey"), guide=NULL)+
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  plot.ms.prev<- ggplot(data = us.ms, aes(x = Var, y = Prevalence)) + geom_point(aes(fill=Label), pch=21) + theme_bw() + 
    ylab(paste("Disease pravelence rate")) + xlab("SNP Hits/Pesticide") + geom_text_repel(aes(label=StateLabel), size=1.5) +
    facet_wrap(.~Disease, scales="free") + 
    geom_smooth(method = "rlm", se=FALSE) + 
    scale_fill_manual(values = c("2"="red", "1"="blue", "0"="grey"), guide=NULL)+
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  #Parkinson disease
  us.pd<- in.plot.dat[Disease %in% "Parkinson Disease"]
  
  plot.pd.inc<- ggplot(data = us.pd, aes(x = Var, y = Incidence)) + geom_point(aes(fill=Label), pch=21) + theme_bw() + 
    ylab(paste("Disease incidence rate")) + xlab("SNP Hits/Pesticide") +  geom_text_repel(aes(label=StateLabel), size=1.5) +
    facet_wrap(.~Disease, scales="free") + 
    geom_smooth(method = "rlm", se=FALSE) + 
    scale_fill_manual(values = c("2"="red", "1"="blue", "0"="grey"), guide=NULL)+
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  plot.pd.prev<- ggplot(data = us.pd, aes(x = Var, y = Prevalence)) + geom_point(aes(fill=Label), pch=21) + theme_bw() + 
    ylab(paste("Disease pravelence rate")) +  xlab("SNP Hits/Pesticide") + geom_text_repel(aes(label=StateLabel), size=1.5) +
    facet_wrap(.~Disease, scales="free") + 
    geom_smooth(method = "rlm", se=FALSE) + 
    scale_fill_manual(values = c("2"="red", "1"="blue", "0"="grey"), guide=NULL)+
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  grid.arrange(plot.alz.inc, plot.ms.inc, plot.pd.inc, plot.alz.prev, plot.ms.prev, plot.pd.prev, nrow=2)
}

########################################################################################
###Distinguish top-bottom states, plot figures w/ labels/trendlines

#Prepare input data
mean.snp.dat<- mean.2018[Disease %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]
mean.snp.dat<- dcast(mean.snp.dat, location_name+Disease+Var~measure_name, value.var="val")

dis.list<- unique(mean.snp.dat$Disease)
top.bot.dat<- NULL

#Prepare separate list per disease
for(dis.in in dis.list){
  mean.snp.dat.in<- mean.snp.dat[Disease %in% dis.in]
  
  #Rescale values
  mean.snp.dat.in$Var<- rescale(mean.snp.dat.in$Var, to=c(0.1, 1))
  mean.snp.dat.in$Incidence<- rescale(mean.snp.dat.in$Incidence, to=c(0.1, 1))
  mean.snp.dat.in$Prevalence<- rescale(mean.snp.dat.in$Prevalence, to=c(0.1, 1))
  
  #Combine incidence and prevalence
  mean.snp.dat.in$Comb<- mean.snp.dat.in$Prevalence + mean.snp.dat.in$Incidence
  
  #Multiply values together
  mean.snp.dat.in$RankVal<- mean.snp.dat.in$Var * mean.snp.dat.in$Comb
  
  top.bot.dat<- rbind(top.bot.dat, mean.snp.dat.in)
}

#Order by top/bottom ranks per disease
top.bot.dat<- top.bot.dat[order(Disease, -RankVal)]
top.bot.dat$Label<- 0

#Identify top ranks
plot.top<- rbind(head(top.bot.dat[Disease %in% "Alzheimer's Disease"], 10),
                 head(top.bot.dat[Disease %in% "Multiple Sclerosis"], 10),
                 head(top.bot.dat[Disease %in% "Parkinson Disease"], 10))
plot.top$Label<- 2

#Identify bottom ranks
plot.bot<- rbind(tail(top.bot.dat[Disease %in% "Alzheimer's Disease"], 10),
                 tail(top.bot.dat[Disease %in% "Multiple Sclerosis"], 10),
                 tail(top.bot.dat[Disease %in% "Parkinson Disease"], 10))
plot.bot$Label<- 1

#Create overall dataset
top.bot.dat<- rbind(top.bot.dat, plot.top, plot.bot)
top.bot.dat<- top.bot.dat[ , .(Label=max(Label)), by=.(location_name, Disease, RankVal, Var, Incidence, Prevalence)]

top.bot.dat[Disease %in% "Alzheimer's Disease" & Label != 0]
top.bot.dat[Disease %in% "Multiple Sclerosis" & Label != 0]
top.bot.dat[Disease %in% "Parkinson Disease" & Label != 0]

###NOTE: SAVE DATA FOR FUTURE ANALYSIS
#Copy of dataset available currently as "Top_Bottom_States_MeanSNPs.csv" used in subsequent steps
setwd("")
fwrite(top.bot.dat, file="Generated_Data/FILENAME.csv")

###Prepare data for plot
snp.plot<- merge(mean.snp.dat, top.bot.dat[, list(location_name, Disease, Label, RankVal)], by=c("location_name", "Disease"))
top.bot.plot<- trend.plot(snp.plot)

setwd("")
jpeg("SITables_Figures/FigureS21.jpeg", width=6, height=3.5, units="in", res=300)
grid.draw(top.bot.plot)
graphics.off()

setwd("")
fwrite(snp.plot, "SITables_Figures/SI_FigureS21.csv")

