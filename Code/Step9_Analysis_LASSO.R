############################################################################

### LASSO
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library(plyr)
library(MASS)
library(readxl)
library(grid)
library(gridExtra)
library(glmnet)
library(caret)
library(c060)

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)

setwd("")
load("Workspaces/Analysis_Workspace.RData")
#Recall: gene match describes the actual gene the SNP is in - genesymbol/ID describe the c-g-p-d link


############################################################################
#Function to prep data for a given year
year.select<- function(chem.dat, pest.dat, year.in){
  
  #Select pesticides applied in or before the year of study
  pest.list<- pest.dat[YEAR <= year.in, list(YEAR, COMPOUND, KG_APPLIED, State, CASRN)]
  pest.list<- pest.list[, .(YearsApplied=length(unique(YEAR)), KG_APPLIED=sum(KG_APPLIED)), 
                        by=.(COMPOUND, State, CASRN)]
  pest.list<- pest.list[!duplicated(pest.list)]
  
  in.full.dat<- chem.dat[, .(.N), by=.(CASRN, ChemicalName, UMLS_Name, GeneID, Gene_Match, snpId)]
  usgs.year<- in.full.dat[CASRN %in% unique(pest.list$CASRN)]
  
  #Merge chem-snp-disease links with pesticide data per year
  usgs.year<- merge(usgs.year, pest.list, by="CASRN", allow.cartesian=TRUE)
  usgs.year #Now have Pesticide - SNP - Disease links per state per year
}

############################################################################
#Prepare data for LASSO regression
#Will do sum of SNP hits across years per chemical for each state

#Set up x dataset
#Overall dataset - will split across diseases
x.dat<- year.select(chem.snp.dis.dat, usgs.year.pest, year.in=2018)

#Prep individual datasets with chem values represented by # SNP hits
x.prep<- function(dis.in){
  x.in<- x.dat[UMLS_Name %in% dis.in, .(SNPs=length(unique(snpId))), 
               by=.(CASRN, YearsApplied, ChemicalName, State)] #determine number of SNPs per chemical
  x.in$FillVal<- x.in$YearsApplied * x.in$SNPs #Determine SNP hits based on years applied
  x.in<- dcast(x.in, State ~ ChemicalName, value.var="FillVal", fill=0)
  x.in[order(State)]
}

chems.per.dis<- x.dat[, .(length(unique(CASRN))), by=.(UMLS_Name)]

x.alz<- x.prep("Alzheimer's Disease")
x.ms<- x.prep("Multiple Sclerosis")
x.pd<- x.prep("Parkinson Disease")

#Set up y dataset - disease incidence/prevalence per state in 2018
y.dat<- dcast(gbd.dat[year == 2018], Disease + location_name ~ measure_name, value.var="val")
y.dat<- y.dat[order(Disease, location_name)]

############################################################################
#Incidence LASSO

###Alzheimers incidence
lasso.inc.alz<- stabpath(y=unlist(as.vector(y.dat[Disease %in% "Alzheimer's Disease", Incidence])), 
                         x=x.alz[, 2:ncol(x.alz)], size=.66, steps=10000)

#Select a stable set of variables based on per-comparison error rate
alz.inc.stab.val<- stabsel(lasso.inc.alz, error=0.05, type="pcer")$lpos
alz.inc.stab<- as.data.table(lasso.inc.alz$x[, alz.inc.stab.val], keep.rownames = TRUE)
length(unique(alz.inc.stab[V2 != 0]$V1))

###MS incidence
lasso.inc.ms<- stabpath(y=unlist(as.vector(y.dat[Disease %in% "Multiple Sclerosis", Incidence])), 
                        x=x.ms[, 2:ncol(x.ms)], size=.66, steps=10000)

#Select a stable set of variables based on per-comparison error rate
ms.inc.stab.val<- stabsel(lasso.inc.ms, error=0.05, type="pcer")$lpos
ms.inc.stab<- as.data.table(lasso.inc.ms$x[, ms.inc.stab.val], keep.rownames = TRUE)
length(unique(ms.inc.stab[V2 != 0]$V1))

###PD incidence
lasso.inc.pd<- stabpath(y=unlist(as.vector(y.dat[Disease %in% "Parkinson Disease", Incidence])), 
                        x=x.pd[, 2:ncol(x.pd)], size=.66, steps=10000)

#Select a stable set of variables based on per-comparison error rate
pd.inc.stab.val<- stabsel(lasso.inc.pd, error=0.05, type="pcer")$lpos
pd.inc.stab<- as.data.table(lasso.inc.pd$x[, pd.inc.stab.val], keep.rownames = TRUE)
length(unique(pd.inc.stab[V2 != 0]$V1))

############################################################################
#Prevalence LASSO

###Alzheimers prevalence
lasso.prev.alz<- stabpath(y=unlist(as.vector(y.dat[Disease %in% "Alzheimer's Disease", Prevalence])), 
                          x=x.alz[, 2:ncol(x.alz)], size=.66, steps=10000)

#Select a stable set of variables based on per-comparison error rate
alz.prev.stab.val<- stabsel(lasso.prev.alz, error=0.05, type="pcer")$lpos
alz.prev.stab<- as.data.table(lasso.prev.alz$x[, alz.prev.stab.val], keep.rownames = TRUE)
length(unique(alz.prev.stab[V2 != 0]$V1))

###MS prevalence
lasso.prev.ms<- stabpath(y=unlist(as.vector(y.dat[Disease %in% "Multiple Sclerosis", 4])), 
                         x=x.ms[, 2:ncol(x.ms)], size=.66, steps=10000)

#Select a stable set of variables based on per-comparison error rate
ms.prev.stab.val<- stabsel(lasso.prev.ms, error=0.05, type="pcer")$lpos
ms.prev.stab<- as.data.table(lasso.prev.ms$x[, ms.prev.stab.val], keep.rownames = TRUE)
length(unique(ms.prev.stab[V2 != 0]$V1))

###PD prevalence
lasso.prev.pd<- stabpath(y=unlist(as.vector(y.dat[Disease %in% "Parkinson Disease", 4])), 
                         x=x.pd[, 2:ncol(x.pd)], size=.66, steps=10000)

#Select a stable set of variables based on per-comparison error rate
pd.prev.stab.val<- stabsel(lasso.prev.pd, error=0.05, type="pcer")$lpos
pd.prev.stab<- as.data.table(lasso.prev.pd$x[, pd.prev.stab.val], keep.rownames = TRUE)
length(unique(pd.prev.stab[V2 != 0]$V1))

############################################################################
###Analyze results

#Which pesticides received 0 coefficients
alz.inc.stab[V2 == 0] 
alz.prev.stab[V2 == 0]
x.alz[, which(colnames(x.alz) %in% rownames(lasso.prev.alz$x[which(rowSums(lasso.prev.alz$x)==0),])), with=FALSE]

ms.inc.stab[V2 == 0]
ms.prev.stab[V2 == 0]
x.ms[, which(colnames(x.ms) %in% rownames(lasso.prev.ms$x[which(rowSums(lasso.prev.ms$x)==0),])), with=FALSE]

pd.inc.stab[V2 == 0]
pd.prev.stab[V2 == 0]
x.pd[, which(colnames(x.pd) %in% rownames(lasso.prev.pd$x[which(rowSums(lasso.prev.pd$x)==0),])), with=FALSE]

#See frequently used pesticides with the same application across states or application in just one state

########################################################################################
###Now assess whether the same accuracy is seen with the rlm with the stable set of coefficients
#For each combination of features, build a rlm and compare coefficients
#Prepare functions/datasets

#Set colors based on geographic regions of the US
usgs.year.pest$Region<- usgs.year.pest$State
usgs.year.pest[Region %in% c("Washington", "Oregon", "Idaho", "Montana", "Wyoming", 
                             "California", "Nevada", "Utah", "Colorado")]$Region<- "West"
usgs.year.pest[Region %in% c("Arizona", "New Mexico", "Texas", "Oklahoma")]$Region<- "Southwest"
usgs.year.pest[Region %in% c("North Dakota", "Minnesota", "Wisconsin", "Michigan", "South Dakota", 
                             "Nebraska", "Iowa", "Illinois", "Indiana", "Ohio", "Kansas", "Missouri")]$Region<- "Midwest"
usgs.year.pest[Region %in% c("West Virginia", "Virginia", "Kentucky", "Tennessee", "North Carolina", "Arkansas",
                             "Mississippi", "Alabama", "Georgia", "South Carolina", "Louisiana", "Florida")]$Region<- "Southeast"
usgs.year.pest[State == Region]$Region<- "Northeast"

#Function to generate correlation plot

#Same function as in Step 7
#Shape and color based on geographic region
dis.cor.plot<- function(plot.in, year.in, measure.type, cpd.lab){
  
  plot.out<- ggplot(data = plot.in, aes(x = Var, y = val)) + geom_point(aes(color=Region, shape=Region)) + 
    theme_bw() + 
    ylab(paste("Disease", measure.type, "Rate")) + xlab(cpd.lab) + 
    scale_color_manual(values=c("West"="tomato","Southwest"="yellow", "Midwest"="springgreen3", 
                                "Southeast"="skyblue2", 
                                "Northeast"="mediumorchid2")) + 
    scale_x_continuous(expand=c(.1,.1)) + scale_y_continuous(expand=c(.1,.1)) +
    facet_wrap(.~Disease, scales="free", nrow=1) + geom_smooth(method = "rlm") + #scale_y_continuous(trans="log10") +
    ggtitle(paste(cpd.lab, "vs Disease", measure.type, "Rate")) + 
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))
  
  return(plot.out)
}  

########################################################################################
#Function to prepare stable dataset
lasso.plot.dat<- function(lasso.in.dat, usgs.in.dat, pest.in.dat, gbd.dat.in, measure.type){
  
  #Only include chemicals with LASSO coefficients - excludes many common pesticides
  lasso.in.dat<- lasso.in.dat[V2 != 0]
  usgs.in.dat<- usgs.in.dat[ChemicalName %in% lasso.in.dat$V1]
  usgs.in.dat<- usgs.in.dat[, .(Var=length(unique(na.omit(snpId)))), by=.(State, UMLS_Name, CASRN)]
  pest.in.dat<- pest.in.dat[YEAR <= 2018, list(YEAR, State, CASRN, Region)]
  usgs.in.dat<- merge(usgs.in.dat, pest.in.dat, by=c("State", "CASRN"), allow.cartesian=TRUE)
  
  gbd.dat.plot<- gbd.dat.in[measure_name %in% measure.type & year == 2018]
  
  dat.out<- merge(usgs.in.dat[, .(Var=mean(Var)), by=.(location_name=State, Disease=UMLS_Name, Region)], 
                  gbd.dat.plot, by=c("location_name", "Disease"), allow.cartesian=TRUE)
  
  dat.out
  
}

usgs.map<- year.select(chem.snp.dis.dat, usgs.year.pest, 2018)

########################################################################################
###Now assess whether the same accuracy is seen with the rlm with the stable set of coefficients
#For each combination of features, build a rlm and compare coefficients

rlm.out<- matrix(data=NA, nrow=1, ncol=12, 
                 dimnames = list(NULL, c("Alz_RLM_IncMAE", "Alz_RLM_IncRMSE", "Alz_RLM_PrevMAE", "Alz_RLM_PrevRMSE",
                                         "MS_RLM_IncMAE", "MS_RLM_IncRMSE", "MS_RLM_PrevMAE", "MS_RLM_PrevRMSE",
                                         "PD_RLM_IncMAE", "PD_RLM_IncRMSE", "PD_RLM_PrevMAE", "PD_RLM_PrevRMSE")))
pear.out<- matrix(data=NA, nrow=1, ncol=12, 
                   dimnames = list(NULL, c("Alz_IncRho", "Alz_IncpVal", "Alz_PrevRho", "Alz_PrevpVal",
                                           "MS_IncRho", "MS_IncpVal", "MS_PrevRho", "MS_PrevpVal",
                                           "PD_IncRho", "PD_IncpVal", "PD_PrevRho", "PD_PrevpVal")))

#Prepare stable dataset per disease/occurrence type
alz.inc.plot<- lasso.plot.dat(alz.inc.stab, usgs.map[UMLS_Name %in% "Alzheimer's Disease"], usgs.year.pest,
                              gbd.dat[Disease %in% "Alzheimer's Disease"], measure.type="Incidence")
alz.prev.plot<- lasso.plot.dat(alz.prev.stab, usgs.map[UMLS_Name %in% "Alzheimer's Disease"], usgs.year.pest,
                               gbd.dat[Disease %in% "Alzheimer's Disease"], measure.type="Prevalence")

ms.inc.plot<- lasso.plot.dat(ms.inc.stab, usgs.map[UMLS_Name %in% "Multiple Sclerosis"], usgs.year.pest,
                             gbd.dat[Disease %in% "Multiple Sclerosis"], measure.type="Incidence")
ms.prev.plot<- lasso.plot.dat(ms.prev.stab, usgs.map[UMLS_Name %in% "Multiple Sclerosis"], usgs.year.pest,
                              gbd.dat[Disease %in% "Multiple Sclerosis"], measure.type="Prevalence")

pd.inc.plot<- lasso.plot.dat(pd.inc.stab, usgs.map[UMLS_Name %in% "Parkinson Disease"], usgs.year.pest,
                             gbd.dat[Disease %in% "Parkinson Disease"], measure.type="Incidence")
pd.prev.plot<- lasso.plot.dat(pd.prev.stab, usgs.map[UMLS_Name %in% "Parkinson Disease"], usgs.year.pest,
                              gbd.dat[Disease %in% "Parkinson Disease"], measure.type="Prevalence")

#Generate plots
inc.mean.2018<- dis.cor.plot(rbind(alz.inc.plot, ms.inc.plot, pd.inc.plot), 2018, measure.type="Incidence", "SNP Hits/Pesticide")
prev.mean.2018<- dis.cor.plot(rbind(alz.prev.plot, ms.prev.plot, pd.prev.plot), 2018, measure.type="Prevalence", "SNP Hits/Pesticide")

all.plot.out<- arrangeGrob(inc.mean.2018, prev.mean.2018, nrow=2)

#Determine RLM stats
rlm.inc.alz <- rlm(val~Var, data=alz.inc.plot)
rlm.prev.alz <- rlm(val~Var, data=alz.prev.plot)
rlm.inc.ms <- rlm(val~Var, data=ms.inc.plot)
rlm.prev.ms <- rlm(val~Var, data=ms.prev.plot)
rlm.inc.pd <- rlm(val~Var, data=pd.inc.plot)
rlm.prev.pd <- rlm(val~Var, data=pd.prev.plot)

rlm.out[1,]<- c(MAE(obs=alz.inc.plot$val, pred=predict(rlm.inc.alz)),
                RMSE(obs=alz.inc.plot$val, pred=predict(rlm.inc.alz)),
                MAE(obs=alz.prev.plot$val, pred=predict(rlm.prev.alz)),
                RMSE(obs=alz.prev.plot$val, pred=predict(rlm.prev.alz)),
                
                MAE(obs=ms.inc.plot$val, pred=predict(rlm.inc.ms)),
                RMSE(obs=ms.inc.plot$val, pred=predict(rlm.inc.ms)),
                MAE(obs=ms.prev.plot$val, pred=predict(rlm.prev.ms)),
                RMSE(obs=ms.prev.plot$val, pred=predict(rlm.prev.ms)),
                
                MAE(obs=pd.inc.plot$val, pred=predict(rlm.inc.pd)),
                RMSE(obs=pd.inc.plot$val, pred=predict(rlm.inc.pd)),
                MAE(obs=pd.prev.plot$val, pred=predict(rlm.prev.pd)),
                RMSE(obs=pd.prev.plot$val, pred=predict(rlm.prev.pd)))

#Determine Pearson's correlation stats
pear.inc.alz <- cor.test(alz.inc.plot$val, alz.inc.plot$Var, method="pearson")
pear.prev.alz <- cor.test(alz.prev.plot$val, alz.prev.plot$Var, method="pearson")
pear.inc.ms <- cor.test(ms.inc.plot$val, ms.inc.plot$Var, method="pearson")
pear.prev.ms <- cor.test(ms.prev.plot$val, ms.prev.plot$Var, method="pearson")
pear.inc.pd <- cor.test(pd.inc.plot$val, pd.inc.plot$Var, method="pearson")
pear.prev.pd <- cor.test(pd.prev.plot$val, pd.prev.plot$Var, method="pearson")

pear.out[1,]<- c(pear.inc.alz$estimate, pear.inc.alz$p.value, pear.prev.alz$estimate, pear.prev.alz$p.value,
                  pear.inc.ms$estimate, pear.inc.ms$p.value, pear.prev.ms$estimate, pear.prev.ms$p.value,
                  pear.inc.pd$estimate, pear.inc.pd$p.value, pear.prev.pd$estimate, pear.prev.pd$p.value)


#Compare correlation between reduced coefficents and all (Fig 2)
setwd("")
jpeg("SITables_Figures/FigureS20.jpeg", width=5.5, height=3.5, units="in", res=300)
grid.draw(all.plot.out)
graphics.off()
###Lower compared to original

#Figure S20 underlying data
length(which(alz.inc.plot$location_name != alz.prev.plot$location_name))
length(which(ms.inc.plot$Disease != ms.prev.plot$Disease))
length(which(pd.inc.plot$Var != pd.prev.plot$Var))
figs20.dat<- rbind(alz.inc.plot[, list(State=location_name, Region, Disease, PestHits=Var, Metric=measure_name,
                                        DiseaseOccurrence=val)], 
                   alz.prev.plot[, list(State=location_name, Region, Disease, PestHits=Var, Metric=measure_name,
                                       DiseaseOccurrence=val)],
                   ms.inc.plot[, list(State=location_name, Region, Disease, PestHits=Var, Metric=measure_name,
                                       DiseaseOccurrence=val)], 
                   ms.prev.plot[, list(State=location_name, Region, Disease, PestHits=Var, Metric=measure_name,
                                        DiseaseOccurrence=val)],
                   pd.inc.plot[, list(State=location_name, Region, Disease, PestHits=Var, Metric=measure_name,
                                       DiseaseOccurrence=val)], 
                   pd.prev.plot[, list(State=location_name, Region, Disease, PestHits=Var, Metric=measure_name,
                                        DiseaseOccurrence=val)])


setwd("")
fwrite(figs20.dat, "SITables_Figures/SI_FigureS20.csv")

########################################################################################
#Identify top coefficients per disease

top.alz<- rbind(alz.inc.stab, alz.prev.stab)
top.alz<- top.alz[, .(Mean=mean(V2), Sum=sum(V2)), by=.(V1)]
top.alz.sum<- top.alz[order(-Sum)]
top.alz.mean<- top.alz[order(-Mean)]
length(which(top.alz.sum$V1 != top.alz.mean$V1))

top.ms<- rbind(ms.inc.stab, ms.prev.stab)
top.ms<- top.ms[, .(Mean=mean(V2), Sum=sum(V2)), by=.(V1)]
top.ms.sum<- top.ms[order(-Sum)]
top.ms.mean<- top.ms[order(-Mean)]
length(which(top.ms.sum$V1 != top.ms.mean$V1))

top.pd<- rbind(pd.inc.stab, pd.prev.stab)
top.pd<- top.pd[, .(Mean=mean(V2), Sum=sum(V2)), by=.(V1)]
top.pd.sum<- top.pd[order(-Sum)]
top.pd.mean<- top.pd[order(-Mean)]
length(which(top.pd.sum$V1 != top.pd.mean$V1))

top.alz.sum$Disease<- "Alzheimer's Disease"
top.ms.sum$Disease<- "Multiple Sclerosis"
top.pd.sum$Disease<- "Parkinson Disease"

top.chem.out<- rbind(top.alz.sum, top.ms.sum, top.pd.sum)
colnames(top.chem.out)[1]<- "ChemicalName"
top.chem.out$CASRN<- usgs.map[match(top.chem.out$ChemicalName, usgs.map$ChemicalName)]$CASRN

#Prepare SI table
top.chem.print<- top.chem.out[Sum != 0, list(Disease, CASRN, ChemicalName, Sum)]

setwd("")
fwrite(top.chem.print, "SITables_Figures/Dataset2-LASSO_Coefficients.csv")

########################################################################################
###Assess top coefficients

quantile(top.chem.out[Disease %in% "Alzheimer's Disease"]$Mean)
quantile(top.chem.out[Disease %in% "Multiple Sclerosis"]$Mean)
quantile(top.chem.out[Disease %in% "Parkinson Disease"]$Mean)

top.chem.out[Disease %in% "Alzheimer's Disease" & Mean > 0.5]
top.chem.out[Disease %in% "Multiple Sclerosis" & Mean > 0.5]
top.chem.out[Disease %in% "Parkinson Disease" & Mean > 0.5]

top.chem.out$Rank<- c(rank(1/top.chem.out[Disease %in% "Alzheimer's Disease"]$Sum),
                      rank(1/top.chem.out[Disease %in% "Multiple Sclerosis"]$Sum),
                      rank(1/top.chem.out[Disease %in% "Parkinson Disease"]$Sum))

#Test differences in pesticide ranks
test.alz<- cbind(top.alz.sum, 1:nrow(top.alz.sum))
test.ms<-  cbind(top.ms.sum, 1:nrow(top.ms.sum))

test.alz<- test.alz[V1 %in% test.ms$V1]
test.ms<- test.ms[match(test.alz$V1, test.ms$V1)]

wilcox.test(test.alz$V2, test.ms$V2) 
wilcox.test(test.alz$Sum, test.ms$Sum) 

wilcox.test(test.alz[1:100]$V2, test.ms[1:100]$V2) #Top 100 significantly different

#Test differences in pesticide ranks
test.alz<- cbind(top.alz.sum, 1:nrow(top.alz.sum))
test.pd<-  cbind(top.pd.sum, 1:nrow(top.pd.sum))

test.alz<- test.alz[V1 %in% test.pd$V1]
test.pd<- test.pd[match(test.alz$V1, test.pd$V1)]

wilcox.test(test.alz$V2, test.pd$V2) 
wilcox.test(test.alz$Sum, test.pd$Sum) 

wilcox.test(test.alz[1:100]$V2, test.pd[1:100]$V2) #Full set isn't significantly different

#Test differences in pesticide ranks
test.ms<- cbind(top.ms.sum, 1:nrow(top.ms.sum))
test.pd<-  cbind(top.pd.sum, 1:nrow(top.pd.sum))

test.ms<- test.ms[V1 %in% test.pd$V1]
test.pd<- test.pd[match(test.ms$V1, test.pd$V1)]

wilcox.test(test.ms$V2, test.pd$V2) 
wilcox.test(test.ms$Sum, test.pd$Sum) 

wilcox.test(test.ms[1:100]$V2, test.pd[1:100]$V2) #Full set isn't significantly different

