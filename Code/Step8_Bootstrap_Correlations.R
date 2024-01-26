############################################################################

### Bootstrapping of kg pest applied and SNPs/pesticide-year vs disease occurrence 
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

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)

setwd("")
load("Workspaces/Analysis_Workspace.RData")
#Recall: gene match describes the actual gene the SNP is in - genesymbol/ID describe the c-g-p-d link

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

############################################################################
#Function to prep data for a given year

#Add code to shuffle data
year.select<- function(chem.dat, pest.dat, year.in, shuffle.data=FALSE){
  
  #Select pesticides applied in or before the year of study
  pest.list<- pest.dat[YEAR <= year.in, list(YEAR, COMPOUND, KG_APPLIED, State, CASRN, Region)]
  pest.list<- pest.list[, .(KG_APPLIED=sum(KG_APPLIED)), 
                        by=.(COMPOUND, State, CASRN, YEAR, Region)]
  
  if(shuffle.data == TRUE){
    
    #For each state, keep the number of pesticides the same, but change the 
    #compound and kg applied
    pest.list$Shuffle<- pest.list[sample(nrow(pest.list))]$COMPOUND
    pest.list$CAS_pe<- pest.list[match(pest.list$Shuffle, pest.list$COMPOUND)]$CASRN
    pest.list$KG_APPLIED<- pest.list[match(pest.list$Shuffle, pest.list$COMPOUND)]$KG_APPLIED
    
    in.full.dat<- chem.dat[!is.na(snpId), .(.N), by=.(CASRN, ChemicalName, UMLS_Name, GeneID, Gene_Match, snpId)]
    usgs.year<- in.full.dat[CASRN %in% unique(pest.list$CASRN)]
    usgs.year<- merge(usgs.year, pest.list, by.x="CASRN", by.y="CAS_pe", allow.cartesian=TRUE)
    
  } else{
    
    in.full.dat<- chem.dat[, .(.N), by=.(CASRN, ChemicalName, UMLS_Name, GeneID, Gene_Match, snpId)]
    usgs.year<- in.full.dat[CASRN %in% unique(pest.list$CASRN)]
    usgs.year<- merge(usgs.year, pest.list, by="CASRN", allow.cartesian=TRUE)
    
  }
  usgs.year
}

############################################################################
###Function for correlation data prep

dis.cor.prep<- function(cpd.dat, pest.in.dat, gbd.dat.in, year.in, measure.type, var.type="SNP", method.type="Mean", shuffle.data=FALSE){
  
  if(shuffle.data == TRUE){
    cpd.dat<- year.select(cpd.dat, pest.in.dat, year.in, shuffle.data=TRUE)
  } else{
    cpd.dat<- year.select(cpd.dat, pest.in.dat, year.in)
  }
  
  if(var.type == "SNP"){
    cpd.dat<- cpd.dat[, .(Var=length(unique(na.omit(snpId)))), by=.(State, UMLS_Name, CASRN, Region, YEAR)]
  }
  
  if(var.type == "Gene"){
    cpd.dat<- cpd.dat[, .(Var=length(unique(GeneID))), by=.(State, UMLS_Name, CASRN, Region, YEAR)]
  }
  
  if(var.type == "KG"){

    cpd.dat<- cpd.dat[, .(.N), by=.(Var=KG_APPLIED, State, CASRN, UMLS_Name, Region, YEAR)]
    cpd.dat<- cpd.dat[, .(Var=sum(Var)), by=.(State, UMLS_Name, Region, YEAR)]
  }
  
  if(method.type=="Sum"){
    cpd.dat<- cpd.dat[, .(Var=sum(Var)), by=.(location_name=State, Disease=UMLS_Name, Region)]
  } 
  
  if(method.type=="Mean"){
    cpd.dat<- cpd.dat[, .(Var=mean(Var)), by=.(location_name=State, Disease=UMLS_Name, Region)]
  }
  
  gbd.dat.plot<- gbd.dat.in[measure_name %in% measure.type & year == year.in]
  
  plot.dat<- merge(cpd.dat, gbd.dat.plot, by=c("location_name", "Disease"), allow.cartesian=TRUE)
  plot.dat
}

############################################################################
###Function for correlation calculations

cor.calc<- function(in.dat.inc, in.dat.prev){
  
  pe.dat<- matrix(data=NA, nrow=length(unique(in.dat.inc$Disease)), ncol=4, 
                  dimnames = list(NULL, c("IncCoeff_pe", "PrevCoeff_pe", "Incp_pe", "Prevp_pe"
                                          )))
  dis.list<- unique(in.dat.inc$Disease)
  for(dis.in in dis.list){
    
    pe.inc.int <- cor.test(in.dat.inc[Disease == dis.in]$Var, 
                           in.dat.inc[Disease == dis.in]$val, method="pearson")
    pe.prev.int <- cor.test(in.dat.prev[Disease == dis.in]$Var, 
                            in.dat.prev[Disease == dis.in]$val, method="pearson")
    
    pe.dat[which(dis.list %in% dis.in),]<- c(pe.inc.int$estimate, pe.prev.int$estimate, 
                                             pe.inc.int$p.value, pe.prev.int$p.value
                                             )
  }
  
  pe.dat<- as.data.table(pe.dat)
  pe.dat[, 1:4]<- lapply(pe.dat[,1:4], function(x) as.numeric(as.character(x)))
  
  lm.dat<- matrix(data=NA, nrow=length(unique(in.dat.inc$Disease)), ncol=9, 
                  dimnames = list(NULL, c("RLM_IncSlope", "RLM_PrevSlope", "RLM_IncInt", "RLM_PrevInt", 
                                          "RLM_IncMAE", "RLM_PrevMAE", "RLM_IncRMSE", "RLM_PrevRMSE",
                                          "Disease")))
  
  for(dis.in in dis.list){
    
    rlm.inc.int <- rlm(val~Var, data=in.dat.inc[Disease %in% dis.in])
    rlm.prev.int <- rlm(val~Var, data=in.dat.prev[Disease %in% dis.in])
    
    lm.dat[which(dis.list %in% dis.in),]<- c(rlm.inc.int$coefficients[2], rlm.prev.int$coefficients[2], 
                                             rlm.inc.int$coefficients[1], rlm.prev.int$coefficients[1], 
                                             MAE(obs=in.dat.inc[Disease %in% dis.in]$val, pred=predict(rlm.inc.int)),
                                             MAE(obs=in.dat.prev[Disease %in% dis.in]$val, pred=predict(rlm.prev.int)),
                                             RMSE(obs=in.dat.inc[Disease %in% dis.in]$val, pred=predict(rlm.inc.int)),
                                             RMSE(obs=in.dat.prev[Disease %in% dis.in]$val, pred=predict(rlm.prev.int)),
                                             
                                             dis.in)
  }
  
  lm.dat<- as.data.table(lm.dat)
  lm.dat[, 1:8]<- lapply(lm.dat[,1:8], function(x) as.numeric(as.character(x)))
  
  cbind(pe.dat, lm.dat)
}


############################################################################
###Shuffle pesticide-SNP data and assess change in correlation strength

inc.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence")
prev.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence")

orig.stats<- cor.calc(inc.mean.2018, prev.mean.2018)

#Try shuffling
inc.shuf.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence", shuffle.data=TRUE)
prev.shuf.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence", shuffle.data=TRUE)

shuf.stats<- cor.calc(inc.shuf.2018, prev.shuf.2018)

############################################################################
###Shuffle kg applied data and assess change in correlation strength

inc.kg.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence", var.type="KG", method.type="Sum")
prev.kg.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence", var.type="KG", method.type="Sum")

orig.kg.stats<- cor.calc(inc.kg.2018, prev.kg.2018)

inc.shuf.kg.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence", shuffle.data=TRUE, var.type="KG", method.type="Sum")
prev.shuf.kg.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence", shuffle.data=TRUE, var.type="KG", method.type="Sum")

shuf.kg.stats<- cor.calc(inc.shuf.kg.2018, prev.shuf.kg.2018)

############################################################################
###Try running 1000 times w/ pesticide shuffling, collect slopes, RMSE, MAE, and pearson's correlations for each run
#Analyze change in relationship for both SNPs/pesticide and KG pesticide

shuffle.dat.allpest.snp<- NULL
shuffle.dat.allpest.kg<- NULL

for(i in c(1:1000)){
  
  #Shuffle SNPs/pesticide
  inc.shuf.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence", shuffle.data=TRUE)
  prev.shuf.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence", shuffle.data=TRUE)
  shuf.stats<- cor.calc(inc.shuf.2018, prev.shuf.2018)
  
  shuf.stats$Rep<- i
  shuffle.dat.allpest.snp<- rbind(shuffle.dat.allpest.snp, shuf.stats)
  
  #Shuffle kg pesticides applied
  inc.shuf.kg.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence", shuffle.data=TRUE, var.type="KG", method.type="Sum")
  prev.shuf.kg.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence", shuffle.data=TRUE, var.type="KG", method.type="Sum")
  shuf.kg.stats<- cor.calc(inc.shuf.kg.2018, prev.shuf.kg.2018)
  
  shuf.kg.stats$Rep<- i
  shuffle.dat.allpest.kg<- rbind(shuffle.dat.allpest.kg, shuf.kg.stats)
  
  if(i %in% seq(1, 1000, 10)){
    print(i)
  }
}

shuffle.dat.allpest.snp<- as.data.table(shuffle.dat.allpest.snp)
shuffle.dat.allpest.kg<- as.data.table(shuffle.dat.allpest.kg)

###NOTE: SAVE WORKSPACE FOR FUTURE ANALYSIS

setwd("")
save.image(file="Workspaces/Step8_Post_Bootstrapping_Runs.RData")
#load("Workspaces/Step8_Post_Bootstrapping_Runs.RData")

############################################################################
###Assess significant runs

#Significant random associations SNPs/Pesticide
length(which(shuffle.dat.allpest.snp[Disease %in% "Alzheimer's Disease"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.snp[Disease %in% "Alzheimer's Disease"]$Prevp_pe < 0.05))

length(which(shuffle.dat.allpest.snp[Disease %in% "Brain Neoplasms"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.snp[Disease %in% "Epilepsy"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.snp[Disease %in% "Migraine Disorders"]$Prevp_pe < 0.05))

length(which(shuffle.dat.allpest.snp[Disease %in% "Multiple Sclerosis"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.snp[Disease %in% "Multiple Sclerosis"]$Prevp_pe < 0.05))

length(which(shuffle.dat.allpest.snp[Disease %in% "Parkinson Disease"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.snp[Disease %in% "Parkinson Disease"]$Prevp_pe < 0.05))

#Correlation coeffs
length(which(shuffle.dat.allpest.snp[Disease %in% "Alzheimer's Disease"]$IncCoeff_pe > 
               orig.stats[Disease %in% "Alzheimer's Disease"]$IncCoeff_pe))
length(which(shuffle.dat.allpest.snp[Disease %in% "Alzheimer's Disease"]$PrevCoeff_pe > 
               orig.stats[Disease %in% "Alzheimer's Disease"]$PrevCoeff_pe))

length(which(shuffle.dat.allpest.snp[Disease %in% "Brain Neoplasms"]$IncCoeff_pe > 
               orig.stats[Disease %in% "Brain Neoplasms"]$IncCoeff_pe)) #15% - won't include
length(which(shuffle.dat.allpest.snp[Disease %in% "Epilepsy"]$IncCoeff_pe > 
               orig.stats[Disease %in% "Epilepsy"]$IncCoeff_pe)) #50% -won't include
length(which(shuffle.dat.allpest.snp[Disease %in% "Migraine Disorders"]$PrevCoeff_pe < 
               orig.stats[Disease %in% "Migraine Disorders"]$PrevCoeff_pe)) #60%

length(which(shuffle.dat.allpest.snp[Disease %in% "Multiple Sclerosis"]$IncCoeff_pe > 
               orig.stats[Disease %in% "Multiple Sclerosis"]$IncCoeff_pe))
length(which(shuffle.dat.allpest.snp[Disease %in% "Multiple Sclerosis"]$PrevCoeff_pe > 
               orig.stats[Disease %in% "Multiple Sclerosis"]$PrevCoeff_pe))

length(which(shuffle.dat.allpest.snp[Disease %in% "Parkinson Disease"]$IncCoeff_pe > 
               orig.stats[Disease %in% "Parkinson Disease"]$IncCoeff_pe))
length(which(shuffle.dat.allpest.snp[Disease %in% "Parkinson Disease"]$PrevCoeff_pe > 
               orig.stats[Disease %in% "Parkinson Disease"]$PrevCoeff_pe))

#MAE comparison random associations SNPs/Pesticide
length(which(shuffle.dat.allpest.snp[Disease %in% "Alzheimer's Disease"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Alzheimer's Disease"]$RLM_IncRMSE))
length(which(shuffle.dat.allpest.snp[Disease %in% "Alzheimer's Disease"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Alzheimer's Disease"]$RLM_PrevRMSE)) #<7%

length(which(shuffle.dat.allpest.snp[Disease %in% "Brain Neoplasms"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Brain Neoplasms"]$RLM_IncRMSE)) #>20%
length(which(shuffle.dat.allpest.snp[Disease %in% "Epilepsy"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Epilepsy"]$RLM_IncRMSE)) #almost all
length(which(shuffle.dat.allpest.snp[Disease %in% "Migraine Disorders"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Migraine Disorders"]$RLM_PrevRMSE)) #almost all

length(which(shuffle.dat.allpest.snp[Disease %in% "Multiple Sclerosis"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Multiple Sclerosis"]$RLM_IncRMSE))
length(which(shuffle.dat.allpest.snp[Disease %in% "Multiple Sclerosis"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Multiple Sclerosis"]$RLM_PrevRMSE))

length(which(shuffle.dat.allpest.snp[Disease %in% "Parkinson Disease"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Parkinson Disease"]$RLM_IncRMSE))
length(which(shuffle.dat.allpest.snp[Disease %in% "Parkinson Disease"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Parkinson Disease"]$RLM_PrevRMSE))


#Significant random associations kg pesticide
length(which(shuffle.dat.allpest.kg[Disease %in% "Alzheimer's Disease"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.kg[Disease %in% "Alzheimer's Disease"]$Prevp_pe < 0.05))

length(which(shuffle.dat.allpest.kg[Disease %in% "Multiple Sclerosis"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.kg[Disease %in% "Multiple Sclerosis"]$Prevp_pe < 0.05))

length(which(shuffle.dat.allpest.kg[Disease %in% "Parkinson Disease"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.kg[Disease %in% "Parkinson Disease"]$Prevp_pe < 0.05))

#MAE comparison random associations SNPs/Pesticide
length(which(shuffle.dat.allpest.kg[Disease %in% "Alzheimer's Disease"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Alzheimer's Disease"]$RLM_IncRMSE))
length(which(shuffle.dat.allpest.kg[Disease %in% "Alzheimer's Disease"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Alzheimer's Disease"]$RLM_PrevRMSE))

length(which(shuffle.dat.allpest.kg[Disease %in% "Multiple Sclerosis"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Multiple Sclerosis"]$RLM_IncRMSE))
length(which(shuffle.dat.allpest.kg[Disease %in% "Multiple Sclerosis"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Multiple Sclerosis"]$RLM_PrevRMSE))

length(which(shuffle.dat.allpest.kg[Disease %in% "Parkinson Disease"]$RLM_IncRMSE < 
               orig.stats[Disease %in% "Parkinson Disease"]$RLM_IncRMSE))
length(which(shuffle.dat.allpest.kg[Disease %in% "Parkinson Disease"]$RLM_PrevRMSE < 
               orig.stats[Disease %in% "Parkinson Disease"]$RLM_PrevRMSE))

######################################################################################
#Generate plots

shuffle.dat.allpest.snp[, 1:12]<- lapply(shuffle.dat.allpest.snp[,1:12], function(x) as.numeric(as.character(x)))
#shuffle.dat.allpest.kg[, 1:12]<- lapply(shuffle.dat.allpest.kg[,1:12], function(x) as.numeric(as.character(x)))

#Prepare stats for plotting histograms
dat.org<- function(in.dat, rlm.in=TRUE){
  
  #melt data.table so variables and values are listed in two columns
  out.dat<- melt(in.dat, id.vars="Disease")
  if(rlm.in == TRUE){
    out.dat<- out.dat[grep("RLM_", out.dat$variable)]
  }
  #out.dat<- out.dat[Disease %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]
  out.dat$Facet<- out.dat$Disease #add a column to facet incidence and prevalence for each disease
  out.dat[grep("Inc", out.dat$variable)]$Facet<- paste0(out.dat[grep("Inc", out.dat$variable)]$Facet, "_Incidence")
  out.dat[grep("Prev", out.dat$variable)]$Facet<- paste0(out.dat[grep("Prev", out.dat$variable)]$Facet, "_Prevalence")
  out.dat
}

#Plot random SNP association data 
shuffle.dat.mae<- dat.org(shuffle.dat.allpest.snp[, 5:ncol(shuffle.dat.allpest.snp)]) #RLM values
shuffle.dat.pear<- dat.org(shuffle.dat.allpest.snp[, c(1:4, (ncol(shuffle.dat.allpest.snp)-1)), with=FALSE], 
                           rlm.in=FALSE) #correlation values

lm.orig<- dat.org(orig.stats) #Prepare original stats the same way
lm.orig$value<- as.numeric(lm.orig$value)
pe.orig<- dat.org(orig.stats[, c(1:2, ncol(orig.stats)), with=FALSE], rlm.in=FALSE)
pe.orig$value<- as.numeric(pe.orig$value)

#Pearson's correlation coefficient
r.plot<- ggplot() + geom_histogram(data=shuffle.dat.pear[grep("Coeff", shuffle.dat.pear$variable)], aes(x=value)) + 
  geom_vline(data=pe.orig[grep("Coeff", pe.orig$variable)], aes(xintercept=value)) + xlab("Pearson's r") +
  ylab("Distribution from 1000 random runs") + facet_wrap(.~Facet, nrow=6)

setwd("")
jpeg(file="SITables_Figures/FigureS18.jpeg", width=4.5, height=6, units="in", res=300)
r.plot
graphics.off()

#Figure S18 underlying data
figs18.dat<- shuffle.dat.pear[grep("Coeff", shuffle.dat.pear$variable)]
colnames(figs18.dat)[3]<- "Pearsons_Correlation_Coeff"
figs18.dat<- figs18.dat[order(Facet)]
figs18.dat$Rep<- rep(1:1000, 12)

setwd("")
fwrite(figs18.dat, "SITables_Figures/SI_FigureS18.csv")

###MAE/RMSE plots
mae.plot<- ggplot() + geom_histogram(data=shuffle.dat.mae[grep("MAE", shuffle.dat.mae$variable)], aes(x=value)) + 
  geom_vline(data=lm.orig[grep("MAE", lm.orig$variable)], aes(xintercept=value)) + xlab("MAE") +
  ylab("Distribution from 1000 random runs") + facet_wrap(.~Facet, nrow=6, scales="free")
rmse.plot<- ggplot() + geom_histogram(data=shuffle.dat.mae[grep("RMSE", shuffle.dat.mae$variable)], aes(x=value)) + 
  geom_vline(data=lm.orig[grep("RMSE", lm.orig$variable)], aes(xintercept=value)) + xlab("RMSE") +
  ylab("Distribution from 1000 random runs") + facet_wrap(.~Facet, nrow=6, scales="free")

setwd("")
jpeg(file="SITables_Figures/FigureS19.jpeg", width=8.5, height=7, units="in", res=300)
grid.arrange(mae.plot, rmse.plot, ncol=2)
graphics.off()

#Figure S19 underlying data
length(which(shuffle.dat.mae[grep("MAE", shuffle.dat.mae$variable)]$Facet != shuffle.dat.mae[grep("RMSE", shuffle.dat.mae$variable)]$Facet))
figs19.dat<- cbind(shuffle.dat.mae[grep("MAE", shuffle.dat.mae$variable), list(Disease, Facet, MAE=value)], 
                   shuffle.dat.mae[grep("RMSE", shuffle.dat.mae$variable), list(RMSE=value)])
figs19.dat<- figs19.dat[order(Facet)]
figs19.dat$Rep<- rep(1:1000, 12)

setwd("")
fwrite(figs19.dat, "SITables_Figures/SI_FigureS19.csv")

########################################################################################
#Assess random kg relationship

#Plot random kg applied data 
shuffle.dat.mae<- dat.org(shuffle.dat.allpest.kg[, 5:ncol(shuffle.dat.allpest.kg)]) #RLM values
shuffle.dat.pear<- dat.org(shuffle.dat.allpest.kg[, c(1:4, (ncol(shuffle.dat.allpest.kg)-1)), with=FALSE], rlm.in=FALSE)

lm.orig<- dat.org(orig.kg.stats) #Prepare original stats the same way
lm.orig$value<- as.numeric(lm.orig$value)
pe.orig<- dat.org(orig.kg.stats[, c(1:2, ncol(orig.kg.stats)), with=FALSE], rlm.in=FALSE)
pe.orig$value<- as.numeric(pe.orig$value)

mae.plot<- ggplot() + geom_histogram(data=shuffle.dat.mae[grep("MAE", shuffle.dat.mae$variable)], aes(x=value)) + 
  geom_vline(data=lm.orig[grep("MAE", lm.orig$variable)], aes(xintercept=value)) + xlab("MAE") +
  ylab("Distribution from 1000 random runs") + facet_wrap(.~Facet, nrow=3, scales="free")
rmse.plot<- ggplot() + geom_histogram(data=shuffle.dat.mae[grep("RMSE", shuffle.dat.mae$variable)], aes(x=value)) + 
  geom_vline(data=lm.orig[grep("RMSE", lm.orig$variable)], aes(xintercept=value)) + xlab("RMSE") +
  ylab("Distribution from 1000 random runs") + facet_wrap(.~Facet, nrow=3, scales="free")

#Mostly significant >80% runs
length(which(shuffle.dat.allpest.kg[Disease %in% "Multiple Sclerosis"]$Incp_pe < 0.05))
length(which(shuffle.dat.allpest.kg[Disease %in% "Multiple Sclerosis"]$Prevp_pe < 0.05))
