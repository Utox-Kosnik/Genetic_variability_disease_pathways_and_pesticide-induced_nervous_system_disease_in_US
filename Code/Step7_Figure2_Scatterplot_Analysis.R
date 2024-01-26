############################################################################

### Determine association b/w disease occurrence, SNPs, and pesticides
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
require(data.table)
require(tidyverse)
require(plyr)
require(grid)
require(gridExtra)
require(MASS)
require(caret)

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
year.select<- function(chem.dat, pest.dat, year.in){
  
  #Select pesticides applied in or before the year of study
  pest.list<- pest.dat[YEAR <= year.in, list(YEAR, COMPOUND, KG_APPLIED, State, CASRN, Region)]
  pest.list<- pest.list[, .(YearsApplied=length(unique(YEAR)), KG_APPLIED=sum(KG_APPLIED)), 
                        by=.(COMPOUND, State, CASRN, Region)]
  pest.list<- pest.list[!duplicated(pest.list)]
  
  in.full.dat<- chem.dat[, .(.N), by=.(CASRN, ChemicalName, UMLS_Name, GeneID, Gene_Match, snpId)]
  usgs.year<- in.full.dat[CASRN %in% unique(pest.list$CASRN)]
  
  #Merge chem-snp-disease links with pesticide data per year
  usgs.year<- merge(usgs.year, pest.list, by="CASRN", allow.cartesian=TRUE)
  usgs.year #Now have Pesticide - SNP - Disease links per state per year
}

############################################################################
###Function for correlation data prep

dis.cor.prep<- function(cpd.dat, pest.in.dat, gbd.dat.in, year.in, measure.type, var.type="SNP", method.type="Mean"){
  
  cpd.dat<- year.select(cpd.dat, pest.in.dat, year.in)
  
  if(var.type == "SNP"){
    
    #How many SNPs implicated per pesticide
    cpd.dat<- cpd.dat[, .(Var=length(unique(na.omit(snpId)))), by=.(State, UMLS_Name, CASRN, Region)] #SNPs per pest/state
    pest.in.dat<- pest.in.dat[YEAR <= year.in, list(YEAR, State, CASRN, Region)]
    #Add years
    cpd.dat<- merge(cpd.dat, pest.in.dat, by=c("State", "CASRN", "Region"), allow.cartesian=TRUE)
    #How many SNPs implicated per pesticide-year
  }
  
  if(var.type == "Gene"){
    
    #How many gene hits per pesticide - Disease per year?
    cpd.dat<- cpd.dat[, .(Var=length(unique(GeneID))), by=.(State, UMLS_Name, CASRN, Region)]
    pest.in.dat<- pest.in.dat[YEAR <= year.in, list(YEAR, State, CASRN, Region)]
    #Add years
    cpd.dat<- merge(cpd.dat, pest.in.dat, by=c("State", "CASRN", "Region"), allow.cartesian=TRUE)
  }
  
  if(var.type == "Chem"){
    #How many pesticide hits per disease per year?
    pest.in.dat<- pest.in.dat[YEAR <= year.in, list(YEAR, State, CASRN, Region)]
    cpd.dat<- merge(cpd.dat, pest.in.dat, by=c("State", "CASRN", "Region"), allow.cartesian=TRUE)
    cpd.dat<- cpd.dat[, .(Var=length(unique(CASRN))), by=.(State, UMLS_Name, Region, YEAR)]
  }
  
  if(var.type == "KG"){
    
    #How many kg pesticide applied per state per year
    pest.in.dat<- pest.in.dat[YEAR <= year.in, list(YEAR, State, CASRN, Var=KG_APPLIED, Region)]
    cpd.dat<- merge(cpd.dat, pest.in.dat, by=c("State", "CASRN", "Region"), allow.cartesian=TRUE)
    cpd.dat<- cpd.dat[, .(.N), by=.(Var, State, CASRN, UMLS_Name, Region, YEAR)]
    cpd.dat<- cpd.dat[, .(Var=sum(Var)), by=.(State, UMLS_Name, Region, YEAR)]
  }
  
  if(method.type=="Sum"){
    
    #Add values (e.g., total SNP hits)
    cpd.dat<- cpd.dat[, .(Var=sum(Var)), by=.(location_name=State, Disease=UMLS_Name, Region)]
    
  } 
  
  if(method.type=="Mean"){
    
    #Average values (e.g., SNP hits/pesticide-year)
    cpd.dat<- cpd.dat[, .(Var=mean(Var)), by=.(location_name=State, Disease=UMLS_Name, Region)]
  }
  
  #Join with GBD disease occurrence data
  gbd.dat.plot<- gbd.dat.in[measure_name %in% measure.type & year == year.in]
  
  plot.dat<- merge(cpd.dat, gbd.dat.plot, by=c("location_name", "Disease"), allow.cartesian=TRUE)
  plot.dat
}

############################################################################
###Function for plotting correlations

#Shape and color based on geographic region
dis.cor.plot<- function(plot.in, year.in, measure.type, cpd.lab, vert.plot=FALSE){
  
  if(vert.plot==TRUE){
    #Orient plots vertically - as in figure 2
    plot.out<- ggplot(data = plot.in, aes(x = Var, y = val)) + geom_point(aes(color=Region, shape=Region)) + 
      theme_bw() + 
      ylab(paste("Disease", measure.type, "Rate")) + xlab(cpd.lab) + 
      scale_color_manual(values=c("West"="tomato","Southwest"="yellow", "Midwest"="springgreen3", 
                                  "Southeast"="skyblue2", 
                                  "Northeast"="mediumorchid2")) + guides(color=guide_legend(nrow=2)) +
      scale_x_continuous(expand=c(.1,.1)) + scale_y_continuous(expand=c(.1,.1)) +
      facet_wrap(.~Disease, scales="free", ncol=1) + geom_smooth(method = "rlm") + #scale_y_continuous(trans="log10") +
      #ggtitle(paste(cpd.lab, "vs Disease", measure.type, "Rate")) + 
      theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"),
            panel.grid = element_blank(), legend.position = "bottom")
    
  } else{
    #Orient plots horizontally - as in SI
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
    
  }
  
  return(plot.out)
}  

############################################################################
###Function for pearson's correlation calculations

#Input incidence and prevalence and determine relationships statistically
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
###SNP hits per pesticide per year

#Prepare data for correlation calculation
inc.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence")
prev.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence")

cor.calc(inc.mean.2018, prev.mean.2018)
#Reports pearson's correlation coefficient and p value (first 4 columns)
#RLM slope and intercept, MAE, RMSE (columns 5 - 12)
#Inc = incidence, Prev = prevalence

#Prepare plots of data
inc.mean.plot.2018<- dis.cor.plot(inc.mean.2018, year.in=2018, measure.type="Incidence", "SNP Hits/Pesticide")
prev.mean.plot.2018<- dis.cor.plot(prev.mean.2018, year.in=2018, measure.type="Prevalence", "SNP Hits/Pesticide")

setwd("")
jpeg("SITables_Figures/FigureS11.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(inc.mean.plot.2018, prev.mean.plot.2018, nrow=2)
graphics.off()

#Figure S11 underlying data
length(which(inc.mean.2018$location_name != prev.mean.2018$location_name))
length(which(inc.mean.2018$Disease != prev.mean.2018$Disease))
length(which(inc.mean.2018$Var != prev.mean.2018$Var))
figs11.dat<- cbind(inc.mean.2018[, list(State=location_name, Region, Disease, SNPHits_Pest=Var, 
                                      Incidence=val)], prev.mean.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs11.dat, "SITables_Figures/SI_FigureS11.csv")

#Figure 2
inc.mean.2018<- inc.mean.2018[Disease %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]
prev.mean.2018<- prev.mean.2018[Disease %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]
inc.mean.plot.2018<- dis.cor.plot(inc.mean.2018, year.in=2018, measure.type="Incidence", "SNP Hits/Pesticide",
                                  vert.plot=TRUE)
prev.mean.plot.2018<- dis.cor.plot(prev.mean.2018, year.in=2018, measure.type="Prevalence", "SNP Hits/Pesticide",
                                   vert.plot=TRUE)

#Removing title for ease of plotting
setwd("")
jpeg("Figure2.jpeg", width=3.25, height=5, units="in", res=300)
grid.arrange(inc.mean.plot.2018, prev.mean.plot.2018, ncol=2)
graphics.off()

#Figure 2 underlying data
length(which(inc.mean.2018$location_name != prev.mean.2018$location_name))
length(which(inc.mean.2018$Disease != prev.mean.2018$Disease))
length(which(inc.mean.2018$Var != prev.mean.2018$Var))
fig2.dat<- cbind(inc.mean.2018[, list(State=location_name, Region, Disease, SNPHits_Pest=Var, 
                                      Incidence=val)], prev.mean.2018[, list(Prevalence=val)])


setwd("")
fwrite(fig2.dat, "SITables_Figures/Figure2.csv")

############################################################################
###Gene hits per pesticide

g.inc.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", var.type="Gene")
g.prev.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", var.type="Gene")

cor.calc(g.inc.mean.2018, g.prev.mean.2018)

g.inc.mean.plot.2018<- dis.cor.plot(g.inc.mean.2018, year.in=2018, measure.type="Incidence", "Gene Hits/Pesticide")
g.prev.mean.plot.2018<- dis.cor.plot(g.prev.mean.2018, year.in=2018, measure.type="Prevalence", "Gene Hits/Pesticide")

setwd("")
jpeg("SITables_Figures/FigureS12.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(g.inc.mean.plot.2018, g.prev.mean.plot.2018, nrow=2)
graphics.off()

#Figure S12 underlying data
length(which(g.inc.mean.2018$location_name != g.prev.mean.2018$location_name))
length(which(g.inc.mean.2018$Disease != g.prev.mean.2018$Disease))
length(which(g.inc.mean.2018$Var != g.prev.mean.2018$Var))
figs12.dat<- cbind(g.inc.mean.2018[, list(State=location_name, Region, Disease, GeneHits_Pest=Var, 
                                        Incidence=val)], g.prev.mean.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs12.dat, "SITables_Figures/SI_FigureS12.csv")

############################################################################
###KG pesticide applied

kg.inc.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", var.type="KG", method.type="Sum")
kg.prev.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", var.type="KG", method.type="Sum")

cor.calc(kg.inc.2018, kg.prev.2018)

kg.inc.plot.2018<- dis.cor.plot(kg.inc.2018, year.in=2018, measure.type="Incidence", "Total KG Applied")
kg.prev.plot.2018<- dis.cor.plot(kg.prev.2018, year.in=2018, measure.type="Prevalence", "Total KG Applied")

setwd("")
jpeg("SITables_Figures/FigureS10.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(kg.inc.plot.2018, kg.prev.plot.2018, nrow=2)
graphics.off()

#Figure S10 underlying data
length(which(kg.inc.2018$location_name != kg.prev.2018$location_name))
length(which(kg.inc.2018$Disease != kg.prev.2018$Disease))
length(which(kg.inc.2018$Var != kg.prev.2018$Var))
figs10.dat<- cbind(kg.inc.2018[, list(State=location_name, Region, Disease, KG_Pest=Var, 
                                          Incidence=val)], kg.prev.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs10.dat, "SITables_Figures/SI_FigureS10.csv")

############################################################################
###SNPs/pesticide, not hits

cpd.ot.dat<- year.select(chem.snp.dis.dat, usgs.year.pest, 2018)
cpd.ot.dat<- cpd.ot.dat[, .(Var=length(unique(na.omit(snpId)))), by=.(CASRN, location_name=State, Region, Disease=UMLS_Name)]
cpd.ot.dat<- cpd.ot.dat[, .(Var=mean(Var)), by=.(location_name, Region, Disease)]

snp.inc.2018<- merge(cpd.ot.dat, gbd.dat[measure_name %in% "Incidence" & year == "2018"],
                     by=c("location_name", "Disease"), allow.cartesian=TRUE)
snp.prev.2018<- merge(cpd.ot.dat, gbd.dat[measure_name %in% "Prevalence" & year == "2018"],
                      by=c("location_name", "Disease"), allow.cartesian=TRUE)

cor.calc(snp.inc.2018, snp.prev.2018)

snp.inc.plot.2018<- dis.cor.plot(snp.inc.2018, year.in=2018, measure.type="Incidence", "SNPs/Pesticide")
snp.prev.plot.2018<- dis.cor.plot(snp.prev.2018, year.in=2018, measure.type="Prevalence", "SNPs/Pesticide")

setwd("")
jpeg("SITables_Figures/FigureS13.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(snp.inc.plot.2018, snp.prev.plot.2018, nrow=2)
graphics.off()

#Figure S13 underlying data
length(which(snp.inc.2018$location_name != snp.prev.2018$location_name))
length(which(snp.inc.2018$Disease != snp.prev.2018$Disease))
length(which(snp.inc.2018$Var != snp.prev.2018$Var))
figs13.dat<- cbind(snp.inc.2018[, list(State=location_name, Region, Disease, SNPs_Pest=Var, 
                                          Incidence=val)], snp.prev.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs13.dat, "SITables_Figures/SI_FigureS13.csv")

############################################################################
###Total SNP hits

inc.2018<-dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", method.type="Sum")
prev.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", method.type="Sum")

cor.calc(inc.2018, prev.2018)

inc.plot.2018<- dis.cor.plot(inc.2018, year.in=2018, measure.type="Incidence", "Total SNP Hits")
prev.plot.2018<- dis.cor.plot(prev.2018, year.in=2018, measure.type="Prevalence", "Total SNP Hits")

setwd("")
jpeg("SITables_Figures/FigureS14.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(inc.plot.2018, prev.plot.2018, nrow=2)
graphics.off()

#Figure S14 underlying data
length(which(inc.2018$location_name != prev.2018$location_name))
length(which(inc.2018$Disease != prev.2018$Disease))
length(which(inc.2018$Var != prev.2018$Var))
figs14.dat<- cbind(inc.2018[, list(State=location_name, Region, Disease, SNPHits=Var, 
                                       Incidence=val)], prev.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs14.dat, "SITables_Figures/SI_FigureS14.csv")

############################################################################
###Total pesticide hits 

c.inc.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", var.type="Chem", method.type="Sum")
c.prev.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", var.type="Chem", method.type="Sum")

cor.calc(c.inc.2018, c.prev.2018)
#Exact p-values can't be computed with ties

c.inc.plot.2018<- dis.cor.plot(c.inc.2018, year.in=2018, measure.type="Incidence", "Total Pesticide Hits")
c.prev.plot.2018<- dis.cor.plot(c.prev.2018, year.in=2018, measure.type="Prevalence", "Total Pesticide Hits")

setwd("")
jpeg("SITables_Figures/FigureS15.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(c.inc.plot.2018, c.prev.plot.2018, nrow=2)
graphics.off()

#Figure S15 underlying data
length(which(c.inc.2018$location_name != c.prev.2018$location_name))
length(which(c.inc.2018$Disease != c.prev.2018$Disease))
length(which(c.inc.2018$Var != c.prev.2018$Var))
figs15.dat<- cbind(c.inc.2018[, list(State=location_name, Region, Disease, PestHits=Var, 
                                   Incidence=val)], c.prev.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs15.dat, "SITables_Figures/SI_FigureS15.csv")

#############################################################################
###Number pesticides per type

###NOTE: COLLECT HERBICIDE, INSECTICIDE, AND FUNGICIDE DATA FROM HRAC, IRAC, AND FRAC

#Import data on pesticide type
setwd("")
pest.type.dat<- fread("Pesticide_Type_CAS_Match.csv")

###Insecticides
ins.pest<- chem.snp.dis.dat[CASRN %in% pest.type.dat[Pest_Type=="Insecticide"]$CASRN]

pest.in.dat<- pest.in.dat[YEAR <= year.in, list(YEAR, State, CASRN, Region)]
#Add years of data
cpd.dat<- merge(cpd.dat, pest.in.dat, by=c("State", "CASRN", "Region"), allow.cartesian=TRUE)
cpd.dat<- cpd.dat[, .(Var=length(unique(CASRN))), by=.(location_name=State, Disease=UMLS_Name, Region)]

ins.inc.nchem.2018<- dis.cor.prep(ins.pest, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", var.type="Chem", method.type="Unique")
ins.prev.nchem.2018<- dis.cor.prep(ins.pest, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", var.type="Chem", method.type="Unique")

cor.calc(ins.inc.nchem.2018, ins.prev.nchem.2018)

ins.inc.nchem.plot.2018<- dis.cor.plot(ins.inc.nchem.2018, year.in=2018, measure.type="Incidence", "Total Insecticides")
ins.prev.nchem.plot.2018<- dis.cor.plot(ins.prev.nchem.2018, year.in=2018, measure.type="Prevalence", "Total Insecticides")

#Herbicides
herb.pest<- chem.snp.dis.dat[CASRN %in% pest.type.dat[Pest_Type=="Herbicides"]$CASRN]

herb.inc.nchem.2018<- dis.cor.prep(herb.pest, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", var.type="Chem", method.type="Unique")
herb.prev.nchem.2018<- dis.cor.prep(herb.pest, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", var.type="Chem", method.type="Unique")

cor.calc(herb.inc.nchem.2018, herb.prev.nchem.2018)

herb.inc.nchem.plot.2018<- dis.cor.plot(herb.inc.nchem.2018, year.in=2018, measure.type="Incidence", "Total Herbicides")
herb.prev.nchem.plot.2018<- dis.cor.plot(herb.prev.nchem.2018, year.in=2018, measure.type="Prevalence", "Total Herbicides")

#Fungicides
fung.pest<- chem.snp.dis.dat[CASRN %in% pest.type.dat[Pest_Type=="Fungicides"]$CASRN]

fung.inc.nchem.2018<- dis.cor.prep(fung.pest, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Incidence", var.type="Chem", method.type="Unique")
fung.prev.nchem.2018<- dis.cor.prep(fung.pest, usgs.year.pest, gbd.dat, year.in=2018, measure.type="Prevalence", var.type="Chem", method.type="Unique")

cor.calc(fung.inc.nchem.2018, fung.prev.nchem.2018)

fung.inc.nchem.plot.2018<- dis.cor.plot(fung.inc.nchem.2018, year.in=2018, measure.type="Incidence", "Total Fungicides")
fung.prev.nchem.plot.2018<- dis.cor.plot(fung.prev.nchem.2018, year.in=2018, measure.type="Prevalence", "Total Fungicides")

#Plot all pesticide types together
setwd("")
jpeg("SITables_Figures/FigureS16.jpeg", width=9, height=10, units="in", res=300)
grid.arrange(ins.inc.nchem.plot.2018, ins.prev.nchem.plot.2018, 
             herb.inc.nchem.plot.2018, herb.prev.nchem.plot.2018,
             fung.inc.nchem.plot.2018, fung.prev.nchem.plot.2018, nrow=6)
graphics.off()

#Figure S16 underlying data


length(which(ins.inc.nchem.2018$location_name != ins.prev.nchem.2018$location_name))
length(which(herb.inc.nchem.2018$Disease != herb.prev.nchem.2018$Disease))
length(which(fung.inc.nchem.2018$Var != fung.prev.nchem.2018$Var))
figs16.dat<- rbind(cbind(ins.inc.nchem.2018[, list(State=location_name, Region, Disease, PestHits=Var, 
                                           Incidence=val)], ins.prev.nchem.2018[, list(Prevalence=val, Type="Insecticide")]),
                   cbind(herb.inc.nchem.2018[, list(State=location_name, Region, Disease, PestHits=Var, 
                                           Incidence=val)], herb.prev.nchem.2018[, list(Prevalence=val, Type="Herbicide")]),
                   cbind(fung.inc.nchem.2018[, list(State=location_name, Region, Disease, PestHits=Var, 
                                           Incidence=val)], fung.prev.nchem.2018[, list(Prevalence=val, Type="Fungicide")]))


setwd("")
fwrite(figs16.dat, "SITables_Figures/SI_FigureS16.csv")

##############################################################################
###Reducing to only same pesticides implicated in low pesticide states

#To ensure data bias in pesticides implicated in high states compared to low states
#is not the driver of relationship, limit all states to the same pesticides
test<- usgs.year.pest[, .(length(unique(CASRN))), by=.(State, Region)]
test<- test[order(V1)] #Identify states with fewest pesticides
test<- usgs.year.pest[State %in% test[1:5]$State] 
usgs.year.pest<- usgs.year.pest[CASRN %in% test$CASRN]
# Reduce pesticides to only those present in the states with the fewest pesticides

length(which(unique(usgs.year.pest$CASRN) %in% unique(chem.snp.dis.dat$CASRN))) #285 -> 183 pesticides

#Now check SNP hits/pesticide-year
inc.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Incidence")
prev.mean.2018<- dis.cor.prep(chem.snp.dis.dat, usgs.year.pest, gbd.dat, year.in="2018", measure.type="Prevalence")

cor.calc(inc.mean.2018, prev.mean.2018)

inc.mean.plot.2018<- dis.cor.plot(inc.mean.2018, year.in=2018, measure.type="Incidence", "SNP Hits/Pesticide")
prev.mean.plot.2018<- dis.cor.plot(prev.mean.2018, year.in=2018, measure.type="Prevalence", "SNP Hits/Pesticide")

setwd("")
jpeg("SITables_Figures/FigureS17.jpeg", width=9, height=3.5, units="in", res=300)
grid.arrange(inc.mean.plot.2018, prev.mean.plot.2018, nrow=2)
graphics.off()

#Figure S17 underlying data
length(which(inc.mean.2018$location_name != prev.mean.2018$location_name))
length(which(inc.mean.2018$Disease != prev.mean.2018$Disease))
length(which(inc.mean.2018$Var != prev.mean.2018$Var))
figs17.dat<- cbind(inc.mean.2018[, list(State=location_name, Region, Disease, PestHits=Var, 
                                     Incidence=val)], prev.mean.2018[, list(Prevalence=val)])


setwd("")
fwrite(figs17.dat, "SITables_Figures/SI_FigureS17.csv")
