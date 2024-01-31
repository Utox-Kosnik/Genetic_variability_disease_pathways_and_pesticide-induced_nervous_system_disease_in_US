############################################################################

### Comparison of modes of action across association rules implicated in top/bottom states
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library (plyr)
library(readxl)

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/Analysis_Workspace.RData")
usgs.year.pest$Compound_Simp<-  gsub(",", "", usgs.year.pest$COMPOUND) #Needed to match with rule data

###NOTE: LOAD OWN DATA FROM STEP 10, OR WORK FROM EXISTING DATASET (BELOW)
top.bot.dat<- fread("Generated_Data/Top_Bottom_States_MeanSNPs.csv")

mba.rules<- readRDS("Generated_Data/MBA_Reduced_PestYears.rds") #from step 11

#Separate rules into individual compounds to match with diseases
mba.rules$Indiv_Pest<- mba.rules$Pest_Set
mba.rules<- as.data.table(separate_rows(mba.rules, Indiv_Pest, sep=","))
length(unique(mba.rules$Pest_Set)) #rules across all 48 states

#Depends on MBA data, IDs unique pesticides in both top and bottom states
mba.dis.sep<- function(top.bot.in, dis.in){
  
  #For each disease, identify top and bottom state rules
  dat.r<- mba.rules[State %in% top.bot.in[Disease %in% dis.in & Label ==1]$location_name]
  dat.r$Label<- 1
  dat.out<- mba.rules[State %in% top.bot.in[Disease %in% dis.in & Label ==2]$location_name]
  dat.out$Label<- 2
  
  dat.out<- rbind(dat.out, dat.r)
  
  #Determine how many rules implicate each pesticide per state
  un.dat<- dat.out[, .(Rules=.N), by=.(State, Indiv_Pest, Label)]
  
  #Determine if these pesticides are uniquely implicated in top or bottom states 
  un.dat.out<- un.dat[, .(State=length(unique(State)), Label=length(unique(Label))), by=.(Indiv_Pest)]
  un.dat.out<- un.dat.out[Label == 1]
  un.dat<- un.dat[Indiv_Pest %in% un.dat.out$Indiv_Pest]
  
  #For each pesticide, note if it is unique to top or bottom states
  dat.out$Un_Chem<- un.dat[match(dat.out$Indiv_Pest, un.dat$Indiv_Pest)]$Label
  dat.out[is.na(Un_Chem)==TRUE]$Un_Chem<- "X"
  
  #Determine how many states that pesticide has rules in
  dat.out$Un_States<- un.dat.out[match(dat.out$Indiv_Pest, un.dat.out$Indiv_Pest)]$State
  dat.out[is.na(Un_States)==TRUE]$Un_States<- 0
  dat.out$Disease<- dis.in
  dat.out$RuleLength<- as.character(dat.out$RuleLength)
  dat.out
}

#Now have rules matched to diseases based on state (note: pesticides may not be implicated in diseases)
mba.alz<- mba.dis.sep(top.bot.dat, "Alzheimer's Disease")
mba.ms<- mba.dis.sep(top.bot.dat, "Multiple Sclerosis")
mba.pd<- mba.dis.sep(top.bot.dat, "Parkinson Disease")

#Pesticides implicated only in top states
unique(mba.alz[Un_Chem == 2]$Indiv_Pest)
unique(mba.ms[Un_Chem == 2]$Indiv_Pest)
unique(mba.pd[Un_Chem == 2]$Indiv_Pest)

##################################################################################
###Now join to actual pesticides implicated in these diseases

chem.match<- chem.snp.dis.dat[, .(.N), by=.(UMLS_Name, CASRN, ChemicalName)]
chem.match$Indiv_Pest<- usgs.year.pest[match(chem.match$CASRN, usgs.year.pest$CASRN)]$Compound_Simp

mba.alz<- merge(mba.alz, chem.match[UMLS_Name %in% "Alzheimer's Disease"], by="Indiv_Pest")
mba.ms<- merge(mba.ms, chem.match[UMLS_Name %in% "Multiple Sclerosis"], by="Indiv_Pest")
mba.pd<- merge(mba.pd, chem.match[UMLS_Name %in% "Parkinson Disease"], by="Indiv_Pest")

##################################################################################
###Generate Venn diagrams of rules

require(ggvenn)

test<- top.bot.dat[Label != 0]
test<- test[, .(length(unique(Label))), by=.(location_name)]
#One state in high one disease, low another - add twice for this analysis across diseases

all.venn<- mba.rules[State %in% test$location_name]
all.venn<- all.venn[!State %in% "Utah"]
label.match<- top.bot.dat[Label != 0]
all.venn$Label<- label.match[match(all.venn$State, label.match$location_name)]$Label
all.venn.res<- all.venn[, .(.N), by=.(Pest_Set, Label)]

utah.dat<- rbind(cbind(mba.rules[State %in% "Utah", .(.N), by=.(Pest_Set)], 2), 
                 cbind(mba.rules[State %in% "Utah", .(.N), by=.(Pest_Set)], 1))
all.venn<- rbind(all.venn.res, utah.dat[, list(Pest_Set, Label=V2, N)])
all.venn<- all.venn[, .(.N), by=.(Pest_Set, Label)]
all.venn[, .(length(unique(Pest_Set))), by=.(Label)]


#Compare association rules in high probability (Label = 2) and low probability (Label = 1) states 
all.venn.plot<- list(A=all.venn[Label == "2"]$Pest_Set, B=all.venn[Label == "1"]$Pest_Set)
names(all.venn.plot)<- c("High States", "Low States")
all.venn.plot<- ggvenn(all.venn.plot, fill_color = c("red", "blue"), set_name_size = 4)

setwd("")
jpeg("Figure3a.jpeg", width=5, height=2.5, units="in", res=300)
all.venn.plot
graphics.off()

#Data underlying Figure 3a
all.venn.print<- all.venn[, .(Label = sum(unique(Label))), by=.(Pest_Set)]

###Calculate overlap
all.venn.calc<- all.venn[, .(sum(unique(Label))), by=.(Pest_Set)]

#Set background as all unique association rules
phyper(q=length(which(all.venn.calc$V1==3)), 
       m=length(which(all.venn.calc$V1 == 1 | all.venn.calc$V1 == 3)),
       n=(length(unique(mba.rules$Pest_Set)) - length(which(all.venn.calc$V1 == 1 | all.venn.calc$V1 == 3))), 
       k=(length(which(all.venn.calc$V1 == 2 | all.venn.calc$V1 == 3))), lower.tail = FALSE)
#Do not overlap (p == 1)

##################################################################################
###Now generate per disease venn diagram

#AD
alz.venn<- mba.alz[, .(.N), by=.(Pest_Set, Label)]

#Compare association rules implicating pesticides in Alzheimer's 
#in high probability (Label = 2) and low probability (Label = 1) states 
alz.venn.plot<- list(A=alz.venn[Label == "2"]$Pest_Set, B=alz.venn[Label == "1"]$Pest_Set)
names(alz.venn.plot)<- c("High States", "Low States")
alz.venn.plot<- ggvenn(alz.venn.plot, fill_color = c("red", "blue"), set_name_size = 4)

alz.venn.calc<- alz.venn[, .(sum(unique(Label))), by=.(Pest_Set)]
phyper(q=length(which(alz.venn.calc$V1==3)), 
       m=length(which(alz.venn.calc$V1 == 1 | alz.venn.calc$V1 == 3)),
       n=(length(unique(mba.rules$Pest_Set)) - length(which(alz.venn.calc$V1 == 1 | alz.venn.calc$V1 == 3))), 
       k=(length(which(alz.venn.calc$V1 == 2 | alz.venn.calc$V1 == 3))), lower.tail = FALSE)
#Do not overlap (p == 1)


#MS
ms.venn<- mba.ms[, .(.N), by=.(Pest_Set, Label)]

#Compare association rules implicating pesticides in MS 
#in high probability (Label = 2) and low probability (Label = 1) states 
ms.venn.plot<- list(A=ms.venn[Label == "2"]$Pest_Set, B=ms.venn[Label == "1"]$Pest_Set)
names(ms.venn.plot)<- c("High States", "Low States")
ms.venn.plot<- ggvenn(ms.venn.plot, fill_color = c("red", "blue"), set_name_size = 8)

ms.venn.calc<- ms.venn[, .(sum(unique(Label))), by=.(Pest_Set)]
phyper(q=length(which(ms.venn.calc$V1==3)), 
       m=length(which(ms.venn.calc$V1 == 1 | ms.venn.calc$V1 == 3)),
       n=(length(unique(mba.rules$Pest_Set)) - length(which(ms.venn.calc$V1 == 1 | ms.venn.calc$V1 == 3))), 
       k=(length(which(ms.venn.calc$V1 == 2 | ms.venn.calc$V1 == 3))), lower.tail = FALSE)

#PD
pd.venn<- mba.pd[, .(.N), by=.(Pest_Set, Label)]

#Compare association rules implicating pesticides in PD
#in high probability (Label = 2) and low probability (Label = 1) states 
pd.venn.plot<- list(A=pd.venn[Label == "2"]$Pest_Set, B=pd.venn[Label == "1"]$Pest_Set)
names(pd.venn.plot)<- c("High States", "Low States")
pd.venn.plot<- ggvenn(pd.venn.plot, fill_color = c("red", "blue"), set_name_size = 8)

pd.venn.calc<- pd.venn[, .(sum(unique(Label))), by=.(Pest_Set)]
phyper(q=length(which(pd.venn.calc$V1==3)), 
       m=length(which(pd.venn.calc$V1 == 1 | pd.venn.calc$V1 == 3)),
       n=(length(unique(mba.rules$Pest_Set)) - length(which(pd.venn.calc$V1 == 1 | pd.venn.calc$V1 == 3))), 
       k=(length(which(pd.venn.calc$V1 == 2 | pd.venn.calc$V1 == 3))), lower.tail = FALSE)

#Will focus on per-disease for rules in multiple states

##################################################################################
###Look at pesticide type

###NOTE: COLLECT HERBICIDE, INSECTICIDE, AND FUNGICIDE DATA FROM HRAC, IRAC, AND FRAC
setwd("")
pest.type.dat<- fread("Pesticide_Type_CAS_Match.csv")

mba.casrn<- mba.rules
mba.casrn$CASRN<- usgs.year.pest[match(mba.casrn$Indiv_Pest, usgs.year.pest$Compound_Simp)]$CASRN
mba.casrn<- merge(mba.casrn, label.match[, list(location_name, Label)], by.x="State", by.y="location_name", 
                  allow.cartesian=TRUE)

pest.type.test<- merge(mba.casrn, pest.type.dat, by="CASRN", allow.cartesian=TRUE)
pest.type.test<- pest.type.test[, .(.N), by=.(CASRN, Indiv_Pest, Pest_Set, Pest_Type, Label)]

###Generate barplot of different pesticide types vs rule overlap
pest.type.plot.full<- mba.casrn[, .(Label = sum(unique(Label))), by=.(Pest_Set)]
pest.type.plot.full<- merge(pest.type.plot.full, pest.type.test[, list(Pest_Set, Indiv_Pest, Pest_Type)], by="Pest_Set")
pest.type.plot.full<- pest.type.plot.full[, .(Count=length(unique(Indiv_Pest))), by=.(Pest_Type, Label)]

pest.type.plot.full$Category<- "All"

##################################################################################
###Repeat w/ requirement of >=5 states per rule

###First compare venn diagrams
#AD
alz.venn<- mba.alz[, .(.N, States=length(unique(State))), by=.(Pest_Set, Label)]
alz.venn<- alz.venn[States >=5]
alz.venn.plot<- list(A=alz.venn[Label == "2"]$Pest_Set, B=alz.venn[Label == "1"]$Pest_Set)
names(alz.venn.plot)<- c("High States", "Low States")
alz.venn.plot<- ggvenn(alz.venn.plot, fill_color = c("red", "blue"), set_name_size = 8)

#MS
ms.venn<- mba.ms[, .(.N, States=length(unique(State))),, by=.(Pest_Set, Label)]
ms.venn<- ms.venn[States >=5]
ms.venn.plot<- list(A=ms.venn[Label == "2"]$Pest_Set, B=ms.venn[Label == "1"]$Pest_Set)
names(ms.venn.plot)<- c("High States", "Low States")
ms.venn.plot<- ggvenn(ms.venn.plot, fill_color = c("red", "blue"), set_name_size = 8)

#PD
pd.venn<- mba.pd[, .(.N, States=length(unique(State))), by=.(Pest_Set, Label)]
pd.venn<- pd.venn[States >=5]
pd.venn.plot<- list(A=pd.venn[Label == "2"]$Pest_Set, B=pd.venn[Label == "1"]$Pest_Set)
names(pd.venn.plot)<- c("High States", "Low States")
pd.venn.plot<- ggvenn(pd.venn.plot, fill_color = c("red", "blue"), set_name_size = 8)

setwd("")
jpeg("Figure3b.jpeg", width=5, height=7.5, units="in", res=300)
grid.arrange(alz.venn.plot, ms.venn.plot, pd.venn.plot, nrow=3)
graphics.off()

###Generate barplot of different pesticide types vs rule overlap
pest.type.plot.alz<- mba.alz[, .(.N), by=.(Pest_Set, Label)]
pest.type.plot.alz<- merge(pest.type.plot.alz, alz.venn, by=c("Pest_Set", "Label"))
pest.type.plot.alz<- pest.type.plot.alz[, .(Label = sum(unique(Label))), by=.(Pest_Set)]
pest.type.plot.alz<- merge(pest.type.plot.alz, pest.type.test[, list(Pest_Set, Indiv_Pest, Pest_Type)], by="Pest_Set")
pest.type.plot.alz<- pest.type.plot.alz[, .(Count=length(unique(Indiv_Pest))), by=.(Pest_Type, Label)]
pest.type.plot.alz$Category<- "AD"

###Barplot of different pesticide types vs rule overlap
pest.type.plot.ms<- mba.ms[, .(.N), by=.(Pest_Set, Label)]
pest.type.plot.ms<- merge(pest.type.plot.ms, ms.venn, by=c("Pest_Set", "Label"))
pest.type.plot.ms<- pest.type.plot.ms[, .(Label = sum(unique(Label))), by=.(Pest_Set)]
pest.type.plot.ms<- merge(pest.type.plot.ms, pest.type.test[, list(Pest_Set, Indiv_Pest, Pest_Type)], by="Pest_Set")
pest.type.plot.ms<- pest.type.plot.ms[, .(Count=length(unique(Indiv_Pest))), by=.(Pest_Type, Label)]
pest.type.plot.ms$Category<- "MS"

###Barplot of different pesticide types vs rule overlap
pest.type.plot.pd<- mba.pd[, .(.N), by=.(Pest_Set, Label)]
pest.type.plot.pd<- merge(pest.type.plot.pd, pd.venn, by=c("Pest_Set", "Label"))
pest.type.plot.pd<- pest.type.plot.pd[, .(Label = sum(unique(Label))), by=.(Pest_Set)]
pest.type.plot.pd<- merge(pest.type.plot.pd, pest.type.test[, list(Pest_Set, Indiv_Pest, Pest_Type)], by="Pest_Set")
pest.type.plot.pd<- pest.type.plot.pd[, .(Count=length(unique(Indiv_Pest))), by=.(Pest_Type, Label)]
pest.type.plot.pd$Category<- "PD"

pest.type.bar.plot<- rbind(pest.type.plot.full, pest.type.plot.alz, pest.type.plot.ms, pest.type.plot.pd)
pest.type.bar.plot$Label<- as.character(pest.type.bar.plot$Label)
pest.type.bar.plot$Label<- factor(pest.type.bar.plot$Label, levels=c("2", "3", "1"), ordered=TRUE)

#All pesticides across rules
all.bar<- ggplot(pest.type.bar.plot[Category == "All"]) + geom_col(aes(x=Pest_Type, y=Count, fill=Pest_Type)) +
  facet_wrap(~Label, scales="free") + 
  scale_fill_manual(values=c("Fungicides"="orange", "Herbicides"="forestgreen", "Insecticide"="midnightblue")) +
  scale_y_continuous(expand=c(0, 0), limits = c(0, 135), breaks=c(0, 50, 100)) + 
  xlab(NULL) + ylab(NULL)+
  theme_bw() + theme(text=element_text(size=8), panel.grid = element_blank(), legend.position="none", 
                     strip.background = element_blank(), strip.text.x = element_blank(),
                     axis.ticks.x = element_blank(), axis.text.x=element_blank(), 
                     axis.text.y=element_text(colour="black", face="bold"))

#Alzheimers associated pesticides across rules
alz.bar<- ggplot(pest.type.bar.plot[Category == "AD"]) + geom_col(aes(x=Pest_Type, y=Count, fill=Pest_Type)) +
  facet_wrap(~Label, scales="free") + 
  scale_fill_manual(values=c("Fungicides"="orange", "Herbicides"="forestgreen", "Insecticide"="midnightblue")) +
  scale_y_continuous(expand=c(0, 0), limits = c(0, 110), breaks=c(0, 50, 100)) + 
  xlab(NULL) + ylab(NULL)+
  theme_bw() + theme(text=element_text(size=8), panel.grid = element_blank(), legend.position="none", 
                     strip.background = element_blank(), strip.text.x = element_blank(),
                     axis.ticks.x = element_blank(), axis.text.x=element_blank(), 
                     axis.text.y=element_text(colour="black", face="bold"))

#MS associated pesticides across rules
ms.bar<- ggplot(pest.type.bar.plot[Category == "MS"]) + geom_col(aes(x=Pest_Type, y=Count, fill=Pest_Type)) +
  facet_wrap(~Label, scales="free") + 
  scale_fill_manual(values=c("Fungicides"="orange", "Herbicides"="forestgreen", "Insecticide"="midnightblue")) +
  scale_y_continuous(expand=c(0, 0), limits = c(0, 110), breaks=c(0, 50, 100)) + 
  xlab(NULL) + ylab(NULL)+
  theme_bw() + theme(text=element_text(size=8), panel.grid = element_blank(), legend.position="none", 
                     strip.background = element_blank(), strip.text.x = element_blank(),
                     axis.ticks.x = element_blank(), axis.text.x=element_blank(), 
                     axis.text.y=element_text(colour="black", face="bold"))

#PD associated pesticides across rules
pd.bar<- ggplot(pest.type.bar.plot[Category == "PD"]) + geom_col(aes(x=Pest_Type, y=Count, fill=Pest_Type)) +
  facet_wrap(~Label, scales="free") + 
  scale_fill_manual(values=c("Fungicides"="orange", "Herbicides"="forestgreen", "Insecticide"="midnightblue")) +
  scale_y_continuous(expand=c(0, 0), limits = c(0, 110), breaks=c(0, 50, 100)) + 
  xlab(NULL) + ylab(NULL)+
  theme_bw() + theme(text=element_text(size=8), panel.grid = element_blank(), legend.position="none", 
                     strip.background = element_blank(), strip.text.x = element_blank(),
                     axis.ticks.x = element_blank(), axis.text.x=element_blank(), 
                     axis.text.y=element_text(colour="black", face="bold"))


setwd("")
jpeg("Figure3_Bars.jpeg", width=2, height=2.5, units="in", res=300)
grid.arrange(all.bar, alz.bar, ms.bar, pd.bar, nrow=4)
graphics.off()

#Barplot data

pest.type.bar.plot$Label<- as.character(pest.type.bar.plot$Label)
pest.type.bar.plot[Label == "1"]$Label<- "Low States"
pest.type.bar.plot[Label == "2"]$Label<- "High States"
pest.type.bar.plot[Label == "3"]$Label<- "Both"
pest.type.bar.plot[Category == "AD"]$Category<- "Alzheimer's Disease"
pest.type.bar.plot[Category == "MS"]$Category<- "Multiple Sclerosis"
pest.type.bar.plot[Category == "PD"]$Category<- "Parkinson Disease"
colnames(pest.type.bar.plot)[3]<- "Pesticides"

setwd("")
fwrite(pest.type.bar.plot, "SITables_Figures/Figure3.csv")
