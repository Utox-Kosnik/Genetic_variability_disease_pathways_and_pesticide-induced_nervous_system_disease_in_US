############################################################################

### Pesticide prioritization, bubble plots, and application over time
### Bottom states fixed
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())

library(data.table)
library(tidyverse)
library(grid)
library(gridExtra)
library (plyr)
library(readxl)
library(usa)

###NOTE: LOAD OWN WORKSPACE FROM STEP 12b, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/Fixed_PriorityList_Workspace.RData")

rm(list=setdiff(ls(), c("pest.high.cv", "snp.high.cv")))

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/Analysis_Workspace.RData")
usgs.year.pest$Compound_Simp<- gsub(",", "", usgs.year.pest$COMPOUND)
chem.snp.dis.dat<- chem.snp.dis.dat[UMLS_Name %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]

###NOTE: LOAD OWN DATA FROM STEP 10, OR WORK FROM EXISTING DATASET (BELOW)
top.bot.dat<- fread("Generated_Data/Top_Bottom_States_MeanSNPs.csv")

mba.rules<- readRDS("Generated_Data/MBA_Reduced_PestYears.rds") #from step 11

length(unique(pest.high.cv$CASRN))

#Need state size for kg/mi2
state.size<- as.data.table(cbind(state.name, state.area, state.x19))
usgs.year.pest$Area<- as.numeric(state.size[match(usgs.year.pest$State, state.size$state.name)]$state.area)

############################################################################
###Generate plots using data across all ten top states

#Because priority SNPs/genes were determined per run - can't use those in pest.high.cv
#Determine priority SNPs/genes and rules implicated by each pesticide

#First determine associated SNPs/genes
pest.rank.plot<- merge(pest.high.cv[, list(CASRN, Disease, Reps, RankOrder)],
                       chem.snp.dis.dat[, list(CASRN, ChemicalName, Disease=UMLS_Name, snpId, GeneID)],
                       by=c("CASRN", "Disease"), allow.cartesian=TRUE)
#Modified one chemical name in pesticide prioritization code - wait til end to fix 

#Determine which SNPs/Genes are priority
pest.rank.plot<- merge(pest.rank.plot, snp.high.cv[, list(snpId, GeneID)], by=c("snpId", "GeneID"))

#Determine number of SNPs and genes per pesticide
pest.rank.plot<- pest.rank.plot[, .(SNPs=length(unique(na.omit(snpId))), Genes=length(unique(GeneID))),
                                by=.(Disease, CASRN, ChemicalName, RankOrder, Reps)]

############################################################################
###Determine associated rules

#Separate rules into individual compounds to match with diseases
mba.rules$Indiv_Pest<- mba.rules$Pest_Set
mba.rules<- as.data.table(separate_rows(mba.rules, Indiv_Pest, sep=","))

#ID rules in top states
mba.rules$CASRN<- usgs.year.pest[match(mba.rules$Indiv_Pest, usgs.year.pest$Compound_Simp)]$CASRN
mba.rules<- merge(mba.rules[, list(State, Pest_Set, Indiv_Pest, CASRN)], 
                  top.bot.dat[Label == 2, list(State=location_name, Disease)], by="State")

#ID top state rules implicating priority pesticides
pest.rank.plot<- merge(pest.rank.plot, mba.rules[, list(Pest_Set, CASRN, Disease)], by=c("Disease", "CASRN"))
pest.rank.plot<- pest.rank.plot[, .(Rules=length(unique(Pest_Set))), 
                                by=.(Disease, CASRN, ChemicalName, RankOrder, SNPs, Genes, Reps)]
pest.rank.plot[CASRN %in% "120116-88-3"]$ChemicalName<- "Cyazofamid" #Give shorter name
#Fix long chemical name

############################################################################
###Prepare output file

#Clean columns for output
pest.rank.out<- pest.rank.plot[order(Disease, RankOrder)]
pest.rank.out$Value<- 1
pest.rank.out[Disease %in% "Alzheimer's Disease"]$Value<- 
  rank(pest.rank.out[Disease %in% "Alzheimer's Disease"]$RankOrder, ties.method=c("average"))
pest.rank.out[Disease %in% "Multiple Sclerosis"]$Value<- 
  rank(pest.rank.out[Disease %in% "Multiple Sclerosis"]$RankOrder, ties.method=c("average"))
pest.rank.out[Disease %in% "Parkinson Disease"]$Value<- 
  rank(pest.rank.out[Disease %in% "Parkinson Disease"]$RankOrder, ties.method=c("average"))

pest.rank.out<- pest.rank.out[, list(Disease, Rank=Value, CASRN, ChemicalName, SNPs, Genes, Rules, Reps)]
pest.rank.out<- pest.rank.out[order(Disease, Rank)]

setwd("")
fwrite(pest.rank.out, "SITables_Figures/Dataset7_Pesticide_Ranks.csv")

############################################################################
###Compare to other preparation

pest.rank.comp<- fread("SITables_Figures/Dataset6_Pesticide_Ranks.csv")

wilcox.comp<- function(cv.dat, full.dat, dis.in){
  
  #adjust rank of values to consecutive numbers
  cv.dat<- cv.dat[Disease %in% dis.in]
  cv.dat$Rank<- 1:nrow(cv.dat)
  
  #adjust rank of values to consecutive numbers
  full.dat<- full.dat[Disease %in% dis.in]
  full.dat$Rank<- 1:nrow(full.dat)
  
  test.orig<- cv.dat[variable %in% full.dat$variable]
  test.cv<- full.dat[match(cv.dat$variable, full.dat$variable)]
  
  wilcox.test(test.orig$Rank, test.cv$Rank)
  
}

wilcox.comp(pest.rank.out[, list(variable=CASRN, Disease)],
            pest.rank.comp[, list(variable=CASRN, Disease)], "Alzheimer's Disease")
wilcox.comp(pest.rank.out[, list(variable=CASRN, Disease)],
            pest.rank.comp[, list(variable=CASRN, Disease)], "Multiple Sclerosis")
wilcox.comp(pest.rank.out[, list(variable=CASRN, Disease)],
            pest.rank.comp[, list(variable=CASRN, Disease)], "Parkinson Disease")
