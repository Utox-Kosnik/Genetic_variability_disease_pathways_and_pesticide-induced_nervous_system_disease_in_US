############################################################################

### Frequent itemset mining of pesticides applied over years in US states
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library (plyr)
library(readxl)
library(arules)
library(arulesViz)

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)

setwd("")
load("Workspaces/Analysis_Workspace.RData")

#Only run for top/bottom states

###NOTE: LOAD OWN DATA FROM STEP 10, OR WORK FROM EXISTING DATASET (BELOW)
setwd("")
top.bot.dat<- fread("Generated_Data/Top_Bottom_States_MeanSNPs.csv")
state.run<- unique(top.bot.dat$location_name)

#Remove commas from pesticides names - helps build comma separated data
usgs.year.pest$Compound_Simp<- gsub(",", "", usgs.year.pest$COMPOUND)
length(unique(usgs.year.pest$COMPOUND))==length(unique(usgs.year.pest$Compound_Simp))

#Create the "transaction format" of data - prep one file per state
#Line up pesticide use per year as transactions per state 

setwd("/MBA_Data")

for(state.in in state.run){
  
  #For each state, list the pesticides per year
  tr.dat.int<- usgs.year.pest[State %in% state.in]
  
  #Paste pesticides per year together as one column
  tr.dat.int<- tr.dat.int[, .(paste(unique(Compound_Simp), sep="", collapse=",")), by=.(YEAR)]
  
  tr.dat.out<- as.data.table(tr.dat.int$V1)
  
  #Output transaction data
  fwrite(tr.dat.out, file=paste0(state.in, "_USGS_Pesticide_Transactions.csv"), quote=FALSE, col.names=FALSE)
}

#Clean workspace
rm(list=setdiff(ls(), "state.run"))

#Now develop one set of rules per state

mba.dat<- NULL

setwd("/MBA_Data")

#For each state, determine set of association rules
for(state.in in state.run){
  
  #Read in transaction file per state
  tr.dat<- read.transactions(paste0(state.in, "_USGS_Pesticide_Transactions.csv"), format="basket", sep=",")
  
  #Run frequent itemset mining
  ar.int<- apriori(tr.dat, parameter = list(supp=0.18, conf=0.8, minlen=2, maxlen=3, maxtime=0),
                   control = list(memopt = TRUE, load = FALSE))
  ar.int<- ar.int[!is.redundant(ar.int)] #Remove redundant rules
  
  #Arrange output
  assoc.table<- as(ar.int, "data.frame"); assoc.table<- as.data.table(assoc.table)
  
  #Organize table - remove excess characters
  assoc.table$rules<- gsub("\\{", "", assoc.table$rules); assoc.table$rules<- gsub("\\}", "", assoc.table$rules)
  assoc.table$Pests_L<- gsub(" =>.*", "", assoc.table$rules)
  assoc.table$Pests_R<- gsub(".*=> ", "", assoc.table$rules)
  assoc.table$State<- state.in
  
  #Join to overall table with each set of rules labeled by state
  mba.dat<- rbind(mba.dat, assoc.table)
  print(state.in)
}

###NOTE: SAVE DATA TO KEEP FULL COPY
setwd("")
saveRDS(mba.dat, file="Generated_Data/All_MBA_PestYears.rds")

#Analyze the overall linkages
#Remove directionality from rules (i.e. A => B,C is the same as B => A,C for our purposes)
mba.dat$Pest_Set<- gsub(" => ", ",", mba.dat$rules)

#For each rule, unlist pesticides, order them alphabetically, and rejoin
mba.dat$Pest_Set<- lapply(mba.dat$Pest_Set, function(x) as.character(unlist(strsplit(x, split=","))))
mba.dat$RuleLength<- lapply(mba.dat$Pest_Set, function(x) length(unlist(x)))
mba.dat$Pest_Set<- lapply(mba.dat$Pest_Set, function(x) paste(x[order(x)], collapse=","))
mba.dat$Pest_Set<- as.character(mba.dat$Pest_Set)
mba.dat$RuleLength<- as.character(mba.dat$RuleLength)

#Reduce rules that were kept distinct because of directionality, keep max values for each association rule
mba.dat<- mba.dat[, .(support=max(support), confidence = max(confidence), coverage=max(coverage),
                      lift=max(lift)), by=.(State, Pest_Set, count, RuleLength)]

###NOTE: SAVE DATA FOR FUTURE ANALYSIS
setwd("")
saveRDS(mba.dat, file="Generated_Data/MBA_Reduced_PestYears.rds")

############################################################################
###Prepare output table 

setwd("")
mba.dat<- readRDS("Generated_Data/MBA_Reduced_PestYears.rds")#From step 11

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)
###NOTE: LOAD OWN DATA FROM STEP 10, OR WORK FROM EXISTING DATASET (BELOW)
load("Workspaces/Analysis_Workspace.RData")
top.bot.dat<- fread("Generated_Data/Top_Bottom_States_MeanSNPs.csv")
top.bot.dat<- as.data.table(rbind(head(top.bot.dat[Label == 2 & Disease == "Alzheimer's Disease"], 5),
                                  head(top.bot.dat[Label == 2 & Disease == "Multiple Sclerosis"], 5),
                                  head(top.bot.dat[Label == 2 & Disease == "Parkinson Disease"], 5)))

usgs.year.pest$Compound_Simp<- gsub(",", "", usgs.year.pest$COMPOUND)
usgs.year.pest<- merge(usgs.year.pest, top.bot.dat[, list(location_name, Disease, Label)], 
                       by.x="State", by.y="location_name", allow.cartesian = TRUE)

chem.snp.dis.out<- chem.snp.dis.dat[UMLS_Name %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"), .(.N), by=.(Disease=UMLS_Name, CASRN, ChemicalName)]
chem.snp.dis.out<- merge(usgs.year.pest, chem.snp.dis.out, by=c("CASRN", "Disease"), allow.cartesian=TRUE)

#Break up rules into individual chemical associations
rules.out<- mba.dat[State %in% top.bot.dat$location_name]
rules.out$Compound_Simp<- rules.out$Pest_Set
rules.out<- as.data.table(separate_rows(rules.out, Compound_Simp, sep=","))

#Merge chemicals in our dataset with the rules to see how many have overlap
rules.out<- merge(rules.out, chem.snp.dis.out[, .(.N), by=.(Disease, State, Compound_Simp, ChemicalName)], 
                  by=c("State", "Compound_Simp"), allow.cartesian=TRUE)

rules.out<- rules.out[, .(Dataset_Rule=paste0(ChemicalName, collapse=";")), by=.(Disease, Rule=Pest_Set, State, Count=count,
                                                                                 RuleLength, Support=support, Confidence=confidence,
                                                                                 Lift=lift)]

rules.out<- rules.out[order(Disease)]

setwd("")
fwrite(rules.out, file="SITables_Figures/Dataset3-Association_Rules.csv")
