############################################################################

### Develop comparison datasets of states for cross validation
### Vary top and bottom states
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library (plyr)
library(readxl)
library(usa)

#For the top 10 and bottom 10 states, select 5 and repeat analysis
#Then, will compare output datasets to the original dataset

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/Analysis_Workspace.RData")
usgs.year.pest$Compound_Simp<- gsub(",", "", usgs.year.pest$COMPOUND)
chem.snp.dis.dat<- chem.snp.dis.dat[UMLS_Name %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]

state.size<- as.data.table(cbind(state.name, state.area, state.x19))
usgs.year.pest$Area<- as.numeric(state.size[match(usgs.year.pest$State, state.size$state.name)]$state.area)

###NOTE: LOAD OWN DATA FROM STEP 10, OR WORK FROM EXISTING DATASET (BELOW)
top.bot.dat<- fread("Generated_Data/Top_Bottom_States_MeanSNPs.csv")

#Identify set of top/bottom states to use for cross-validation
alz<- top.bot.dat[Disease %in% "Alzheimer's Disease" & Label != 0]
ms<- top.bot.dat[Disease %in% "Multiple Sclerosis" & Label != 0]
pd<- top.bot.dat[Disease %in% "Parkinson Disease" & Label != 0]

setwd("")
mba.rules<- readRDS("Generated_Data/MBA_Reduced_PestYears.rds") #from step 11

###Prep SNP severity table
#https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
snp.severity<- data.table(c("splice acceptor variant", "splice donor variant", "stop gained", "frameshift variant", 
                              "stop lost", "start lost",
                              "inframe insertion", "inframe deletion", "missense variant", "protein altering variant",
                              "splice region variant", "synonymous variant",
                              "coding sequence variant", "mature miRNA variant", "5 prime UTR variant", "3 prime UTR variant", 
                              "non coding transcript exon variant", "intron variant", "non coding transcript variant",
                              "upstream gene variant", "downstream gene variant", "Intergenic variant"))
snp.severity$Severity<- c("High", "High", "High", "High", "High", "High", 
                            "Moderate", "Moderate", "Moderate", "Moderate", 
                            "Low", "Low", 
                            "Modifier", "Modifier", "Modifier", "Modifier", "Modifier", "Modifier", "Modifier", "Modifier",
                            "Modifier", "Modifier")

############################################################################
###1) Select 5 random top and bottom states out of 10

state.select<- function(){
  state.int<- rbind(alz[location_name %in% c(sample(alz[Label ==2]$location_name, 5), 
                                             sample(alz[Label ==1]$location_name, 5))],
                    ms[location_name %in% c(sample(ms[Label ==2]$location_name, 5), 
                                            sample(ms[Label ==1]$location_name, 5))],
                    pd[location_name %in% c(sample(pd[Label ==2]$location_name, 5), 
                                            sample(pd[Label ==1]$location_name, 5))])
  state.int
}

############################################################################
###2) Organize rules

#Function to organize rules (SI table creation in Step 11)
rules.select<- function(state.set, usgs.dat, pest.dat, mba.dat){
  
  pest.dat<- pest.dat[State %in% state.set[Label == 2]$location_name]
  usgs.dat<- usgs.dat[, .(.N), by=.(Disease=UMLS_Name, CASRN, ChemicalName)]
  usgs.dat<- merge(pest.dat, usgs.dat, by="CASRN", allow.cartesian=TRUE)
  
  #Break up rules into individual chemical associations
  rules.int<- mba.dat
  rules.int$Compound_Simp<- rules.int$Pest_Set
  rules.int<- as.data.table(separate_rows(rules.int, Compound_Simp, sep=","))
  
  #Merge chemicals in our dataset with the rules to see how many have overlap
  rules.int<- merge(rules.int, usgs.dat[, .(.N), by=.(Disease, State, Compound_Simp, ChemicalName)], 
                    by=c("State", "Compound_Simp"), allow.cartesian=TRUE)
  
  rules.int<- rules.int[, .(Dataset_Rule=paste0(ChemicalName, collapse=";")), by=.(Disease, Rule=Pest_Set, State, Count=count,
                                                                                   RuleLength, Support=support, Confidence=confidence,
                                                                                   Lift=lift)]
  rules.int
}

############################################################################
###3) SNP priority list

#Function to sort rules per disease and distinguish top/bottom states
mba.dis.sep<- function(top.bot.in, mba.dat, dis.in){
  
  #For each disease, identify top and bottom state rules
  dat.r<- mba.dat[State %in% top.bot.in[Disease %in% dis.in & Label ==1]$location_name]
  dat.r$Label<- "Bottom"
  dat.out<- mba.dat[State %in% top.bot.in[Disease %in% dis.in & Label ==2]$location_name]
  dat.out$Label<- "Top"
  
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

#Function to ID SNP hit differences
high.snp.assess<- function(usgs.dat, year.pest, top.bot.5, dis.in){
  
  #For disease, reduce input data columns
  snp.dat.in<- usgs.dat[UMLS_Name %in% dis.in, list(CASRN, ChemicalName, GeneID, GeneSymbol, Gene_Match,
                                                    snpId, chromosome, position, variant_consequence_type)]
  snp.dat.in<- snp.dat.in[!duplicated(snp.dat.in)]
  
  #Collect only top/bottom state data
  pest.dat.in<- year.pest[State %in% top.bot.5[Label != 0 & Disease %in% dis.in]$location_name]
  
  #Merge SNP data/years
  snp.dat.in<- merge(snp.dat.in, pest.dat.in, by="CASRN", allow.cartesian=TRUE)
  
  #How many chemicals/how much KG linked with each SNP?
  snp.dat.in<- snp.dat.in[, .(Chems=length(unique(CASRN)), KG=sum(KG_APPLIED)), by=.(State, YEAR, GeneID, GeneSymbol, Gene_Match,
                                                                                     snpId, chromosome, position, variant_consequence_type)]
  
  #Label states
  snp.dat.in<- merge(snp.dat.in, top.bot.5[Disease %in% dis.in & Label != 0], 
                     by.x="State", by.y="location_name")
  
  #Determine values across years - hits = number pesticide applications across years
  snp.dat.in<- snp.dat.in[, .(Hits=sum(Chems)), by=.(Disease, Label, State, GeneID, GeneSymbol, Gene_Match,
                                                     snpId, chromosome, position, variant_consequence_type)]
  
  #Find greatest differences between top/bottom states
  
  #Determine difference across states within each category (top vs bottom)
  dif.dat<- snp.dat.in[, .(Hits=sum(Hits)),
                       by=.(Label, GeneID, GeneSymbol, Gene_Match, snpId)]
  
  #Determine differences in values between top and bottom 
  dif.dat<- dif.dat[, .(HighHits=(Hits[Label == 2]-Hits[Label ==1])),
                    by=.(GeneID, GeneSymbol, Gene_Match, snpId)]
  
  #Which case had more hits in high states than low states
  dif.dat<- dif.dat[HighHits > 0]
  
  snp.dat.out<- snp.dat.in[GeneID %in% dif.dat$GeneID] #Covers both genes and SNPs
  
  snp.dat.out<- merge(snp.dat.out, dif.dat, by=c("GeneID", "GeneSymbol", "Gene_Match", "snpId"))
  snp.dat.out[order(-Label)] #Sort with top states first
}

###Function to determine how many rules implicate ranked SNPs
rule.dat.assess<- function(dis.in, mba.in, snp.in, usgs.dat, pest.dat){
  
  gene.set<- unique(snp.in[Disease %in% dis.in]$GeneID)
  dat.out<- NULL
  
  #For each high-ranked gene, 
  for(gene.in in gene.set){
    
    #Collect gene data
    rank.in<- snp.in[Disease %in% dis.in & GeneID %in% gene.in]
    
    #How many SNPs per gene
    rank.in<- rank.in[, .(SNPs=length(unique(snpId))), by=.(GeneID, GeneSymbol, Gene_Match, Disease, Label, State,
                                                            Hits, HighHits)]
    
    #Merge with original data to collect chemical information for rules matching
    rank.in<- merge(rank.in, usgs.dat[UMLS_Name %in% dis.in, list(Disease=UMLS_Name, CASRN, ChemicalName, GeneID, GeneSymbol, Gene_Match)],
                    by=c("Disease", "GeneID", "GeneSymbol", "Gene_Match"), allow.cartesian=TRUE)
    rank.in<- rank.in[!duplicated(rank.in)]
    
    #Add matching USGS name to merge with MBA results
    rank.in$Indiv_Pest<- pest.dat[match(rank.in$CASRN, pest.dat$CASRN)]$Compound_Simp
    rank.in<- merge(rank.in, mba.in[Label == "Top", list(State, Pest_Set, count, RuleLength, Indiv_Pest, Un_Chem, 
                                                         Un_States)], by=c("State", "Indiv_Pest"), allow.cartesian=TRUE)
    
    #Now add up rules
    rule.count<- rank.in[, .(TotChems=length(unique(CASRN))), by=.(GeneID, GeneSymbol, Gene_Match, HighHits, 
                                                                   Pest_Set, count, RuleLength, Un_Chem, Un_States, State)]
    
    #FIRST
    #Rules with complete chemical coverage in SNP integrated dataset
    top.test<- rule.count[RuleLength == TotChems]
    #How many states
    top.test<- top.test[, .(State= length(unique(State))), by=.(GeneID, GeneSymbol, Gene_Match, HighHits, 
                                                                Pest_Set, RuleLength, Un_Chem)]
    top.test$Match<-paste(top.test$Pest_Set, top.test$GeneID, top.test$Gene_Match, sep="_")
    rule.count$Match<-paste(rule.count$Pest_Set, rule.count$GeneID, rule.count$Gene_Match, sep="_")
    
    rule.count$AllChem_Rule<- 0
    rule.count[Match %in% top.test[State >=3]$Match]$AllChem_Rule<- 5 #All chems in rule, many states
    
    #Second
    #Rules with incomplete chemical coverage, but at least one chemical
    ot.rule<- rule.count[AllChem_Rule == 0, .(State= length(unique(State))), by=.(GeneID, GeneSymbol, Gene_Match, HighHits, 
                                                                                  Pest_Set, RuleLength, Un_Chem, Match)]
    
    rule.count[Match %in% ot.rule[State >=3]$Match]$AllChem_Rule<- 1
    
    rule.count<- rule.count[, .(Pest_Set=length(unique(Pest_Set))), by=.(GeneID, GeneSymbol, Gene_Match, HighHits,
                                                                         TotChems, AllChem_Rule)]
    
    dat.out<- rbind(dat.out, rule.count)
    gc()
    
    if(gene.in %in% gene.set[seq(1, 3000, 50)]){
      print(which(gene.set %in% gene.in))
    }
  }
  
  return(dat.out)
}
#Top rules have all chemicals in the rule associated w/ the disease and are in 3+ states
#Middle rules have at least one chemical in the rule associated w/ the disease and are in 3+ states

snp.hit.select<- function(state.set.in, usgs.dat.in, pest.dat.in){
  
  high.snp.all<- rbind(high.snp.assess(usgs.dat.in, pest.dat.in, state.set.in, "Alzheimer's Disease"), 
                       high.snp.assess(usgs.dat.in, pest.dat.in, state.set.in, "Multiple Sclerosis"), 
                       high.snp.assess(usgs.dat.in, pest.dat.in, state.set.in, "Parkinson Disease"))
  
  high.snp.all$Severity<- snp.severity[match(high.snp.all$variant_consequence_type, snp.severity$V1)]$Severity
  high.snp.all$Severity<- factor(high.snp.all$Severity, levels=c("High", "Moderate", "Modifier", "Low"), ordered=TRUE)
  
  high.snp.all<- high.snp.all[order(Disease, Severity, -HighHits, -Label)]
  high.snp.all$Match<- paste(high.snp.all$Disease, high.snp.all$GeneID, high.snp.all$snpId, sep="_")
  
  rank.dat<- high.snp.all[, .(HighState=length(unique(State[Label == 2])),
                              LowState=length(unique(State[Label == 1]))), 
                          by=.(Disease, GeneID, GeneSymbol, Gene_Match, snpId, chromosome, position, 
                               variant_consequence_type, Severity,
                               HighHits)]
  rank.dat$Match<- paste(rank.dat$Disease, rank.dat$GeneID, rank.dat$snpId, sep="_")
  
  length(which(rank.dat$HighState > rank.dat$LowState))
  length(which(rank.dat$HighState == rank.dat$LowState))
  
  rank.dat<- rank.dat[HighState >= LowState]
  
  merge(rank.dat[, list(Disease, GeneID, GeneSymbol, Gene_Match, snpId, chromosome, position,
                        variant_consequence_type, Severity, HighHits, HighState, LowState, Match)], 
        high.snp.all[, list(Label, State, Hits, Match)], by="Match", allow.cartesian=TRUE)
}

snp.rule.select<- function(state.set.in, usgs.dat.in, pest.dat.in, mba.int, high.snp.in){
  #Separate rules into individual pesticides
  mba.int$Indiv_Pest<- mba.int$Pest_Set
  mba.int<- as.data.table(separate_rows(mba.int, Indiv_Pest, sep=","))
  
  mba.alz<- mba.dis.sep(state.set.in, mba.int, "Alzheimer's Disease")
  mba.ms<- mba.dis.sep(state.set.in, mba.int, "Multiple Sclerosis")
  mba.pd<- mba.dis.sep(state.set.in, mba.int, "Parkinson Disease")
  
  alz.rule.rank<- rule.dat.assess("Alzheimer's Disease", mba.alz, high.snp.in, usgs.dat.in, pest.dat.in)
  ms.rule.rank<- rule.dat.assess("Multiple Sclerosis", mba.ms, high.snp.in, usgs.dat.in, pest.dat.in)
  pd.rule.rank<- rule.dat.assess("Parkinson Disease", mba.pd, high.snp.in, usgs.dat.in, pest.dat.in)
  
  alz.rule.rank<- merge(alz.rule.rank, usgs.dat.in[UMLS_Name %in% "Alzheimer's Disease", (.N), by=.(GeneID, Gene_Match, snpId, chromosome, position, variant_consequence_type)],
                        by=c("GeneID", "Gene_Match"), allow.cartesian=TRUE)
  ms.rule.rank<- merge(ms.rule.rank, usgs.dat.in[UMLS_Name %in% "Multiple Sclerosis", (.N), by=.(GeneID, Gene_Match, snpId, chromosome, position, variant_consequence_type)],
                       by=c("GeneID", "Gene_Match"), allow.cartesian=TRUE)
  pd.rule.rank<- merge(pd.rule.rank, usgs.dat.in[UMLS_Name %in% "Parkinson Disease", (.N), by=.(GeneID, Gene_Match, snpId, chromosome, position, variant_consequence_type)],
                       by=c("GeneID", "Gene_Match"), allow.cartesian=TRUE)
  
  alz.rule.rank$Disease<- "Alzheimer's Disease"
  ms.rule.rank$Disease<- "Multiple Sclerosis"
  pd.rule.rank$Disease<- "Parkinson Disease"
  
  as.data.table(rbind(alz.rule.rank, ms.rule.rank, pd.rule.rank))
  
}

#Function to develop overall numerical value per SNP/gene rank
fin.snp.calc<- function(rule.dat.in, snp.dat.in){
  
  rule.rank.fin<- rule.dat.in[AllChem_Rule %in% c(1, 5)]
  rule.rank.fin<- rule.rank.fin[, .(RuleRank=sum(AllChem_Rule * Pest_Set), Rules=sum(Pest_Set)), 
                                by=.(snpId, GeneID, GeneSymbol, Gene_Match, HighHits, variant_consequence_type, Disease)]
  #Create ranking based on presence in several top states
  
  ###Adjust hit rank
  #Current rank dat are ordered by number of hits and severity of hit
  #Add rank value that accounts for ties and severity
  
  rule.rank.fin$HitRank<- as.numeric(as.character(rule.rank.fin$HighHits))
  rule.rank.fin$Severity<- snp.dat.in[match(rule.rank.fin$variant_consequence_type, snp.dat.in$variant_consequence_type)]$Severity
  table(rule.rank.fin$Severity) #Adjust ranks
  
  #Rather than rescaling, weight - keep gene hits the same
  #Leave this and modifier alone, reduce low, up moderate and high
  rule.rank.fin[is.na(Severity) == TRUE & is.na(snpId) == FALSE]$HitRank<- rule.rank.fin[is.na(Severity) == TRUE & is.na(snpId) == FALSE]$HitRank
  rule.rank.fin[Severity == "Modifier"]$HitRank<- rule.rank.fin[Severity == "Modifier"]$HitRank
  rule.rank.fin[Severity == "Low"]$HitRank<- rule.rank.fin[Severity == "Low"]$HitRank * 0.5
  rule.rank.fin[Severity == "Moderate"]$HitRank<- rule.rank.fin[Severity == "Moderate"]$HitRank * 2
  rule.rank.fin[Severity == "High"]$HitRank<- rule.rank.fin[Severity == "High"]$HitRank * 5
  
  #Overall ranking is ranking based on rules multiplied by ranking based on hits
  rule.rank.fin<- rule.rank.fin[, .(RankValue=RuleRank * HitRank), by=.(Disease, snpId, GeneID, GeneSymbol, Gene_Match, 
                                                                        Rules, HighHits, RuleRank, HitRank)]
  rule.rank.fin<- rule.rank.fin[order(Disease, -RankValue)]
  
  snp.dat.out<- merge(rule.rank.fin, snp.dat.in[, list(Disease, GeneID, Gene_Match, snpId, chromosome, position, 
                                                         SeverityOrig=Severity, HighState)], 
                      by=c("Disease", "GeneID", "Gene_Match", "snpId"))
  snp.dat.out<- snp.dat.out[!duplicated(snp.dat.out)]
  
  snp.dat.out[is.na(snpId) == FALSE & Gene_Match == ""]$GeneID<- ""
  snp.dat.out[is.na(snpId) == FALSE & Gene_Match == ""]$GeneSymbol<- ""
  
  snp.dat.out$Value<- 1
  snp.dat.out[Disease %in% "Alzheimer's Disease"]$Value<- 
    rank(1/(snp.dat.out[Disease %in% "Alzheimer's Disease"]$RankValue), ties.method=c("average"))
  snp.dat.out[Disease %in% "Multiple Sclerosis"]$Value<- 
    rank(1/(snp.dat.out[Disease %in% "Multiple Sclerosis"]$Rank), ties.method=c("average"))
  snp.dat.out[Disease %in% "Parkinson Disease"]$Value<- 
    rank(1/(snp.dat.out[Disease %in% "Parkinson Disease"]$Rank), ties.method=c("average"))
  snp.dat.out<- snp.dat.out[order(Disease, Value)]
  
  snp.dat.out[, list(Rank=Value, Disease, snpId, GeneID, GeneSymbol, HighHits, Rules)]
}

############################################################################
###4) Pathway priority list
#Not proceeding to ranked list as end output - need full set for Figure 4

path.select<- function(state.set.in, usgs.dat.in, pest.dat.in){
  
  pest.tmp<- pest.dat.in[, .(YearsApplied=length(unique(YEAR)), KG_APPLIED=sum(KG_APPLIED)), 
                            by=.(State, CASRN)]
  
  usgs.year<- usgs.dat.in[, .(.N), by=.(CASRN, ChemicalName, UMLS_Name, Gene_Match, 
                                        snpId, GeneID, GeneSymbol, PathwayID, PathwayName)]
  usgs.year<- merge(usgs.year, pest.tmp, by="CASRN", allow.cartesian=TRUE)
  
  high.path<- merge(usgs.year, state.set.in[Label == 2], by.x= c("UMLS_Name", "State"), 
                    by.y=c("Disease", "location_name"), allow.cartesian=TRUE)
  low.path<- merge(usgs.year, state.set.in[Label == 1], by.x= c("UMLS_Name", "State"), 
                   by.y=c("Disease", "location_name"), allow.cartesian=TRUE)
  
  #Per SNP/gene, determine how many pesticides and how many hits
  high.path<- high.path[, .(Hits=sum(YearsApplied), Chems=length(unique(CASRN))), 
                        by=.(Disease=UMLS_Name, State, Gene_Match, snpId, GeneID, 
                             GeneSymbol, PathwayID, PathwayName)]
  
  low.path<- low.path[, .(Hits=sum(YearsApplied), Chems=length(unique(CASRN))), 
                      by=.(Disease=UMLS_Name, State, Gene_Match, snpId, GeneID, 
                           GeneSymbol, PathwayID, PathwayName)]
  
  #Now assess total hit differences b/w high/low states
  high.path.sum<- high.path[, .(TotHits=sum(Hits), meanChems=mean(Chems), States=length(unique(State))), 
                            by=.(Disease, Gene_Match, snpId, GeneID, GeneSymbol, 
                                 PathwayID, PathwayName)]
  
  low.path.sum<- low.path[, .(TotHits=sum(Hits), meanChems=mean(Chems), States=length(unique(State))), 
                          by=.(Disease, Gene_Match, snpId, GeneID, GeneSymbol, 
                               PathwayID, PathwayName)]
  
  #Now merge the two and determine which has more snp/gene hits per pathway
  high.path.sum$Label<- "High"
  low.path.sum$Label<- "Low"
  
  path.sum<- rbind(high.path.sum, low.path.sum)
  
  path.sum<- path.sum[, .(FoldHit=(TotHits[Label == "High"]/TotHits[Label == "Low"]),
                          HighStates=States[Label == "High"], LowStates=States[Label == "Low"]),
                      by=.(Disease, Gene_Match, snpId, GeneID, GeneSymbol, PathwayID, PathwayName)]
  path.sum<- path.sum[Gene_Match != ""] #not in path
  path.sum<- path.sum[!is.na(HighStates)]
  path.sum<- path.sum[HighStates >= 3]
  
  #Determine how many SNPs in high pathways
  snp.val<- path.sum[, .(TotalSNPs=length(unique(na.omit(snpId)))), by=.(Disease, PathwayID, PathwayName)]
  
  path.sum.plot<- path.sum[, .(SNPs=length(unique(snpId))), by=.(Disease, GeneID, GeneSymbol, 
                                                                 PathwayID, PathwayName, FoldHit)]
  path.sum.plot<- merge(path.sum.plot, snp.val, by=c("Disease", "PathwayID", "PathwayName"))
  
  path.sum.plot
}

############################################################################
###5) Pesticide priority list

pest.rank.select<- function(state.set.in, usgs.dat.in, pest.dat.in, snp.dat.in, mba.dat.in){
  
  kg.per.pest<- pest.dat.in[, .(KG=sum(KG_APPLIED), Years=.N), by=.(CASRN, COMPOUND, State, Area)]
  
  pest.rank.dat<- usgs.dat.in[, .(.N), by=.(CASRN, ChemicalName, GeneID, GeneSymbol, snpId, Disease=UMLS_Name)]
  pest.rank.dat<- merge(pest.rank.dat, kg.per.pest[State %in% unique(state.set.in$location_name)], 
                        by="CASRN", allow.cartesian=TRUE)
  
  pest.rank.dat<- merge(pest.rank.dat, state.set.in, by.x=c("Disease", "State"),
                        by.y=c("Disease", "location_name"), allow.cartesian=TRUE)
  
  #Reduce to top SNPs
  pest.rank.dat<- merge(pest.rank.dat, snp.dat.in, by=c("Disease", "snpId", "GeneSymbol", "GeneID"))
  
  #Determine how many SNPs and genes per chemical
  pest.rank.summary<- pest.rank.dat[, .(SNPs=length(unique(na.omit(snpId))), Genes=length(unique(GeneID))), 
                                    by=.(Disease, State, Area, KG, Years,
                                         CASRN, ChemicalName, Label)]
  #Determine KG applied 
  pest.rank.summary$KGArea<- pest.rank.summary$KG/pest.rank.summary$Area
  
  #Per state, label, determine sum and mean of values
  pest.rank.high<- pest.rank.summary[, .(Years=sum(Years), KG=sum(KG), KGArea=sum(KGArea),
                                         meanYears=mean(Years), meanKG=mean(KG), 
                                         meanKGArea=mean(KGArea),
                                         meanGenes=mean(Genes), States=length(unique(State))),
                                     by=.(Disease, Label, CASRN, ChemicalName, SNPs, Genes)]
  
  #Determine when high states have larger values than low states
  pest.rank.high<- pest.rank.high[, .(Years=(Years[Label == 2] - Years[Label == 1]), 
                                      KG=(KG[Label == 2] - KG[Label == 1]), 
                                      KGArea=(KGArea[Label == 2] - KGArea[Label == 1]), 
                                      meanYears=(meanYears[Label == 2] - meanYears[Label == 1]), 
                                      meanKG=(meanKG[Label == 2] - meanKG[Label == 1]), 
                                      meanKGArea=(meanKGArea[Label == 2] - meanKGArea[Label == 1]),
                                      StatesHigh=States[Label == 2], StatesLow=States[Label==1]),
                                  by=.(Disease, CASRN, ChemicalName, SNPs, Genes)]
  pest.rank.high<- pest.rank.high[StatesHigh >= StatesLow]
  
  #Determine which high state pesticides have more KG applied
  pest.rank.fin<- pest.rank.high[KG > 0 | KGArea > 0]
  
  pest.rank.plot<- merge(pest.rank.summary[Label == 2, list(Disease, CASRN, Label, State)], 
                         pest.rank.fin, by=c("Disease","CASRN"))
  pest.rank.plot$Indiv_Pest<- pest.dat.in[match(pest.rank.plot$CASRN, pest.dat.in$CASRN)]$Compound_Simp
  
  #Determine which high state pesticides are in more rules
  #Separate rules into individual pesticides
  mba.int<- mba.dat.in
  mba.int$Indiv_Pest<- mba.int$Pest_Set
  mba.int<- as.data.table(separate_rows(mba.int, Indiv_Pest, sep=","))
  
  mba.alz<- mba.dis.sep(state.set.in, mba.int, "Alzheimer's Disease")
  mba.ms<- mba.dis.sep(state.set.in, mba.int, "Multiple Sclerosis")
  mba.pd<- mba.dis.sep(state.set.in, mba.int, "Parkinson Disease")
  
  rule.dat<- rbind(mba.alz[Label == "Top"], mba.ms[Label == "Top"], mba.pd[Label == "Top"])
  pest.rank.plot<- merge(pest.rank.plot, rule.dat, by=c("State", "Disease", "Indiv_Pest"), allow.cartesian=TRUE)
  
  pest.rank.plot<- pest.rank.plot[, .(State=length(unique(State))), 
                                  by=.(Disease, CASRN, ChemicalName, KG, KGArea, meanYears, meanKG,
                                       SNPs, Genes, Pest_Set, Un_Chem, Un_States)]
  pest.rank.plot<- pest.rank.plot[, .(State=mean(State), Rules=length(unique(Pest_Set))),
                                  by=.(Disease, CASRN, ChemicalName, KG, KGArea, meanYears, meanKG,
                                       SNPs, Genes, Un_Chem, Un_States)]
  
  #Make sure pesticides applied in at least 3 states
  test<- merge(pest.rank.plot, pest.rank.summary[, list(Disease, State, CASRN, Label)],
               by=c("Disease", "CASRN"), allow.cartesian=TRUE)
  test<- test[Label==2, .(StateCount=length(unique(State.y))), by=.(Disease, CASRN)]
  
  pest.rank.plot<- merge(pest.rank.plot, test, by=c("Disease", "CASRN"))
  pest.rank.plot<- pest.rank.plot[StateCount >=3]
  pest.rank.plot<- pest.rank.plot[order(Disease, Un_Chem, -Un_States, -StateCount, -SNPs, -Genes, -KGArea, -Rules)]
  pest.rank.plot[CASRN %in% "120116-88-3"]$ChemicalName<- "Cyazofamid" #Give shorter name
  
  pest.rank.out<- cbind(c(1:length(which(pest.rank.plot$Disease == "Alzheimer's Disease")), 
                          1:length(which(pest.rank.plot$Disease == "Multiple Sclerosis")), 
                          1:length(which(pest.rank.plot$Disease == "Parkinson Disease"))), 
                        pest.rank.plot[, list(Disease, CASRN, ChemicalName, SNPs, Genes, Rules)])
  
  #Add additional details about each chemical
  add.dat<- pest.rank.summary[, .(High_States=length(unique(State[Label == 2])),
                                  Low_States=length(unique(State[Label == 1])),
                                  Avg_YrsHigh=mean(Years[Label==2]), Avg_YrsLow=mean(Years[Label==1]),
                                  Avg_KgArea=mean(KGArea[Label ==2])), by=.(Disease, CASRN)]
  pest.rank.out<- merge(pest.rank.out, add.dat, by=c("Disease", "CASRN"))
  
  #Clean columns for output
  pest.rank.out<- pest.rank.out[, list(Rank=V1, Disease, CASRN, ChemicalName, SNPs, Genes, Rules,
                                       Avg_KgArea, High_States, Low_States, Avg_YrsHigh, Avg_YrsLow)]
  pest.rank.out<- pest.rank.out[order(Disease, Rank)]
  
  pest.rank.out$Avg_KgArea<- format(pest.rank.out$Avg_KgArea, digits=2)
  pest.rank.out$Avg_YrsHigh<- format(pest.rank.out$Avg_YrsHigh, digits=2)
  pest.rank.out$Avg_YrsLow<- format(pest.rank.out$Avg_YrsLow, digits=2)
  
  pest.rank.out
}

############################################################################
###Cycle through repetitions and build priority lists (5 random top vs 5 random bottom)

for(i in c(1:10)){
  
  state.list<- state.select()
  rules.out<- rules.select(state.list, chem.snp.dis.dat, usgs.year.pest, mba.rules)
  print("Rules prepared")
  
  hit.rank.dat<- snp.hit.select(state.list, chem.snp.dis.dat, usgs.year.pest)
  print("SNP hits prepared")
  
  rule.rank.dat<- snp.rule.select(state.list, chem.snp.dis.dat, usgs.year.pest, mba.rules, hit.rank.dat)
  print("SNP rules prepared")
  
  snp.rank.list<- fin.snp.calc(rule.rank.dat, hit.rank.dat)
  print("SNP rank prepared")
  
  path.rank.list<- path.select(state.list, chem.snp.dis.dat, usgs.year.pest)
  print("Path rank prepared")
  
  pest.rank.list<- pest.rank.select(state.list, chem.snp.dis.dat, usgs.year.pest, snp.rank.list, mba.rules)
  print("Pest rank prepared")
  
  assign(paste0("state.list", i), state.list)
  assign(paste0("rules.out", i), rules.out)
  assign(paste0("snp.rank.list", i), snp.rank.list)
  assign(paste0("path.rank.list", i), path.rank.list)
  assign(paste0("pest.rank.list", i), pest.rank.list)
  
  print(paste(i, "repetition completed"))
}

#Clean workspace
val.clear<- c(paste0("state.list", c(1:10)), paste0("rules.out",  c(1:10)), paste0("snp.rank.list", c(1:10)), 
              paste0("path.rank.list", c(1:10)), paste0("pest.rank.list", c(1:10)))

rm(list=setdiff(ls(), c(val.clear, "top.bot.dat", "chem.snp.dis.dat", "mba.rules")))

###SAVE WORKSPACE FOR FUTURE ANALYSIS
setwd("")
save.image(file="Workspaces/Step12a-Post_10State_CrossValidation_Trial.RData")
#load("Workspaces/Step12a-Post_10State_CrossValidation_Trial.RData")

##################################################################################
#Join individual reps together to make on overall table

bind.func<- function(base.name){
  dat.out<- NULL
  for(i in c(1:10)){
    dat.int<- get(paste0(base.name, i))
    dat.int$Rep<- i
    dat.out<- rbind(dat.out, dat.int)
  }
  dat.out
}

full.state.cv<- bind.func("state.list")
full.rule.cv<- bind.func("rules.out")
full.snp.cv<- bind.func("snp.rank.list")
full.path.cv<- bind.func("path.rank.list")
full.pest.cv<- bind.func("pest.rank.list")

#Clean workspace
rm(list=setdiff(ls(), c("full.state.cv", "full.rule.cv", "full.snp.cv", "full.path.cv",
                        "full.pest.cv", "chem.snp.dis.dat", "mba.rules")))

#############################################################################
###Build priority lists
#Based on ranks in individual repetitions, identify which values repeatedly rank highly
#Keep values present in >=6 repetitions, and order based on the number of reps and overall rank

###Priority list of rules
rule.high.cv<- full.rule.cv[, .(Reps=length(unique(Rep))), by=.(Rule, Dataset_Rule, Disease)]
rule.high.cv<- rule.high.cv[Reps >=6]
rule.high.cv<- rule.high.cv[order(Disease, -Reps)]

###Priority list of SNPs
snp.high.cv<- full.snp.cv[, .(Rank=sum(Rank), Reps=length(unique(Rep)), HighHits=mean(HighHits), 
                              Rules=mean(Rules)), by=.(snpId, GeneID, GeneSymbol, Disease)]
snp.high.cv<- snp.high.cv[Reps >=6]
snp.high.cv$RankOrder<- snp.high.cv$Rank/snp.high.cv$Reps
snp.high.cv<- snp.high.cv[order(Disease, Rank=RankOrder)]

###Priority list of pathways
#NOTE: need full path list for figure 4 - reserve
full.path.out<- full.path.cv

#Total fold hit across SNPs and genes
full.path.cv<- full.path.cv[FoldHit >1]
full.path.cv<- full.path.cv[, .(TotalGenes=length(unique(GeneID)), TotalHit=sum(FoldHit)), by=.(Disease, PathwayID, PathwayName, TotalSNPs, Rep)]
full.path.cv<- full.path.cv[order(Rep, Disease, -TotalHit, -TotalSNPs, -TotalGenes)]
full.path.cv$RankValue<- full.path.cv$TotalHit + full.path.cv$TotalSNPs + full.path.cv$TotalGenes
full.path.cv$Value<- 1

#For each rep, determine rank of pathways
path.out<- NULL
for(i in c(1:10)){
  path.int<- full.path.cv[Rep == i]
  path.int[Disease %in% "Alzheimer's Disease"]$Value<- 
    rank(1/(path.int[Disease %in% "Alzheimer's Disease"]$RankValue), ties.method=c("average"))
  path.int[Disease %in% "Multiple Sclerosis"]$Value<- 
    rank(1/(path.int[Disease %in% "Multiple Sclerosis"]$RankValue), ties.method=c("average"))
  path.int[Disease %in% "Parkinson Disease"]$Value<- 
    rank(1/(path.int[Disease %in% "Parkinson Disease"]$RankValue), ties.method=c("average"))
  path.int<- path.int[order(Disease, Value)]
  
  path.int<- path.int[, list(Rank=Value, Disease, PathwayID, PathwayName, HighHits=TotalHit, 
                             SNPs=TotalSNPs, Genes=TotalGenes)]
  path.int$Rep<- i
  path.out<- rbind(path.out, path.int)
  
}

full.path.cv<- path.out
#Now combine ranks across pathways as with other values
path.high.cv<- full.path.cv[, .(Rank=sum(Rank), Reps=length(unique(Rep)), HighHits=mean(HighHits), 
                                SNPs=mean(SNPs), Genes=mean(Genes)), by=.(PathwayID, PathwayName, Disease)]
path.high.cv<- path.high.cv[Reps >=6]
path.high.cv$RankOrder<- path.high.cv$Rank/path.high.cv$Reps
path.high.cv<- path.high.cv[order(Disease, Rank=RankOrder)]


###Priority list of pesticides
full.pest.cv[, c(5:12)]<- lapply(full.pest.cv[, c(5:12)], function(x) as.numeric(as.character(x)))
pest.high.cv<- full.pest.cv[, .(Rank=sum(Rank), Reps=length(unique(Rep)), SNPs=mean(SNPs),
                                Genes=mean(Genes), Rules=mean(Rules), Avg_KgArea=mean(Avg_KgArea),
                                High_States=mean(High_States), Low_States=mean(Low_States),
                                Avg_YrsHigh=mean(Avg_YrsHigh), Avg_YrsLow=mean(Avg_YrsLow)), 
                            by=.(CASRN, ChemicalName, Disease)]
pest.high.cv<- pest.high.cv[Reps >=6]
pest.high.cv$RankOrder<- pest.high.cv$Rank/pest.high.cv$Reps
pest.high.cv<- pest.high.cv[order(Disease, Rank=RankOrder)]

#################################################################################
#Do cross-validation with comparison to the priority lists

#Wilcoxon rank sum test
#Do ranks of individual priority lists match the rank of the overall list?
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


#Hypergeometric test
#Do items in individual priority lists significantly overlap with items in the overall list?
enrich.comp<- function(cv.dat, full.dat, dis.in, background.list){
  
  cv.dat<- cv.dat[Disease %in% dis.in]
  full.dat<- full.dat[Disease %in% dis.in]
  
  phyper(q= length(which(unique(cv.dat$variable) %in% full.dat$variable)),
         m= length(unique(full.dat$variable)),
         n=(length(background.list) - length(unique(full.dat$variable))),
         k= length(unique(cv.dat$variable)), lower.tail=FALSE)
  
}

cv.dat.high<- NULL

for(i in c(1:10)){
  
  #Get variables for a given cross validation rep
  cv.rep.states<- full.state.cv[Rep == i]
  cv.rep.rules<- full.rule.cv[Rep == i]
  cv.rep.snps<- full.snp.cv[Rep == i]
  cv.rep.paths<- full.path.cv[Rep == i]
  cv.rep.pests<- full.pest.cv[Rep == i]
  
  dat.int<- NULL
  
  #Get full list of prioritized items
  full.list.rules<- rule.high.cv
  full.list.snps<- snp.high.cv
  full.list.paths<- path.high.cv
  full.list.pests<- pest.high.cv
  
  #Rank compare SNPs
  alz.snp.w<- wilcox.comp(cv.rep.snps[!is.na(snpId), .(.N), by=.(variable=snpId, Rank, Disease)], 
                          full.list.snps[!is.na(snpId), .(.N), by=.(variable=snpId, Rank, Disease)], "Alzheimer's Disease")
  ms.snp.w<- wilcox.comp(cv.rep.snps[!is.na(snpId), .(.N), by=.(variable=snpId, Rank, Disease)], 
                         full.list.snps[!is.na(snpId), .(.N), by=.(variable=snpId, Rank, Disease)], "Multiple Sclerosis")
  pd.snp.w<- wilcox.comp(cv.rep.snps[!is.na(snpId), .(.N), by=.(variable=snpId, Rank, Disease)], 
                         full.list.snps[!is.na(snpId), .(.N), by=.(variable=snpId, Rank, Disease)], "Parkinson Disease")
  
  dat.fill<- as.data.table(cbind(rep("SNP_Wilcox", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.snp.w$p.value, ms.snp.w$p.value, pd.snp.w$p.value)))
  
  dat.int<- rbind(dat.int, dat.fill)
  
  #Rank compare paths
  alz.path<- wilcox.comp(cv.rep.paths[, list(variable=PathwayID, Rank, Disease)], 
                         full.list.paths[, list(variable=PathwayID, Rank, Disease)], "Alzheimer's Disease")
  ms.path<- wilcox.comp(cv.rep.paths[, list(variable=PathwayID, Rank, Disease)], 
                        full.list.paths[, list(variable=PathwayID, Rank, Disease)], "Multiple Sclerosis")
  pd.path<- wilcox.comp(cv.rep.paths[, list(variable=PathwayID, Rank, Disease)], 
                        full.list.paths[, list(variable=PathwayID, Rank, Disease)], "Parkinson Disease")
  
  dat.fill<- as.data.table(cbind(rep("Path_Wilcox", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.path$p.value, ms.path$p.value, pd.path$p.value)))
  
  dat.int<- rbind(dat.int, dat.fill)
  
  #Rank compare pesticides
  alz.pest<- wilcox.comp(cv.rep.pests[, list(variable=CASRN, Rank, Disease)], 
                         full.list.pests[, list(variable=CASRN, Rank, Disease)], "Alzheimer's Disease")
  ms.pest<- wilcox.comp(cv.rep.pests[, list(variable=CASRN, Rank, Disease)], 
                        full.list.pests[, list(variable=CASRN, Rank, Disease)], "Multiple Sclerosis")
  pd.pest<- wilcox.comp(cv.rep.pests[, list(variable=CASRN, Rank, Disease)], 
                        full.list.pests[, list(variable=CASRN, Rank, Disease)], "Parkinson Disease")
  
  dat.fill<- as.data.table(cbind(rep("Pest_Wilcox", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.pest$p.value, ms.pest$p.value, pd.pest$p.value)))
  
  dat.int<- rbind(dat.int, dat.fill)
  
  #Enrichment test rules
  alz.rules.e<- enrich.comp(cv.rep.rules[, list(variable=Rule, Disease)],
                            full.list.rules[, list(variable=Rule, Disease)],
                            "Alzheimer's Disease", background.list=unique(mba.rules$Pest_Set))
  ms.rules.e<- enrich.comp(cv.rep.rules[Rule!="", list(variable=Rule, Disease)],
                           full.list.rules[!is.na(Rule), list(variable=Rule, Disease)],
                           "Multiple Sclerosis", background.list=unique(mba.rules$Pest_Set))
  pd.rules.e<- enrich.comp(cv.rep.rules[Rule!="", list(variable=Rule, Disease)],
                           full.list.rules[!is.na(Rule), list(variable=Rule, Disease)],
                           "Parkinson Disease", background.list=unique(mba.rules$Pest_Set))

  dat.fill<- as.data.table(cbind(rep("Rules_Enrich", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.rules.e, ms.rules.e, pd.rules.e)))

  dat.int<- rbind(dat.int, dat.fill)
  
  #Enrichment test SNPs
  alz.snp.e<- enrich.comp(cv.rep.snps[!is.na(snpId), list(variable=snpId, Disease)], 
                          full.list.snps[!is.na(snpId), list(variable=snpId, Disease)], 
                          "Alzheimer's Disease", background.list=unique(na.omit(chem.snp.dis.dat$snpId)))
  ms.snp.e<- enrich.comp(cv.rep.snps[snpId!="", list(variable=snpId, Disease)], 
                         full.list.snps[!is.na(snpId), list(variable=snpId, Disease)], 
                         "Multiple Sclerosis", background.list=unique(na.omit(chem.snp.dis.dat$snpId)))
  pd.snp.e<- enrich.comp(cv.rep.snps[snpId!="", list(variable=snpId, Disease)], 
                         full.list.snps[!is.na(snpId), list(variable=snpId, Disease)], 
                         "Parkinson Disease", background.list=unique(na.omit(chem.snp.dis.dat$snpId)))
  
  dat.fill<- as.data.table(cbind(rep("SNP_Enrich", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.snp.e, ms.snp.e, pd.snp.e)))
  
  dat.int<- rbind(dat.int, dat.fill)
  
  #Enrichment test pathways
  alz.path.e<- enrich.comp(cv.rep.paths[, list(variable=PathwayID, Disease)], 
                           full.list.paths[, list(variable=PathwayID, Disease)], 
                           "Alzheimer's Disease", background.list=unique(na.omit(chem.snp.dis.dat$PathwayID)))
  ms.path.e<- enrich.comp(cv.rep.paths[, list(variable=PathwayID, Disease)], 
                          full.list.paths[, list(variable=PathwayID, Disease)], 
                          "Multiple Sclerosis", background.list=unique(na.omit(chem.snp.dis.dat$PathwayID)))
  pd.path.e<- enrich.comp(cv.rep.paths[, list(variable=PathwayID, Disease)], 
                          full.list.paths[, list(variable=PathwayID, Disease)], 
                          "Parkinson Disease", background.list=unique(na.omit(chem.snp.dis.dat$PathwayID)))
  
  dat.fill<- as.data.table(cbind(rep("Path_Enrich", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.path.e, ms.path.e, pd.path.e)))
  
  dat.int<- rbind(dat.int, dat.fill)
  
  #Enrichment test pesticides
  alz.pest.e<- enrich.comp(cv.rep.pests[, list(variable=CASRN, Disease)], 
                           full.list.pests[, list(variable=CASRN, Disease)], 
                           "Alzheimer's Disease", background.list=unique(na.omit(chem.snp.dis.dat$CASRN)))
  ms.pest.e<- enrich.comp(cv.rep.pests[, list(variable=CASRN, Disease)], 
                          full.list.pests[, list(variable=CASRN, Disease)], 
                          "Multiple Sclerosis", background.list=unique(na.omit(chem.snp.dis.dat$CASRN)))
  pd.pest.e<- enrich.comp(cv.rep.pests[, list(variable=CASRN, Disease)], 
                          full.list.pests[, list(variable=CASRN, Disease)], 
                          "Parkinson Disease", background.list=unique(na.omit(chem.snp.dis.dat$CASRN)))
  
  dat.fill<- as.data.table(cbind(rep("Pest_Enrich", 3),
                                 c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"),
                                 c(alz.pest.e, ms.pest.e, pd.pest.e)))
  
  dat.int<- rbind(dat.int, dat.fill)
  
  dat.int$Rep<- i
  colnames(dat.int)[1:3]<- c("Method", "Disease", "Value")
  dat.int$Value<- as.numeric(as.character(dat.int$Value))
  
  dat.int$State_High<- "X"
  dat.int$State_Low<- "X"
  dat.int[Disease == "Alzheimer's Disease"]$State_High<- paste0(cv.rep.states[Disease == "Alzheimer's Disease" & Label ==2]$location_name,
                                                                collapse=";")
  dat.int[Disease == "Multiple Sclerosis"]$State_High<- paste0(cv.rep.states[Disease == "Multiple Sclerosis" & Label ==2]$location_name,
                                                               collapse=";")
  dat.int[Disease == "Parkinson Disease"]$State_High<- paste0(cv.rep.states[Disease == "Parkinson Disease" & Label ==2]$location_name,
                                                              collapse=";")

  dat.int[Disease == "Alzheimer's Disease"]$State_Low<- paste0(cv.rep.states[Disease == "Alzheimer's Disease" & Label ==1]$location_name,
                                                               collapse=";")
  dat.int[Disease == "Multiple Sclerosis"]$State_Low<- paste0(cv.rep.states[Disease == "Multiple Sclerosis" & Label ==1]$location_name,
                                                              collapse=";")
  dat.int[Disease == "Parkinson Disease"]$State_Low<- paste0(cv.rep.states[Disease == "Parkinson Disease" & Label ==1]$location_name,
                                                             collapse=";")
  
  cv.dat.high<- as.data.table(rbind(cv.dat.high, dat.int))
}

#Prepare SI Table S2
cv.dat.out<- cv.dat.high[!Method %in% "Rules_Enrich"]
cv.dat.out[Value < 1e-16]$Value<- 0
cv.dat.out$Value<- format(cv.dat.out$Value, digits=2)
cv.dat.out[Value == "0.0e+00"]$Value<- " <1e-16"
cv.dat.out<- dcast(cv.dat.out, Rep+Disease~Method, value.var="Value")

setwd("")
fwrite(cv.dat.out, file="SITables_Figures/SITable2_CV_10State_Stats.csv")

#Significant Wilcox p-value means the list ranks are different
#Significant enrichment p-value means the list items overlap

#6, 1, 6
length(which(cv.dat.high[Method == "SNP_Wilcox" & Disease == "Alzheimer's Disease"]$Value >=0.05)) 
length(which(cv.dat.high[Method == "SNP_Wilcox" & Disease == "Multiple Sclerosis"]$Value >=0.05)) 
length(which(cv.dat.high[Method == "SNP_Wilcox" & Disease == "Parkinson Disease"]$Value >=0.05)) 

#10, 9, 10
length(which(cv.dat.high[Method == "Path_Wilcox" & Disease == "Alzheimer's Disease"]$Value >=0.05)) 
length(which(cv.dat.high[Method == "Path_Wilcox" & Disease == "Multiple Sclerosis"]$Value >=0.05)) 
length(which(cv.dat.high[Method == "Path_Wilcox" & Disease == "Parkinson Disease"]$Value >=0.05)) 

#8, 7, 8
length(which(cv.dat.high[Method == "Pest_Wilcox" & Disease == "Alzheimer's Disease"]$Value >=0.05)) 
length(which(cv.dat.high[Method == "Pest_Wilcox" & Disease == "Multiple Sclerosis"]$Value >=0.05)) 
length(which(cv.dat.high[Method == "Pest_Wilcox" & Disease == "Parkinson Disease"]$Value >=0.05)) 

#all sig overlap
length(which(cv.dat.high[Method == "Rules_Enrich" & Disease == "Alzheimer's Disease"]$Value < 0.05))
length(which(cv.dat.high[Method == "Rules_Enrich" & Disease == "Multiple Sclerosis"]$Value < 0.05))
length(which(cv.dat.high[Method == "Rules_Enrich" & Disease == "Parkinson Disease"]$Value < 0.05))

#all sig overlap
length(which(cv.dat.high[Method == "SNP_Enrich" & Disease == "Alzheimer's Disease"]$Value < 0.05))
length(which(cv.dat.high[Method == "SNP_Enrich" & Disease == "Multiple Sclerosis"]$Value < 0.05))
length(which(cv.dat.high[Method == "SNP_Enrich" & Disease == "Parkinson Disease"]$Value < 0.05))

#all sig overlap
length(which(cv.dat.high[Method == "Path_Enrich" & Disease == "Alzheimer's Disease"]$Value < 0.05))
length(which(cv.dat.high[Method == "Path_Enrich" & Disease == "Multiple Sclerosis"]$Value < 0.05))
length(which(cv.dat.high[Method == "Path_Enrich" & Disease == "Parkinson Disease"]$Value < 0.05))

#all sig overlap
length(which(cv.dat.high[Method == "Pest_Enrich" & Disease == "Alzheimer's Disease"]$Value < 0.05))
length(which(cv.dat.high[Method == "Pest_Enrich" & Disease == "Multiple Sclerosis"]$Value < 0.05))
length(which(cv.dat.high[Method == "Pest_Enrich" & Disease == "Parkinson Disease"]$Value < 0.05))

#############################################################################
#Prepare workspace of priority lists
#Priority rules will be those in high states, not subsets

rm(list=setdiff(ls(), c("snp.high.cv", "path.high.cv", "pest.high.cv", "full.path.out")))
#Need full path list for figure 4

###NOTE: SAVE WORKSPACE FOR FUTURE ANALYSIS
#Copy of workspace available currently as "PriorityList_Workspace.RData" used in subsequent steps
setwd("")
save.image("Workspaces/FILENAME.RData")
