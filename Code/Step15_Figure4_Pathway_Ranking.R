############################################################################

### Pathway prioritization and fold hit calculation
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library (plyr)
library(readxl)

###NOTE: LOAD OWN WORKSPACE FROM STEP 12a, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/PriorityList_Workspace.RData")

#Only need priority pathways
rm(list=setdiff(ls(), c("path.high.cv", "full.path.out")))

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/Analysis_Workspace.RData")
chem.snp.dis.dat<- chem.snp.dis.dat[UMLS_Name %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]

###NOTE: LOAD OWN DATA FROM STEP 10, OR WORK FROM EXISTING DATASET (BELOW)
top.bot.dat<- fread("Top_Bottom_States_MeanSNPs.csv")

###Combine pathway hits across reps - take average

full.path.plot<- full.path.out[, .(FoldHit=mean(FoldHit)), by=.(Disease, PathwayID, PathwayName, GeneID, GeneSymbol, SNPs, TotalSNPs)]

#How many high/low hits?
test.low<- full.path.plot[FoldHit < 1]
test.high<- full.path.plot[FoldHit > 1]

length(unique(test.high$PathwayID))
length(which(!unique(test.high$PathwayID) %in% test.low$PathwayID))


############################################################################
#Plot of fold hits in gene-pathway

fold.plot<- full.path.plot
fold.plot$Disease<- factor(fold.plot$Disease, levels=c("Parkinson Disease", "Multiple Sclerosis", "Alzheimer's Disease"), ordered=TRUE)
range(fold.plot$FoldHit)

setwd("")
jpeg("Figure4c.jpg", width=2.5, height=1.5, units="in", res=300)
ggplot(fold.plot) + geom_density(aes(x=FoldHit, fill=Disease), alpha=0.7) + #facet_wrap(~Disease) + 
  scale_fill_manual(values=c("Alzheimer's Disease"="mediumpurple3", "Multiple Sclerosis"="goldenrod2", 
                             "Parkinson Disease"="darkcyan"), guide=NULL) + xlab("Fold difference high vs low pathway hits")+
  geom_vline(xintercept = 1, color="grey25", linewidth=1) +
  theme_bw() + scale_x_continuous(expand=c(0, 0), trans="sqrt") + scale_y_continuous(expand=c(0,0), limits=c(0, 5)) +
  theme(panel.grid = element_blank(), text=element_text(size=8))
graphics.off()

setwd("")
fwrite(fold.plot, "SITables_Figures/Figure4c.csv")

############################################################################
#Assess high pathways

path.high.cv

#92, 84, 89
length(unique(path.high.cv[Disease %in% "Alzheimer's Disease"]$PathwayID))/length(unique(chem.snp.dis.dat[UMLS_Name %in% "Alzheimer's Disease"]$PathwayID))
length(unique(path.high.cv[Disease %in% "Multiple Sclerosis"]$PathwayID))/length(unique(chem.snp.dis.dat[UMLS_Name %in% "Multiple Sclerosis"]$PathwayID))
length(unique(path.high.cv[Disease %in% "Parkinson Disease"]$PathwayID))/length(unique(chem.snp.dis.dat[UMLS_Name %in% "Parkinson Disease"]$PathwayID))

#Same top pathways as original approach
path.high.cv[Disease %in% "Alzheimer's Disease"]
path.high.cv[Disease %in% "Multiple Sclerosis"]
path.high.cv[Disease %in% "Parkinson Disease"]

#Top pathways and genes table
high.path.genes<- fold.plot[FoldHit > 1, .(TotalFold=sum(FoldHit), Paths=length(unique(PathwayID))), 
                            by=.(Disease, GeneID, GeneSymbol, SNPs)]
high.path.genes<- high.path.genes[order(Disease, -TotalFold)]
top.gen<- rbind(head(high.path.genes[Disease %in% "Alzheimer's Disease"], 5),
                head(high.path.genes[Disease %in% "Multiple Sclerosis"], 5),
                head(high.path.genes[Disease %in% "Parkinson Disease"], 5))
top.path<- rbind(head(path.high.cv[Disease %in% "Alzheimer's Disease"], 4), #Has replicate names
                 head(path.high.cv[Disease %in% "Multiple Sclerosis"], 3),
                 head(path.high.cv[Disease %in% "Parkinson Disease"], 3))

#Clean output for table

data.out<- as.data.table(rbind("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease"))
data.out<- cbind(data.out, rbind(paste(unique(top.path[Disease %in% "Alzheimer's Disease"]$PathwayName), collapse=", "),
                                 paste(unique(top.path[Disease %in% "Multiple Sclerosis"]$PathwayName), collapse=", "),
                                 paste(unique(top.path[Disease %in% "Parkinson Disease"]$PathwayName), collapse=", ")),
                 rbind(paste(unique(top.gen[Disease %in% "Alzheimer's Disease"]$GeneSymbol), collapse=", "),
                       paste(unique(top.gen[Disease %in% "Multiple Sclerosis"]$GeneSymbol), collapse=", "),
                       paste(unique(top.gen[Disease %in% "Parkinson Disease"]$GeneSymbol), collapse=", ")))


setwd("")
fwrite(data.out, "Figure4B_Pathway_Table.csv")

############################################################################
#Ranked pathways table

path.dat.out<- path.high.cv
path.dat.out$Value<- 1
path.dat.out[Disease %in% "Alzheimer's Disease"]$Value<- 
  rank(path.dat.out[Disease %in% "Alzheimer's Disease"]$RankOrder, ties.method=c("average"))
path.dat.out[Disease %in% "Multiple Sclerosis"]$Value<- 
  rank(path.dat.out[Disease %in% "Multiple Sclerosis"]$RankOrder, ties.method=c("average"))
path.dat.out[Disease %in% "Parkinson Disease"]$Value<- 
  rank(path.dat.out[Disease %in% "Parkinson Disease"]$RankOrder, ties.method=c("average"))
path.dat.out<- path.dat.out[order(Disease, Value)]

path.dat.out<- path.dat.out[, list(Rank=Value, Disease, PathwayID, PathwayName, CV_Reps=Reps, 
                                 Avg_CV_Rank=Rank/Reps, Avg_CV_HighHits=HighHits, Avg_CV_SNPs=SNPs,
                                 Avg_CV_Genes=Genes)]

setwd("")
fwrite(path.dat.out, "SITables_Figures/Dataset5_Pathway_Ranks.csv")
