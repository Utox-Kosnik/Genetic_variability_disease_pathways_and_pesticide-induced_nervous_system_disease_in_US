############################################################################

### Pesticide prioritization, bubble plots, and application over time
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

###NOTE: LOAD OWN WORKSPACE FROM STEP 12a, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/PriorityList_Workspace.RData")

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
###Prep bubble plots of number of SNPs, genes, and rules

pest.rank.plot$Disease<- factor(pest.rank.plot$Disease, levels=unique(pest.rank.plot$Disease), ordered=TRUE)

quantile(pest.rank.plot$SNPs)
quantile(pest.rank.plot$Genes)
quantile(pest.rank.plot$Rules)

gene.snp.plot<- ggplot(pest.rank.plot, aes(x=Disease, y= reorder(ChemicalName, desc(ChemicalName)))) +  
  geom_point(aes(fill=SNPs, size=Genes), shape=21, color="grey50") + 
  scale_fill_gradient(low="white", high="midnightblue", trans="sqrt", breaks=c(0, 10, 50, 150, 300),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_size_continuous(breaks=c(5, 10, 100, 500), range = c(.5, 3)) + 
  scale_x_discrete(position="top", labels=function(x) str_wrap(x, width=10)) +
  ylab(NULL) + xlab(NULL) +
  theme_bw() + 
  theme(panel.grid=element_blank(), #panel.border=element_rect(color=NA, fill=NA), #panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour=c("mediumpurple3", "goldenrod2", "darkcyan"), face="bold"))

rules.plot<- ggplot(pest.rank.plot, aes(x=Disease, y= reorder(ChemicalName, desc(ChemicalName)))) +  
  geom_point(aes(fill=Reps, size=Rules), shape=21, color="grey50") + 
  scale_fill_gradient(low="white", high="steelblue", 
                      guide=guide_colourbar(title=str_wrap("Number repetitions", width=10), frame.colour="black", ticks.colour = "black")) + 
  scale_size_continuous(breaks=c(500, 1000, 5000, 1000), limits=c(200, 10005), range = c(.5, 3)) + 
  scale_x_discrete(position="top", labels=function(x) str_wrap(x, width=10)) +
  ylab(NULL) + xlab(NULL) +
  theme_bw() + 
  theme(panel.grid=element_blank(), #panel.border=element_rect(color=NA, fill=NA), #panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour=c("mediumpurple3", "goldenrod2", "darkcyan"), face="bold"))

setwd("")
fwrite(pest.rank.plot, "SITables_Figures/Figure5ab.csv")

############################################################################
#Plot full set of application data over time per disease 

usgs.year.pest$KG_Area<- usgs.year.pest$KG_APPLIED/usgs.year.pest$Area

#Extract kg/area data for Alzheimer's priority compounds
alz.compounds<- usgs.year.pest[State %in% top.bot.dat[Label != 0 & Disease %in% "Alzheimer's Disease"]$location_name &
                            CASRN %in% pest.rank.plot[Disease %in% "Alzheimer's Disease"]$CASRN]
alz.compounds<- dcast(alz.compounds, CASRN+YEAR~State, value.var="KG_Area")
alz.compounds[is.na(alz.compounds)]<- 0
alz.compounds$ChemicalName<- pest.rank.plot[match(alz.compounds$CASRN, pest.rank.plot$CASRN)]$ChemicalName
alz.compounds$ChemicalName<- factor(alz.compounds$ChemicalName, 
                                    levels=unique(pest.rank.plot[Disease %in% "Alzheimer's Disease"]$ChemicalName), 
                                    ordered=TRUE)

#All priority Alzheimer's Disease compounds
ggplot(data=alz.compounds, aes(x=YEAR)) + geom_line(aes(y=`Rhode Island`), color="red4") +
  geom_line(aes(y=Maine), color="red3") + geom_line(aes(y=`New Hampshire`), color="red2") +
  geom_line(aes(y=Connecticut), color="red1") + geom_line(aes(y=Massachusetts), color="red") +
  geom_line(aes(y=Vermont), color="tomato") + geom_line(aes(y=Florida), color="indianred3") +
  geom_line(aes(y=Delaware), color="indianred2") + geom_line(aes(y=`New York`), color="indianred1") +
  geom_line(aes(y=`West Virginia`), color="indianred") + geom_line(aes(y=Missouri), color="steelblue") +
  #geom_line(aes(y=California), color="steelblue1") +
  geom_line(aes(y=Nevada), color="steelblue2") +
  geom_line(aes(y=Georgia), color="steelblue3") + geom_line(aes(y=Oklahoma), color="royalblue4") +
  geom_line(aes(y=Texas), color="royalblue3") + geom_line(aes(y=Mississippi), color="royalblue2") +
  geom_line(aes(y=Arkansas), color="royalblue1") + geom_line(aes(y=Louisiana), color="royalblue") +
  geom_line(aes(y=Utah), color="midnightblue") + facet_wrap(.~ChemicalName, scales="free_y") + theme_bw() +
  ylab("Alzheimer's Disease") + xlab(NULL) + scale_y_continuous(trans="sqrt") +
  ggtitle("kg pesticide applied/square mile") +
  scale_x_continuous(limits=c(1992, 2018), breaks = c(1992, 2000, 2010, 2018)) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=8, hjust=0.5), 
        text=element_text(size=8, family="sans"), axis.title.y=element_text(colour="mediumpurple3", face="bold"))

#Reduced compound for figure 5
alz.red.compounds<- alz.compounds[ChemicalName %in% c("Mancozeb", "Fenarimol", "Phosmet", "Captan")]

alz.red.plots<- ggplot(data=alz.red.compounds, aes(x=YEAR)) + geom_line(aes(y=`Rhode Island`), color="red4") +
  geom_line(aes(y=Maine), color="red3") + geom_line(aes(y=`New Hampshire`), color="red2") +
  geom_line(aes(y=Connecticut), color="red1") + geom_line(aes(y=Massachusetts), color="red") +
  geom_line(aes(y=Texas), color="royalblue3") + geom_line(aes(y=Mississippi), color="royalblue2") +
  geom_line(aes(y=Arkansas), color="royalblue1") + geom_line(aes(y=Louisiana), color="royalblue") +
  geom_line(aes(y=Utah), color="midnightblue") + facet_wrap(.~ChemicalName, scales="free_y", nrow=1) + theme_bw() +
  ylab("Alzheimer's Disease") + xlab(NULL) + scale_y_continuous(trans="sqrt") +
  ggtitle("kg pesticide applied/square mile") +
  scale_x_continuous(limits=c(1992, 2018), breaks = c(1992, 2000, 2010, 2018)) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=8, hjust=0.5), 
        text=element_text(size=8, family="sans"), axis.title.y=element_text(colour="mediumpurple3", face="bold"))

#All MS priority compounds
ms.compounds<- usgs.year.pest[State %in% top.bot.dat[Label != 0 & Disease %in% "Multiple Sclerosis"]$location_name &
                                 CASRN %in% pest.rank.plot[Disease %in% "Multiple Sclerosis"]$CASRN]
ms.compounds<- dcast(ms.compounds, CASRN+YEAR~State, value.var="KG_Area")
ms.compounds[is.na(ms.compounds)]<- 0
ms.compounds$ChemicalName<- pest.rank.plot[match(ms.compounds$CASRN, pest.rank.plot$CASRN)]$ChemicalName
ms.compounds$ChemicalName<- factor(ms.compounds$ChemicalName, 
                                    levels=unique(pest.rank.plot[Disease %in% "Multiple Sclerosis"]$ChemicalName), 
                                    ordered=TRUE)

ggplot(data=ms.compounds, aes(x=YEAR)) + geom_line(aes(y=`Rhode Island`), color="red4") +
  geom_line(aes(y=`New Hampshire`), color="red3") + geom_line(aes(y=Vermont), color="red2") +
  geom_line(aes(y=Massachusetts), color="red1") + geom_line(aes(y=Connecticut), color="red") +
  geom_line(aes(y=Maine), color="tomato") + geom_line(aes(y=Utah), color="indianred3") +
  geom_line(aes(y=Idaho), color="indianred2") + geom_line(aes(y=`New York`), color="indianred1") +
  geom_line(aes(y=Wyoming), color="indianred") + geom_line(aes(y=`South Carolina`), color="steelblue") +
  #geom_line(aes(y=California), color="steelblue1") + 
  geom_line(aes(y=Georgia), color="steelblue2") +
  geom_line(aes(y=Kansas), color="steelblue3") + geom_line(aes(y=Alabama), color="royalblue4") +
  geom_line(aes(y=Texas), color="royalblue3") + geom_line(aes(y=Oklahoma), color="royalblue2") +
  geom_line(aes(y=Arkansas), color="royalblue1") + geom_line(aes(y=Mississippi), color="royalblue") +
  geom_line(aes(y=Louisiana), color="midnightblue") + facet_wrap(.~ChemicalName, scales="free_y") + theme_bw() +
  ylab("Multiple Sclerosis") + xlab(NULL) + scale_y_continuous(trans="sqrt") +
  scale_x_continuous(limits=c(1992, 2018), breaks = c(1992, 2000, 2010, 2018)) +
  theme(panel.grid=element_blank(), #panel.border=element_rect(color=NA, fill=NA), #panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.title.y=element_text(colour="goldenrod2", face="bold"))

#MS compounds for figure 5
ms.red.compounds<- ms.compounds[ChemicalName %in% c("Zoxamide","Polymarcine", "Dimethomorph", "Triflumizole")]

ms.red.plots<- ggplot(data=ms.red.compounds, aes(x=YEAR)) + geom_line(aes(y=`Rhode Island`), color="red4") +
  geom_line(aes(y=`New Hampshire`), color="red3") + geom_line(aes(y=Vermont), color="red2") +
  geom_line(aes(y=Massachusetts), color="red1") + geom_line(aes(y=Connecticut), color="red") +
  geom_line(aes(y=Texas), color="royalblue3") + geom_line(aes(y=Oklahoma), color="royalblue2") +
  geom_line(aes(y=Arkansas), color="royalblue1") + geom_line(aes(y=Mississippi), color="royalblue") +
  geom_line(aes(y=Louisiana), color="midnightblue") + facet_wrap(.~ChemicalName, scales="free_y", nrow=1) + theme_bw() +
  ylab("Multiple Sclerosis") + xlab(NULL) + scale_y_continuous(trans="sqrt") +
  scale_x_continuous(limits=c(1992, 2018), breaks = c(1992, 2000, 2010, 2018)) +
  theme(panel.grid=element_blank(), #panel.border=element_rect(color=NA, fill=NA), #panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.title.y=element_text(colour="goldenrod2", face="bold"))

#PD priority compounds
pd.compounds<- usgs.year.pest[State %in% top.bot.dat[Label != 0 & Disease %in% "Parkinson Disease"]$location_name &
                                CASRN %in% pest.rank.plot[Disease %in% "Parkinson Disease"]$CASRN]
pd.compounds<- dcast(pd.compounds, CASRN+YEAR~State, value.var="KG_Area")
pd.compounds[is.na(pd.compounds)]<- 0
pd.compounds$ChemicalName<- pest.rank.plot[match(pd.compounds$CASRN, pest.rank.plot$CASRN)]$ChemicalName
pd.compounds$ChemicalName<- factor(pd.compounds$ChemicalName, 
                                   levels=unique(pest.rank.plot[Disease %in% "Parkinson Disease"]$ChemicalName), 
                                   ordered=TRUE)


#All PD compounds
ggplot(data=pd.compounds, aes(x=YEAR)) + geom_line(aes(y=Maine), color="red4") +
  geom_line(aes(y=Vermont), color="red3") + geom_line(aes(y=`New Hampshire`), color="red2") +
  geom_line(aes(y=Massachusetts), color="red1") + geom_line(aes(y=`Rhode Island`), color="red") +
  geom_line(aes(y=Connecticut), color="tomato") + geom_line(aes(y=`New York`), color="indianred3") +
  geom_line(aes(y=Florida), color="indianred2") + geom_line(aes(y=Arizona), color="indianred1") +
  geom_line(aes(y=`West Virginia`), color="indianred") + geom_line(aes(y=Nebraska), color="steelblue") +
  geom_line(aes(y=Virginia), color="steelblue1") + geom_line(aes(y=`South Dakota`), color="steelblue2") +
  geom_line(aes(y=Georgia), color="steelblue3") + geom_line(aes(y=Texas), color="royalblue4") +
  geom_line(aes(y=Oklahoma), color="royalblue3") + geom_line(aes(y=Mississippi), color="royalblue2") +
  geom_line(aes(y=Utah), color="royalblue1") + geom_line(aes(y=Louisiana), color="royalblue") +
  geom_line(aes(y=Arkansas), color="midnightblue") + facet_wrap(.~ChemicalName, scales="free_y") + theme_bw() +
  ylab("Parkinson Disease") + xlab("Year pesticide applied") + scale_y_continuous(trans="sqrt") +
  scale_x_continuous(limits=c(1992, 2018), breaks = c(1992, 2000, 2010, 2018)) +
  theme(panel.grid=element_blank(), #panel.border=element_rect(color=NA, fill=NA), #panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.title.y=element_text(colour="darkcyan", face="bold"))

#Compounds for figure 5
pd.red.compounds<- pd.compounds[ChemicalName %in% c("Cyprodinil", 
                                                    "Copper","Maneb","Ziram")]

pd.red.plots<- ggplot(data=pd.red.compounds, aes(x=YEAR)) + geom_line(aes(y=Maine), color="red4") +
  geom_line(aes(y=Vermont), color="red3") + geom_line(aes(y=`New Hampshire`), color="red2") +
  geom_line(aes(y=Massachusetts), color="red1") + geom_line(aes(y=`Rhode Island`), color="red") +
  geom_line(aes(y=Oklahoma), color="royalblue3") + geom_line(aes(y=Mississippi), color="royalblue2") +
  geom_line(aes(y=Utah), color="royalblue1") + geom_line(aes(y=Louisiana), color="royalblue") +
  geom_line(aes(y=Arkansas), color="midnightblue") + facet_wrap(.~ChemicalName, scales="free_y", nrow=1) + theme_bw() +
  ylab("Parkinson Disease") + xlab("Year pesticide applied") + scale_y_continuous(trans="sqrt") +
  scale_x_continuous(limits=c(1992, 2018), breaks = c(1992, 2000, 2010, 2018)) +
  theme(panel.grid=element_blank(), #panel.border=element_rect(color=NA, fill=NA), #panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.title.y=element_text(colour="darkcyan", face="bold"))



setwd("C:/Users/markos/Documents/Manuscripts/23-CSD_Manuscript/EHP_Submission/Revision/Code")
jpeg("Figure5.jpeg", width=7, height=8.2, units="in", res=150)
grid.arrange(gene.snp.plot, rules.plot, alz.red.plots, ms.red.plots, pd.red.plots, nrow=4,
             heights=c(3.15, 1.1, 0.95, 1.05), layout_matrix=rbind(c(1,2), c(3,3), c(4,4), c(5,5)))
graphics.off()

###Write data underlying trendline figures
alz.red.compounds$Disease<- "Alzheimer's Disease"
ms.red.compounds$Disease<- "Multiple Sclerosis"
pd.red.compounds$Disease<- "Parkinson Disease"

alz.red.compounds<- melt(alz.red.compounds, id.vars=c("Disease", "YEAR", "CASRN", "ChemicalName"))
ms.red.compounds<- melt(ms.red.compounds, id.vars=c("Disease", "YEAR", "CASRN", "ChemicalName"))
pd.red.compounds<- melt(pd.red.compounds, id.vars=c("Disease", "YEAR", "CASRN", "ChemicalName"))
fig5.dat<- rbind(alz.red.compounds, ms.red.compounds, pd.red.compounds)

setwd("")
fwrite(fig5.dat, "SITables_Figures/Figure5c.csv")

############################################################################
###Prepare output file

#Clean columns for output
pest.rank.out<- merge(pest.rank.plot, pest.high.cv[, list(CASRN, Disease, Avg_CV_Rank=Rank/Reps)])
pest.rank.out<- pest.rank.out[order(Disease, RankOrder)]
pest.rank.out$Value<- 1
pest.rank.out[Disease %in% "Alzheimer's Disease"]$Value<- 
  rank(pest.rank.out[Disease %in% "Alzheimer's Disease"]$RankOrder, ties.method=c("average"))
pest.rank.out[Disease %in% "Multiple Sclerosis"]$Value<- 
  rank(pest.rank.out[Disease %in% "Multiple Sclerosis"]$RankOrder, ties.method=c("average"))
pest.rank.out[Disease %in% "Parkinson Disease"]$Value<- 
  rank(pest.rank.out[Disease %in% "Parkinson Disease"]$RankOrder, ties.method=c("average"))

pest.rank.out<- pest.rank.out[, list(Disease, Rank=Value, CASRN, ChemicalName, CV_Reps=Reps, Avg_CV_Rank, SNPs, Genes, Rules)]
pest.rank.out<- pest.rank.out[order(Disease, Rank)]

pest.rank.out$CASRN<- paste0('"', pest.rank.out$CASRN, '"')

setwd("")
fwrite(pest.rank.out, "SITables_Figures/Dataset6_Pesticide_Ranks.csv")
