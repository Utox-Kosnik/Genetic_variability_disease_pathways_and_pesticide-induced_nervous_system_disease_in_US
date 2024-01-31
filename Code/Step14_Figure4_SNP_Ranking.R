############################################################################

### SNP prioritization figure and list
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

#Only need priority SNPs
rm(list=setdiff(ls(), c("snp.high.cv")))

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)
setwd("")
load("Workspaces/Analysis_Workspace.RData")
chem.snp.dis.dat<- chem.snp.dis.dat[UMLS_Name %in% c("Alzheimer's Disease", "Multiple Sclerosis", "Parkinson Disease")]

##########################################################################
#Add SNP data to ranked SNPs

#Prep SNPs not in genes for merging with the SNP data
chem.snp.dis.merge<- chem.snp.dis.dat[, .(.N), by=.(snpId, GeneID, GeneSymbol, Gene_Match, Disease=UMLS_Name, 
                                                    chromosome, position, 
                                                    variant_consequence_type)]
chem.snp.dis.merge[Gene_Match == ""]$GeneID<- NA
chem.snp.dis.merge[Gene_Match == ""]$GeneSymbol<- ""
#In cases where the SNP is intergenic, GeneID/Symbol refer to the gene with the original SNP IDd - fix

snp.high.full.dat<- merge(snp.high.cv, chem.snp.dis.merge[, list(snpId, GeneID, GeneSymbol, Disease,
                                                                 chromosome, position, variant_consequence_type)],
                          by=c("snpId", "GeneID", "GeneSymbol", "Disease"))

snp.high.full.dat<- snp.high.full.dat[!duplicated(snp.high.full.dat)] #Some duplicates from gene_match

#https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
var.class<- data.table(c("splice acceptor variant", "splice donor variant", "stop gained", "frameshift variant", 
                         "stop lost", "start lost",
                         "inframe insertion", "inframe deletion", "missense variant", "protein altering variant",
                         "splice region variant", "synonymous variant",
                         "coding sequence variant", "mature miRNA variant", "5 prime UTR variant", "3 prime UTR variant", 
                         "non coding transcript exon variant", "intron variant", "non coding transcript variant",
                         "upstream gene variant", "downstream gene variant", "Intergenic variant"))
var.class$Severity<- c("High", "High", "High", "High", "High", "High", 
                       "Moderate", "Moderate", "Moderate", "Moderate", 
                       "Low", "Low", 
                       "Modifier", "Modifier", "Modifier", "Modifier", "Modifier", "Modifier", "Modifier", "Modifier",
                       "Modifier", "Modifier")
snp.high.full.dat$Severity<- var.class[match(snp.high.full.dat$variant_consequence_type, var.class$V1)]$Severity

##########################################################################
###Look at SNP patterns

#62
length(unique(snp.high.full.dat[is.na(snpId) == FALSE]$snpId))/length(unique(chem.snp.dis.dat[is.na(snpId) == FALSE]$snpId))
length(unique(snp.high.full.dat$GeneSymbol))/length(unique(chem.snp.dis.dat$GeneSymbol)) #69

#Number SNPs lower hits in high states
length(which(!unique(chem.snp.dis.dat$snpId) %in% unique(snp.high.full.dat$snpId)))
#Number genes lower hits in high states
length(which(!unique(chem.snp.dis.dat$GeneID) %in% unique(snp.high.full.dat$GeneID)))

###Look at top SNPs in each disease
###Recall: RankOrder = Rank (total rank across 10 reps)/Reps (number of reps)
###Alzheimer's
test<- snp.high.full.dat[Disease %in% "Alzheimer's Disease" & !is.na(snpId)]
test<- test[order(RankOrder)]
test[1:15]
#moderate and modifier

#MS
test<- snp.high.full.dat[Disease %in% "Multiple Sclerosis" & !is.na(snpId)]
test<- test[order(RankOrder)]
test[1:15]
#modifier, moderate, 10th is high

###Parkinson
test<- snp.high.full.dat[Disease %in% "Parkinson Disease" & !is.na(snpId)]
test<- test[order(RankOrder)]
test[1:15]
#All 10 top genes moderate severity

##################################################################################
###Determine overlapping SNPs between diseases
#Overall ranking is ranking based on rules multiplied by ranking based on hits

high.snp.plot<- snp.high.full.dat[is.na(snpId)==FALSE]

#How many diseases implicated per SNP?
dis.fill<- high.snp.plot[, .(length(unique(Disease))), by=.(snpId)]
table(dis.fill$V1) #Most unique, look at overlapping SNPs

test<- high.snp.plot[snpId %in% dis.fill[V1==2]$snpId]
test[, .(length(unique(snpId))), by=.(Disease)] #Mostly Alzheimer's and Parkinson

test<- high.snp.plot[snpId %in% dis.fill[V1==3]$snpId]
test[, .(length(unique(snpId))), by=.(Disease)] #Lower ranks

##################################################################################
###Generate Manhattan plot
#Figure 4a

#Set shape based on number of overlapping diseases
high.snp.plot$Shape<- high.snp.plot$Disease
high.snp.plot[snpId %in% dis.fill[V1 > 1]$snpId]$Shape<- "Many"

quantile(high.snp.plot$RankOrder)

#Assign correct position of SNPs on full axis
high.snp.plot<- high.snp.plot[chromosome %in% c(1:22)]
high.snp.plot$chromosome<- as.numeric(as.character(high.snp.plot$chromosome))
high.snp.plot$BPfull<- 0
s<- 0

#Create continuous set of positions along all chromosomes
for(i in c(1:22)){
  val<- max(high.snp.plot[chromosome==i]$position)
  if(val==-Inf){
    next()
  }
  high.snp.plot[chromosome==i]$BPfull<- high.snp.plot[chromosome==i]$position + s
  s<- s+val
}


high.snp.plot$Disease <- factor(high.snp.plot$Disease, levels = c("Alzheimer's Disease", "Multiple Sclerosis", 
                                                                  "Parkinson Disease"), ordered=TRUE) 

#Set shape of off number of diseases SNP present in 
high.snp.plot[Shape %in% "Many"]$Shape <- "8"
high.snp.plot[!Shape %in% "8"]$Shape <- "19" 
high.snp.plot$Shape<- factor(high.snp.plot$Shape, levels = c(19, 8), ordered=TRUE)

#Set size off of severity
high.snp.plot$SeverityPlot<- 2
high.snp.plot[Severity %in% "Low"]$SeverityPlot<- 1
high.snp.plot[Severity %in% "Moderate"]$SeverityPlot<- 3
high.snp.plot[Severity %in% "High"]$SeverityPlot<- 4
high.snp.plot$SeverityPlot<- factor(high.snp.plot$SeverityPlot, levels=c(1, 2, 3, 4), ordered=TRUE)

#Add labels to top SNPs
label.list<- rbind(head(snp.high.cv[is.na(snpId) == FALSE & Disease %in% "Alzheimer's Disease"], 3),
                   head(snp.high.cv[is.na(snpId) == FALSE & Disease %in% "Multiple Sclerosis"], 3),
                   head(snp.high.cv[is.na(snpId) == FALSE & Disease %in% "Parkinson Disease"], 3))

high.snp.plot$Label<- high.snp.plot$snpId
high.snp.plot[!Label %in% label.list$snpId]$Label<- ""
high.snp.plot[Label != ""] #All different

#Prep position of labels along axis
axis.set <- high.snp.plot [,.(center = (max(BPfull) + min(BPfull)) / 2), by=.(chromosome)]

library(ggrepel)
setwd("")
jpeg("Figure4a.jpg", width=7.5, height=2.7, units="in", res=300)
ggplot(high.snp.plot, aes(x=BPfull, y=1/RankOrder, label=Label)) +  
  geom_point(aes(color=Disease, shape=Shape, size=SeverityPlot)) + geom_text_repel(size=2.25) +
  scale_color_manual(values=c("Alzheimer's Disease"="mediumpurple3", "Multiple Sclerosis"="goldenrod2", 
                              "Parkinson Disease"="darkcyan", "Many"="grey50")) +
  scale_size_manual(values=c(1, 2, 3, 4), labels=c("4"="High", "3"="Moderate", "2"="Modifier",
                                                   "1"="Low")) +
  scale_shape_manual(values=c(19, 8), labels=c("19"="One Disease", "8"="Many Diseases")) +
  guides(color=guide_legend(nrow=2), size=guide_legend(nrow=2), shape=guide_legend(nrow=2)) + 
  scale_y_continuous(trans="log10") +     
  scale_x_continuous(expand=c(.01, .01), label = axis.set$chromosome, breaks = axis.set$center) +
  geom_hline(aes(yintercept=0), color="black", linewidth=0.8, linetype="solid")+
  ylab("Rank values") + xlab("Chromosome") +
  theme_bw() + 
  theme(panel.grid=element_blank(), panel.border=element_rect(color=NA, fill=NA), panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), legend.position="bottom",
        axis.title.x= element_text(vjust=7, hjust=-.085, face="bold"), axis.title.y= element_text(face="bold", vjust=-0.2),
        axis.line = element_line())
graphics.off()

setwd("")
fwrite(high.snp.plot, "SITables_Figures/Figure4a_SNP_Manhattan_Plot.csv")

##################################################################################
###Generate output table

snp.dat.out<- snp.high.full.dat

snp.dat.out$Value<- 1
snp.dat.out[Disease %in% "Alzheimer's Disease"]$Value<- 
  rank(snp.dat.out[Disease %in% "Alzheimer's Disease"]$RankOrder, ties.method=c("average"))
snp.dat.out[Disease %in% "Multiple Sclerosis"]$Value<- 
  rank(snp.dat.out[Disease %in% "Multiple Sclerosis"]$RankOrder, ties.method=c("average"))
snp.dat.out[Disease %in% "Parkinson Disease"]$Value<- 
  rank(snp.dat.out[Disease %in% "Parkinson Disease"]$RankOrder, ties.method=c("average"))
snp.dat.out<- snp.dat.out[order(Disease, Value)]

snp.dat.out<- snp.dat.out[, list(Rank=Value, Disease, snpId, GeneID, GeneSymbol, CV_Reps=Reps, 
                                 Avg_CV_Rank=Rank/Reps, Avg_CV_HighHits=HighHits, Avg_CV_Rules=Rules)]


setwd("")
fwrite(snp.dat.out, "SITables_Figures/Dataset4_SNP_Ranks.csv")
