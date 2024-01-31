############################################################################

### Overview of Pesticide-Pathway-Gene-SNP-Disease linkages
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library (plyr)
library(readxl)
library(grid)
library(gridExtra)
library(rgdal)

###NOTE: LOAD OWN WORKSPACE FROM STEP 5, OR WORK FROM EXISTING WORKSPACE (BELOW)

setwd("")
load("Workspaces/Analysis_Workspace.RData")
#Recall: gene match describes the actual gene the SNP is in - genesymbol/ID describe the c-g-p-d link

############################################################################
#Make barplot of linkage values and barplot of underlying numbers of unique items

item.dat<- chem.snp.dis.dat[, .(Pesticides=length(unique(CASRN)), Genes=length(unique(GeneID)),
                                SNPs=length(unique(na.omit(snpId))), Pathways=length(unique(PathwayID))), 
                            by=.()]
item.dat<- melt(item.dat)
item.dat$variable<- factor(item.dat$variable, levels=unique(item.dat$variable), ordered=TRUE)

#Data values described in text
cd.dat<- chem.snp.dis.dat[, .(.N), by=.(CASRN, UMLS_Name)] #Unique Chemical - Disease links
csd.dat<- chem.snp.dis.dat[!is.na(snpId), .(.N), by=.(CASRN, snpId, UMLS_Name)] 
csd.dat<- chem.snp.dis.dat[is.na(snpId) == FALSE] #unique Chemical - SNP - Disease links
cgd.dat<- chem.snp.dis.dat[, .(.N), by=.(CASRN, GeneID, UMLS_Name)] #unique Chemical - Gene - Disease links
cpgsd.dat<- chem.snp.dis.dat[, .(.N), by=.(CASRN, PathwayID, GeneID, snpId, UMLS_Name)] 
#unique chemical - SNP - Gene - Pathway - Disease links

#Median SNPs per pest - disease linkage
cd.dat<- chem.snp.dis.dat[!is.na(snpId), .(SNPs=length(unique(snpId))), by=.(CASRN, UMLS_Name)] 
quantile(cd.dat$SNPs)

nrow(chem.snp.dis.dat[is.na(snpId), .(.N), by=.(CASRN, GeneID, UMLS_Name)]) 

#Chemicals with most SNPs
c.dat<- chem.snp.dis.dat[!is.na(snpId), .(SNPs=length(unique(snpId))), by=.(CASRN, ChemicalName)]
c.dat<- c.dat[order(-SNPs)]
c.dat[1:10]

#Chemicals and SNPs per disease
d.dat<- chem.snp.dis.dat[, .(SNPs=length(unique(na.omit(snpId))), Chems=length(unique(CASRN))), by=.(UMLS_Name)]

############################################################################
#Number of items per category
#Figure 1

var.plot<- ggplot(item.dat) + geom_col(aes(x=variable, y=value, fill=variable)) +
  scale_fill_manual(values=c("midnightblue", "orange", "darkviolet", "forestgreen"), guide="none") + 
  scale_y_continuous(expand=c(0, 0), breaks=c(0, 300, 750, 1500, 2000)) + theme_bw() + 
  ylab(NULL) + xlab(NULL) +
  theme(panel.grid=element_blank(), panel.border=element_rect(color=NA, fill=NA), panel.background = element_blank(), 
        text=element_text(size=8, family="sans"), axis.ticks.x=element_blank(),
        axis.title.x= element_text(vjust=7, hjust=-.05, face="bold"), axis.title.y= element_text(face="bold", vjust=-0.2),
        axis.line = element_line())

jpeg("Figure1.jpeg", width=2, height=1.1, units="in", res=300)
var.plot
graphics.off()

setwd("")
fwrite(item.dat, "SITables_Figures/Figure1.csv")

############################################################################
#SI Maps

library(usmap)
library(usa)

state.size<- as.data.table(cbind(state.name, state.area, state.x19))

#KG pesticide/mi2
kg.plot.dat<- usgs.year.pest[, .(KG=sum(KG_APPLIED)), by=.(state=State)]
kg.plot.dat$Area<- as.numeric(as.character(state.size[match(kg.plot.dat$state, state.size$state.name)]$state.area))
kg.plot.dat<- kg.plot.dat[, .(Variable=KG/Area, Type="KG/Area"), by=.(state)]

#Pesticides applied per state
chem.plot.dat<- usgs.year.pest[, .(Variable=length(unique(COMPOUND)), Type="Pesticides"), by=.(state=State)]

#SNP hits per state
snp.plot.dat<- chem.snp.dis.dat[is.na(snpId)==FALSE, .(.N), by=.(snpId, CASRN)]
snp.plot.dat<- merge(snp.plot.dat, usgs.year.pest, by="CASRN", allow.cartesian=TRUE)
snp.plot.dat<- snp.plot.dat[, .(Variable=.N, Type="SNP Hits"), by=.(state=State)]

#Gene hits per state
gene.plot.dat<- chem.snp.dis.dat[Gene_Match != "", .(.N), by=.(GeneID, CASRN)]
gene.plot.dat<- merge(gene.plot.dat, usgs.year.pest, by="CASRN", allow.cartesian=TRUE)
gene.plot.dat<- gene.plot.dat[, .(Variable=.N, Type="Gene Hits"), by=.(state=State)]

kg.plot<- plot_usmap(data=kg.plot, values="Variable", exclude=c("Alaska", "Hawaii"), color="white") +
  theme_bw() + ggtitle(paste("KG/mi2 pesticide applied")) +
  scale_fill_gradient(low="seashell", high="red", guide = guide_colorbar(frame.colour = "black"), name="KG/mi2") + #coord_map() +
  theme(axis.text= element_blank(), axis.ticks= element_blank(), axis.title = element_blank(), 
        panel.border = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))

chem.plot<- plot_usmap(data=chem.plot, values="Variable", exclude=c("Alaska", "Hawaii"), color="white") +
  theme_bw() + ggtitle(paste("Number pesticides applied")) +
  scale_fill_gradient(low="seashell", high="red", guide = guide_colorbar(frame.colour = "black")) + #coord_map() +
  theme(axis.text= element_blank(), axis.ticks= element_blank(), axis.title = element_blank(), 
        panel.border = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))

snp.plot<- plot_usmap(data=snp.plot, values="Variable", exclude=c("Alaska", "Hawaii"), color="white") +
  theme_bw() + ggtitle(paste("Number SNP hits")) +
  scale_fill_gradient(low="seashell", high="red", guide = guide_colorbar(frame.colour = "black")) + #coord_map() +
  theme(axis.text= element_blank(), axis.ticks= element_blank(), axis.title = element_blank(), 
        panel.border = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))

gene.plot<- plot_usmap(data=gene.plot, values="Variable", exclude=c("Alaska", "Hawaii"), color="white") +
  theme_bw() + ggtitle(paste("Number gene hits")) +
  scale_fill_gradient(low="seashell", high="red", guide = guide_colorbar(frame.colour = "black")) + #coord_map() +
  theme(axis.text= element_blank(), axis.ticks= element_blank(), axis.title = element_blank(), 
        panel.border = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face="bold"), text=element_text(size=8, family="sans"))


setwd("")
jpeg("SITables_Figures/SIFigure3.jpeg", width=6, height=3.2, units="in", res=300)
grid.arrange(chem.plot, kg.plot,
             gene.plot, snp.plot, nrow=2)
graphics.off()

sifig3.dat<- rbind(kg.plot.dat, chem.plot.dat, snp.plot.dat, gene.plot.dat)

setwd("")
fwrite(sifig3.dat, "SITables_Figures/SI_FigureS3.csv")
