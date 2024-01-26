############################################################################

### Find chemical-pathway-disease associations
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
require(data.table)
require(tidyverse)
require(plyr)
require(readxl)
require(gridExtra)

###NOTE: DOWNLOAD DISGENET GENE ASSOCIATION DATA
setwd("")
gene.dat<- fread(gene_associations...tsv) # (downloaded 221116)

setwd("")
load("Workspaces/Step1_USGS_GBD_Data.RData") #From step 1
chem.gene<- fread("Generated_Data/Integrated_HTS_CTD_ChemGene_Links.csv") #From Step 2
path.dis<- fread("Generated_Data/Full_Enriched_Path_Disease_Links.csv") #From Step 3

path.dis<- path.dis[Count >=10 & FDR <.01] #Based on 1st sensitivity analysis, this will be min req

path.dis<- path.dis[, list(GeneID=shared_geneid, GeneSymbol=shared_symbol, PathwayID, PathwayName, UMLS=ID, 
                           UMLS_Name=Description, OverlapSize=Count, padj=FDR)]

############################################################################
#Now form linkages

cas.list<- unique(c(chem.gene$CASRN))
dat.out<- NULL

#For each chemical set, merge gene associations with genes in a path - dis intersection
#Keeping all linkages for now, will then assess using Fisher's test
for(cas.in in cas.list){
  
  dat.in<- chem.gene[CASRN %in% cas.in]
  dat.inter<- merge(dat.in, path.dis, by=c("GeneID", "GeneSymbol"), allow.cartesian=TRUE)
  dat.out<- rbind(dat.out, dat.inter)
  
  gc()
  
  print(which(cas.list == cas.in))
}

############################################################################
#Now do Fisher's test with these associations

#Have chemical-gene data integrated with pathway-disease data
#Set minimum overlap b/w chemical/pathway to count as chemical/pathway/disease link
#Disease associations needed at least 10 genes overlapping with a pathway

#######g.p.d        ng.p.d
#c.g: dat.out      int.dat - dat.out
#nc.g: dis-dat.out gene.dis-all

#genedis data from DisGeNET = background, all genes implicated across diseases


#Reduce disease data to the 6 diseases present in GBD
dat.out<- dat.out[UMLS_Name %in% c("Alzheimer's Disease", "Multiple Sclerosis", 
                                   "Parkinson Disease", "Brain Neoplasms", "Epilepsy",
                                   "Migraine Disorders")]

dat.out$CD<- paste0(dat.out$CASRN, "_", dat.out$UMLS)
#Set each chemical - disease association as the basis of the enrichment

cd.list<- unique(dat.out$CD)

fish.out<- as.data.table(matrix(data=0, nrow=length(cd.list), ncol=6, 
                                dimnames = list(NULL, c("CD", 
                                                        "CGYGDY", #chem gene and gene dis overlap
                                                        "CGYGDN", #chem gene and no gene dis overlap
                                                        "CGNGDY", #no chem gene and gene dis overlap
                                                        "CGNGDN", #no chem gene and no gene dis overlap
                                                        "FishP"))))
fish.out$CD<- cd.list
for(i in c(1:length(cd.list))){
  #Fishers test for all association
  cd.in<- cd.list[i] 
  
  #chem gene and gene dis overlap
  cg.yes.gd.yes<- dat.out[CD %in% cd.in] 
  
  #chem gene and no gene dis overlap
  cg.yes.gd.no<- chem.gene[CASRN %in% cg.yes.gd.yes$CASRN]
  cg.yes.gd.no<- cg.yes.gd.no[!GeneID %in% cg.yes.gd.yes$GeneID]
  
  #no chem gene and gene dis overlap 
  cg.no.gd.yes<- path.dis[UMLS %in% cg.yes.gd.yes$UMLS]
  cg.no.gd.yes<- cg.no.gd.yes[!GeneID %in% cg.yes.gd.yes$GeneID]
  
  #no chem gene and no gene dis overlap - whole spectrum genes implicated in diseases
  back.dat<- gene.dat[!geneId %in% unique(c(cg.yes.gd.yes$GeneID, cg.yes.gd.no$GeneID, cg.no.gd.yes$GeneID))]
  
  fish.dat<- data.frame(c(length(unique(cg.yes.gd.yes$GeneID)), length(unique(cg.yes.gd.no$GeneID))), 
                        c(length(unique(cg.no.gd.yes$GeneID)), length(unique(back.dat$geneId))))
  
  fish.val<- fisher.test(fish.dat)
  
  fish.out[i]$CGYGDY<- length(unique(cg.yes.gd.yes$GeneID))
  fish.out[i]$CGYGDN<- length(unique(cg.yes.gd.no$GeneID))
  fish.out[i]$CGNGDY<- length(unique(cg.no.gd.yes$GeneID))
  fish.out[i]$CGNGDN<- length(unique(back.dat$geneId))
  
  fish.out[i]$FishP<- fish.val$p.value
  
}

fish.red<- fish.out[FishP < 0.05 & CGYGDY >=3] #min reqs p <0.05 and at least 3 genes overlapping
#

cgpd.dat<- merge(dat.out, fish.red, by="CD", allow.cartesian=TRUE)

###NOTE: SAVE DATA FOR FUTURE STEPS
setwd("")
saveRDS(cgpd.dat, file="Generated_Data/Full_Integrated_Chem_Gene_Path_Dis_Data.rds")

############################################################################
###SI Figure 2

setwd("")
cgpd.dat<- readRDS("Generated_Data/Full_Integrated_Chem_Gene_Path_Dis_Data.rds")

#Plot differences in number of data values based on different restrictions
#pathway-disease p values .01, .001
#pathway-gene-disease overlap = 10, 20
#chemical-disease p values .001, .01, .05
#chemical-gene-disease overlap = 3, 5, 10, 20

#Prepare different datasets and put together

range(cgpd.dat$padj)
range(cgpd.dat$OverlapSize)
range(cgpd.dat$FishP)
range(cgpd.dat$CGYGDY)


#for each combination of p-value and gene count, determine relevant data values
p.count<- function(in.dat, p.value, g.value.1, g.value.2, p.print){
  
  #restrict dataset - input data will be restricted by path-dis p-value condition
  int.val<- in.dat[FishP < p.value & OverlapSize >=g.value.1 & CGYGDY >= g.value.2, 
                      list(UMLS_Name, CASRN, PathwayID, GeneID)]
  
  #counts of unique items
  int.un.count<- int.val[, .(chemicals=length(unique(CASRN)), genes=length(unique(GeneID)), pathways=length(unique(PathwayID)))]
  
  #linkage counts
  cd.link<- int.val[, .(.N), by=.(CASRN, UMLS_Name)]
  cp.link<- int.val[, .(.N), by=.(CASRN, PathwayID)]
  cpd.link<- int.val[, .(.N), by=.(CASRN, PathwayID, UMLS_Name)]
  
  int.un.count<- cbind(int.un.count, data.table(`chem-dis`=nrow(cd.link)), data.table(`chem-path`=nrow(cp.link)), 
                       data.table(`chem-path-dis`=nrow(cpd.link)))
  
  #Put values together
  dat.out<- melt(int.un.count)
  dat.out$variable<- paste("Unique", dat.out$variable)
  
  dat.out$Pvalue<- p.print
  dat.out$DPGenes<- g.value.1
  dat.out$CDGenes<- g.value.2
  
  dat.out
}

full.plot<- NULL
p.list<- c(0.001, 0.01, 0.05)
g.list1<- c(10, 20)
g.list2<- c(3, 5, 10, 20)

for(p.in in p.list){
  for(g.in1 in g.list1){
    for(g.in2 in g.list2){
      
    #For each pvalue - gene overlap combination, compute unique paths, genes, and linkages
    full.plot.p1<- p.count(cgpd.dat[padj< 0.01], p.in, g.in1, g.in2, paste("Fisher p <", p.in))
    full.plot.p1$ORAp<- 0.01
    
    full.plot.p2<- p.count(cgpd.dat[padj< 0.001], p.in, g.in1, g.in2, paste("Fisher p <", p.in))
    full.plot.p2$ORAp<- 0.001
    
    full.plot<- rbind(full.plot, full.plot.p1, full.plot.p2)
    }
  }
  
}

full.plot$Pvalue<- factor(full.plot$Pvalue, levels= unique(full.plot$Pvalue), ordered=TRUE)
full.plot$DPGenes<- as.character(full.plot$DPGenes)

plot.p1<- ggplot(full.plot[ORAp==0.01]) + geom_point(aes(x=CDGenes, y=value, color=CDGenes, shape=DPGenes)) + theme_bw() +
  facet_grid(rows=vars(variable), cols=vars(Pvalue), scales="free_y") + ggtitle("ORA p < 0.01") +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), text=element_text(size=6))

plot.p2<- ggplot(full.plot[ORAp==0.001]) + geom_point(aes(x=CDGenes, y=value, color=CDGenes, shape=DPGenes)) + theme_bw() +
  facet_grid(rows=vars(variable), cols=vars(Pvalue), scales="free_y") + ggtitle("ORA p < 0.001") +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5), text=element_text(size=6))


setwd("")
jpeg("SITables_Figures/FigureS2.jpeg", width=6.5, height=6, units="in", res=300)
grid.arrange(plot.p1, plot.p2, ncol=2)
graphics.off()

setwd("")
fwrite(full.plot, "SITables_Figures/SI_FigureS2.csv")

