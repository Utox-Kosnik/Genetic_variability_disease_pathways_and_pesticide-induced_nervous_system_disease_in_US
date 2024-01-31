############################################################################

### Find chemical-pathway-SNP-disease associations
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
library(data.table)
library(tidyverse)
library (plyr)
library(readxl)


###NOTE: DOWNLOAD DISGENET SNP ASSOCIATION DATA
setwd("")
snp.dat<- fread(snp_associations...tsv) # (downloaded 221116)

setwd("")
chem.path.dis<- readRDS("Generated_Data/Full_Integrated_Chem_Gene_Path_Dis_Data.rds") #From step 4
chem.path.dis<- chem.path.dis[padj < 0.001 & FishP < 0.01 & CGYGDY >=5] #Based on sensitivity analysis


############################################################################
#Now bring in SNPs also 

library(disgenet2r)

#Add disgenet2r api key
disgenet_api_key <- get_disgenet_api_key(
  email = "", 
  password = "" )

Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

#Collect SNP data for each of the 6 diseases under study
dis.list<- as.character(unique(chem.path.dis$UMLS))
dat.dgn<- disease2variant(disease = dis.list, database = "ALL", score = c(0, 1))
dat.out<- dat.dgn@qresult

dat.out$Gene_Match<- dat.out$gene_symbol #Prepare to separate piped genes
dat.out<- as.data.table(separate_rows(dat.out, Gene_Match, sep=";"))

#Add chromosome and coordinate to SNPs
snp.dis<- merge(snp.dat, dat.out, by.x="snpId", by.y="variantid", allow.cartesian=TRUE)

############################################################################
#Now form linkages 

chem.path.dis$Dis_SNP_merge<- paste0(chem.path.dis$UMLS, "_", chem.path.dis$GeneSymbol)
snp.dis$Dis_SNP_merge<- paste0(snp.dis$diseaseid, "_", snp.dis$Gene_Match)

length(unique(chem.path.dis$CASRN)); length(unique(chem.path.dis$UMLS_Name))
length(unique(chem.path.dis$PathwayName)); length(unique(chem.path.dis$GeneSymbol))

snp.merge.dat<- snp.dis[Dis_SNP_merge %in% chem.path.dis$Dis_SNP_merge] #Reduce for "all" merge
chem.snp.dis<- merge(chem.path.dis, snp.merge.dat, by="Dis_SNP_merge", allow.cartesian = TRUE, all=TRUE)

chem.snp.dis[1:10] #Looks correct

############################################################################
##Now identify nearby SNPs not located in genes

length(which(snp.dis$Gene_Match == "")) #SNPs not in genes

#Set matches based on disease chromosome matches
snp.comp.dat<- snp.dis[Gene_Match == ""]
snp.comp.dat$Dis_SNP_merge<- paste0(snp.comp.dat$diseaseid, "_", snp.comp.dat$chromosome)

chem.snp.dis$Dis_SNP_merge<- paste0(chem.snp.dis$UMLS, "_", chem.snp.dis$chromosome)
snp.list<- unique(snp.comp.dat$snpId) #SNPs not in genes

#For each SNP not on a gene, determine if a snp is w/i 10k bases
dat.out<- NULL
for(snp.in in snp.list){
  print(which(snp.list %in% snp.in))
  
  dat.snp<- snp.comp.dat[snpId %in% snp.in]
  dat.chem<- chem.snp.dis[Dis_SNP_merge %in% dat.snp$Dis_SNP_merge]
  
  if(nrow(dat.chem)==0){
    next()
  } else{
    
    dat.snp$coord_high<- dat.snp$position + 10000
    dat.snp$coord_low<- dat.snp$position - 10000
    
    dat.comp<- merge(dat.chem, dat.snp, by="Dis_SNP_merge", allow.cartesian = TRUE)
    dat.comp<- dat.comp[position.x <= coord_high & position.x >= coord_low]
    
    if(nrow(dat.comp) == 0){
      next()
    } else{
      
      chem.fill<- dat.comp[,which(colnames(dat.comp) %in% colnames(chem.path.dis)), with=FALSE]
      colnames(chem.fill)[which(colnames(chem.fill)=="score")]<- "score.x"
      
      snp.fill<- dat.comp[, grep("\\.y$", colnames(dat.comp)), with=FALSE]
      colnames(snp.fill)<- gsub("\\.y$", "", colnames(snp.fill))
      colnames(snp.fill)[which(colnames(snp.fill)=="score")]<- "score.y"
      
      dat.fill<- as.data.table(cbind(chem.fill, snp.fill))
      dat.out<- as.data.table(rbind(dat.out, dat.fill))
    }

  }
  
}

#Now fix colnames to overlap between original and non-gene SNP dataset
length(which(!colnames(chem.snp.dis) %in% colnames(dat.out)))
colnames(dat.out)[which(!colnames(dat.out) %in% colnames(chem.snp.dis))]

colnames(dat.out)[which(colnames(dat.out) == "score.y")]<- "score"
dat.out<- dat.out[, match(colnames(chem.snp.dis), colnames(dat.out)), with=FALSE]

#Join files together
dat.out<- as.data.table(rbind(chem.snp.dis, dat.out))
dat.out<- dat.out[!duplicated(dat.out)]

#Dataset overview
length(unique(dat.out$CASRN)); length(unique(dat.out$UMLS_Name))
length(unique(dat.out$PathwayName)); length(unique(dat.out$GeneSymbol))
length(unique(dat.out$snpId))

test<- dat.out[is.na(snpId) == FALSE]
length(unique(test[Gene_Match == ""]$snpId)) 
length(unique(test[Gene_Match == ""]$CASRN)) 
unique(test[Gene_Match == ""]$UMLS_Name) 

###NOTE: SAVE DATA FOR FUTURE STEPS
setwd("")
saveRDS(dat.out, file="Generated_Data/Full_SNPJoined_Chem_Gene_Path_Dis_Data.rds")

############################################################################
#Build workspace for future analysis

rm(list = ls())

setwd("")
chem.snp.dis.dat<- readRDS("Generated_Data/Full_SNPJoined_Chem_Gene_Path_Dis_Data.rds") #from step 5
load("Workspaces/Step1_USGS_GBD_Data.RData") #from step 1

usgs.year.pest<- usgs.year.pest[KG_APPLIED != 0]
#Remove no application rows

###NOTE: SAVE WORKSPACE FOR FUTURE ANALYSIS
#Copy of workspace available currently as "Analysis_Workspace.RData" used in subsequent steps
setwd("")
save.image("Workspaces/FILENAME.RData")

############################################################################
#Export table of associations

chem.snp.dis.out<- chem.snp.dis.dat[, .(SNPs=paste0(na.omit(unique(snpId)), collapse=";"),
                                        Pathways=paste0(unique(PathwayID), collapse=";")),
                                    by=.(CASRN, ChemicalName, Disease=UMLS_Name, GeneID, GeneSymbol, Gene_Match)]
#Because GeneMatch indicates whether or not a SNP is in a gene, use that as gene reference

chem.snp.dis.out[Gene_Match == "" & SNPs != ""]$GeneID<- NA
chem.snp.dis.out[Gene_Match == "" & SNPs != ""]$GeneSymbol<- NA
#This means the SNP used as reference for IDing this intergenic SNP was in the gene given by GeneID/Symbol
#Remove to clarify it does not mean that SNP was in that gene

chem.snp.dis.out<- chem.snp.dis.out[, -"Gene_Match"]

chem.snp.dis.out$CASRN<- paste0('"', chem.snp.dis.out$CASRN, '"')

setwd("")
fwrite(chem.snp.dis.out, file="SITables_Figures/Dataset1-Chem-SNP-Dis_Table.csv")