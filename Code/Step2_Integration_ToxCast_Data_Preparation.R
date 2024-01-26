############################################################################

### ToxCast/Tox21 active chem-gene endpoint determination 
### Integration with CTD

### Author: Marissa Kosnik

############################################################################

rm(list = ls())
require(data.table)
require(readxl)
require(tidyverse)

###NOTE: DOWNLOAD CTD AND TOX21/TOXCAST DATA AS DESCRIBED

setwd("")
ctd.dat<- fread("CTD_chem_gene_ixns.tsv") #231129 release of CTD chemical - gene interactions data
###CTD chemical - gene interactions
###http://ctdbase.org/downloads

###Tox21/ToxCast data: invitrodb_v4_1 files
###https://www.epa.gov/chemical-research/exploring-toxcast-data-downloadable-data
setwd("")
mod.dat<- fread("mc5-6_winning_model_fits-flags_invitrodb_v4_1_SEPT2023.csv")
cyto.dat<- as.data.table(read_excel("cytotox_invitrodb_v4_1_SEPT2023.xlsx"))
gene.dat<- as.data.table(read_excel("assay_gene_mappings_invitrodb_v4_1_SEPT2023.xlsx")) 
assay.dat<- as.data.table(read_excel("assay_annotations_invitrodb_v4_1_SEPT2023.xlsx"))

############################################################################
#Select active chemicals in Tox21/ToxCast

act.dat<- mod.dat[hitc >= 0.9] #Extract active aeids (active assay - endpoint hits)
act.dat<- act.dat[is.na(chid)==FALSE] #Remove those w/o a chemical ID (DMSO, steroids, etc.)
act.dat$logAC50<- log10(act.dat$ac50) #Determine logAC50

############################################################################
#Sort genes by species
#Bridge across all

head(gene.dat)

#Select for genes
length(unique(gene.dat$aeid)) == nrow(gene.dat) #Some aeids listed for more genes
gene.res<- gene.dat[, .(.N), by=.(entrez_gene_id, gene_symbol)] #use to match symbols later
gene.dat<- gene.dat[, list(aeid, entrez_gene_id)]
test<- gene.dat[, .(.N), by=.(aeid)] #how many genes per assay?
test<- test[N >1] #Assays with more than one gene target
test<- gene.dat[aeid %in% test$aeid] #Get the data for these assays
test<- test[, .(entrez_gene_id=paste(entrez_gene_id, collapse= "|")), by=.(aeid)] #Pipe overlapping genes 
#Replace separated aeid-gene entries with piped entries
gene.dat<- gene.dat[!aeid %in% test$aeid]
gene.dat<- as.data.table(rbind(gene.dat, test))
gene.dat<- gene.dat[!duplicated(gene.dat)]

act.dat<- merge(act.dat, gene.dat, by="aeid") #Determine whether or not the assay was active

#Now split piped gene data
act.dat<- as.data.table(separate_rows(act.dat, entrez_gene_id, sep="[|]"))
act.dat<- act.dat[!duplicated(act.dat)]
length(which(!act.dat$entrez_gene_id %in% gene.res$entrez_gene_id)) # should be 0
act.dat$gene_symbol<- gene.res[match(act.dat$entrez_gene_id, gene.res$entrez_gene_id)]$gene_symbol
#Add gene symbol

#Make gene symbols uniform
act.dat$gene_symbol<- casefold(act.dat$gene_symbol, upper=TRUE)
############################################################################
#Determine general cytotoxicity value

head(cyto.dat)

#New cytotox values don't cover all chemicals - using general cytotox/viability tests as general, backup cytotox values
cyto.assays<- act.dat[aeid %in% assay.dat[assay_function_type == "cytotoxicity"]$aeid]
cell.vi.assays<- act.dat[aeid %in% assay.dat[cell_viability_assay == 1]$aeid]
test<- rbind(cyto.assays, cell.vi.assays)
test<- test[!duplicated(test)]
gen.cyt.val<- median(test$logAC50)

############################################################################
#Remove instances of cytotox activity - determined as cytotoxicity conc lower than active conc

length(unique(cyto.dat$chid)) == nrow(cyto.dat) #should be true - 1 cytotox value per chid

quantile(na.omit(cyto.dat$cytotox_median_log)) #Compared to burst cytotox value, gen value <25th quartile

#cytotox_median_log = The cytotoxicity point, or the value in “med” when “nhit” is at least 5% of the total assay endpoints
test<- cyto.dat[match(act.dat$chid, cyto.dat$chid)]
act.dat$CytologAC50<- test$cytotox_median_log #Add cytotoxicity value
act.dat[is.na(CytologAC50)]$CytologAC50<- gen.cyt.val #When no cytotox value provided, use general value

#Cytotox - chemical-endpoint parameters must be less than cytotox value
act.dat2<- act.dat[logAC50 < CytologAC50]
act.dat2<- act.dat2[is.na(casn) == FALSE]
act.dat2<- act.dat2[is.na(entrez_gene_id) == FALSE]

test<- act.dat2[, .(.N), by=.(casn, entrez_gene_id)]
length(unique(test$casn))
length(unique(test$entrez_gene_id))

############################################################################
#Prep CTD data

ctd.dat<- ctd.dat[CasRN != ""]
ctd.dat<- ctd.dat[Organism == "Homo sapiens"] #Only keep human linkages
length(unique(ctd.dat$CasRN))
length(unique(ctd.dat$GeneID))
############################################################################
#Merge Tox21/ToxCast and CTD

hts.dat<- act.dat

rm(list=setdiff(ls(), c("hts.dat", "ctd.dat")))

length(which(unique(hts.dat$casn) %in% ctd.dat$CasRN)) #Most do not overlap
length(which(unique(hts.dat$entrez_gene_id) %in% ctd.dat$GeneID)) #Most overlap

###For full info on the two, will integrate with original sources kept identifiable
hts.dat$HTS_Source<- "Y"
ctd.dat$CTD_Source<- "Y"

hts.res<- hts.dat #Keep for chemical names
hts.dat<- hts.dat[, list(CASRN = casn, GeneID=entrez_gene_id, GeneSymbol=gene_symbol, HTS_Source)]
ctd.dat<- ctd.dat[, list(CASRN=CasRN, ChemicalName, GeneID, GeneSymbol, CTD_Source)]

hts.dat$GeneID<- as.character(hts.dat$GeneID)
ctd.dat$GeneID<- as.character(ctd.dat$GeneID)

#Fully integrate data
int.dat.full<- merge(hts.dat, ctd.dat, by=c("CASRN", "GeneID", "GeneSymbol"), all=TRUE)
int.dat.full[is.na(HTS_Source)]$HTS_Source<- "X"
int.dat.full[is.na(CTD_Source)]$CTD_Source<- "X"
int.dat.full$Both_Source<- "X"
int.dat.full[HTS_Source=="Y" & CTD_Source=="Y"]$Both_Source<- "Y"

#Add HTS names where needed - priority to CTD names
hts.res<- hts.res[match(int.dat.full$CASRN, hts.res$casn)]
int.dat.full[is.na(ChemicalName)]$ChemicalName<- hts.res[is.na(int.dat.full$ChemicalName)]$chnm

#Data checks
length(which(is.na(int.dat.full$GeneID) == TRUE))
int.dat.full[is.na(GeneID)==TRUE] #Relic of NA being treated as a unique gene
int.dat.full<- int.dat.full[is.na(GeneID)==FALSE]
length(which(is.na(int.dat.full$ChemicalName)==TRUE))
length(unique(int.dat.full$CASRN))
length(unique(int.dat.full$GeneID))

#Check some data linkages

int.dat.full[Both_Source == "Y"]
table(int.dat.full[Both_Source == "Y"]$HTS_Source) #should all be Y
table(int.dat.full[Both_Source == "Y"]$CTD_Source) #should all be Y
both.cont<- int.dat.full[Both_Source == "Y", .(.N), by=.(CASRN, GeneID)]
hts.cont<- int.dat.full[Both_Source == "X" & HTS_Source == "Y", .(.N), by=.(CASRN, GeneID, CTD_Source)]
table(hts.cont$CTD_Source) #should all be x
ctd.cont<- int.dat.full[Both_Source == "X" & CTD_Source == "Y", .(.N), by=.(CASRN, GeneID, HTS_Source)]
table(ctd.cont$HTS_Source) #should all be x

int.dat.full<- int.dat.full[!duplicated(int.dat.full)] 
#Now have integrated CTD and Tox21/ToxCast chemical-gene linkages

############################################################################
#Reduce dataset to pesticides and harmonize chemical names

setwd("")
load("Workspaces/Step1_USGS_GBD_Data.RData")
#Data from workspace step 1
#usgs.year.pest = USGS pesticides used per state per year
#gbd.dat = GBD disease rates per state per year

int.dat.pest<- int.dat.full[CASRN %in% unique(usgs.year.pest$CASRN)]

#Turn all chemical names upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

int.dat.pest$ChemicalName<- firstup(int.dat.pest$ChemicalName)
int.dat.pest<- int.dat.pest[!ChemicalName %in% "Zinc"] #Remove Zinc - too many therapeutic associations
test<- int.dat.pest[, .(length(unique(ChemicalName))), by=.(CASRN)]
test[V1 > 1]

#Fix these pesticides with different chemical names - give priority to CTD/shorter name
#Names mismatched between ToxCast and CTD
int.dat.pest[CASRN %in% "15299-99-7"]$ChemicalName<- "Devrinol"
int.dat.pest[CASRN %in% "27314-13-2"]$ChemicalName<- "Norflurazon"
int.dat.pest[CASRN %in% "298-00-0"]$ChemicalName<- "Methyl Parathion"
int.dat.pest[CASRN %in% "51-03-6"]$ChemicalName<- "Piperonyl Butoxide"
int.dat.pest[CASRN %in% "23950-58-5"]$ChemicalName<- "Pronamide"
int.dat.pest[CASRN %in% "59756-60-4"]$ChemicalName<- "Fluridone"
int.dat.pest[CASRN %in% "68694-11-1"]$ChemicalName<- "Triflumizole"
int.dat.pest[CASRN %in% "9006-42-2"]$ChemicalName<- "Polymarcine"
int.dat.pest[CASRN %in% "82-68-8"]$ChemicalName<- "Quintozene"
int.dat.pest[CASRN %in% "1897-45-6"]$ChemicalName<- "Tetrachloroisophthalonitrile"
int.dat.pest[CASRN %in% "137-26-8"]$ChemicalName<- "Thiram"
int.dat.pest[CASRN %in% "2312-35-8"]$ChemicalName<- "Omite"
int.dat.pest[CASRN %in% "28249-77-6"]$ChemicalName<- "Benthiocarb"
int.dat.pest[CASRN %in% "93-65-2"]$ChemicalName<- "Mecoprop"
int.dat.pest[CASRN %in% "94-75-7"]$ChemicalName<- "2,4-Dichlorophenoxyacetic Acid"
int.dat.pest[CASRN %in% "709-98-8"]$ChemicalName<- "Propanil"
int.dat.pest[CASRN %in% "142459-58-3"]$ChemicalName<- "FOE 5043"
int.dat.pest[CASRN %in% "55283-68-6"]$ChemicalName<- "Ethafluralin"
int.dat.pest[CASRN %in% "7287-19-6"]$ChemicalName<- "Prometryne"
int.dat.pest[CASRN %in% "94-74-6"]$ChemicalName<- "2-Methyl-4-chlorophenoxyacetic Acid"
int.dat.pest[CASRN %in% "123-33-1"]$ChemicalName<- "Maleic Hydrazide"
int.dat.pest[CASRN %in% "23564-05-8"]$ChemicalName<- "Thiophanate"
int.dat.pest[CASRN %in% "112-05-0"]$ChemicalName<- "Pelargonic Acid"
int.dat.pest[CASRN %in% "112-30-1"]$ChemicalName<- "N-decyl Alcohol"
int.dat.pest[CASRN %in% "1610-18-0"]$ChemicalName<- "Prometone"
int.dat.pest[CASRN %in% "2104-64-5"]$ChemicalName<- "Phenylphosphonothioic Acid, 2-Ethyl 2-(4-Nitrophenyl) Ester"
int.dat.pest[CASRN %in% "42509-80-8"]$ChemicalName<- "Isazophos"
int.dat.pest[CASRN %in% "145-73-3"]$ChemicalName<- "Endothall"
int.dat.pest[CASRN %in% "86-50-0"]$ChemicalName<- "Azinphosmethyl"
int.dat.pest[CASRN %in% "78-83-1"]$ChemicalName<- "Isobutyl alcohol"
int.dat.pest[CASRN %in% "886-50-0"]$ChemicalName<- "Terbutryne"
int.dat.pest[CASRN %in% "3383-96-8"]$ChemicalName<- "Temefos"
int.dat.pest[CASRN %in% "26644-46-2"]$ChemicalName<- "Triforine"
int.dat.pest[CASRN %in% "76-03-9"]$ChemicalName<- "Trichloroacetic Acid"
int.dat.pest[CASRN %in% "834-12-8"]$ChemicalName<- "Ametryne"
int.dat.pest[CASRN %in% "1194-65-6"]$ChemicalName<- "Dichlobanil"
int.dat.pest[CASRN %in% "16672-87-0"]$ChemicalName<- "Ethephon"
int.dat.pest[CASRN %in% "2303-17-5"]$ChemicalName<- "Triallate"
int.dat.pest[CASRN %in% "759-94-4"]$ChemicalName<- "EPTC"

length(unique(int.dat.pest$CASRN))
length(unique(int.dat.pest$GeneID))

#Now have integrated CTD and Tox21/ToxCast chemical-gene linkages

###NOTE: SAVE FILE FOR FUTURE STEPS
setwd("")
fwrite(int.dat.pest, file="Generated_Data/Integrated_HTS_CTD_ChemGene_Links.csv")

