############################################################################

### Find pathway-disease associations
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())
require(data.table)
require(tidyverse)
require (plyr)
require(readxl)
require(disgenet2r)
require(pathfindR)
require(KEGGREST)

###NOTE: DOWNLOAD DISGENET GENE ASSOCIATION DATA
setwd("")
gene.dat<- fread(gene_associations...tsv) # (downloaded 221116)

###NOTE: DOWNLOAD DAVID KNOWLEDGEBASE FILES AS DESCRIBED

#https://david.ncifcrf.gov/knowledgebase/DAVID_knowledgebase.html
#homo sapiens, central identifier = Entrez gene ID, reactome and wikipathways, 
setwd("")
david.path.reac<- fread("ENTREZ_GENE_ID2REACTOME_PATHWAY.txt") # (downloaded 231206)
david.path.wiki<- fread("ENTREZ_GENE_ID2WIKIPATHWAYS.txt") # (downloaded 231206)

############################################################################
#Clean DAVID pathway data

david.path.reac$PathwayName<- gsub(".*~", "", david.path.reac$V2)
david.path.reac$V2<- gsub("~.*", "", david.path.reac$V2)
david.path.wiki$PathwayName<- gsub(".*~", "", david.path.wiki$V2)
david.path.wiki$V2<- gsub("~.*", "", david.path.wiki$V2)
colnames(david.path.reac)[1:2]<- c("GeneID", "PathwayID")
colnames(david.path.wiki)[1:2]<- c("GeneID", "PathwayID")


############################################################################
#Import and clean KEGG pathway data

#pathfindR KEGG list
kegg.list<- get_gene_sets_list(source="KEGG", org_code = "hsa")
kegg.res<- kegg.list[2] 

#Join names, paths, genes into a flat format
#Prep path IDs
kegg.nams<- matrix(data=NA, nrow=length(kegg.res[1][[1]]), ncol=2, dimnames = list(NULL, c("PathwayID", "PathwayName")))
for(i in c(1:length(kegg.res[1][[1]]))){
  val<- kegg.res[1][[1]][i]
  kegg.nams[i,]<- c(names(val), val)
}
kegg.nams<- as.data.table(kegg.nams)

#Prep path genes
kegg.res<- kegg.list[1]
kegg.paths<- NULL
for(i in c(1:length(kegg.res[1][[1]]))){
  val<- as.data.table(kegg.res[1][[1]][i])
  val<- as.data.table(cbind(rep(colnames(val, nrow(val))), val))
  colnames(val)<- c("PathwayID", "GeneSymbol")
  kegg.paths<- as.data.table(rbind(kegg.paths, val))
}

kegg.paths$PathwayName<- kegg.nams[match(kegg.paths$PathwayID, kegg.nams$PathwayID)]$PathwayName
kegg.paths$GeneID<- gene.dat[match(kegg.paths$GeneSymbol, gene.dat$geneSymbol)]$geneId

############################################################################
#Now join pathway data into one file

path.enrich<- rbind(kegg.paths[, list(GeneID, PathwayID, PathwayName)], david.path.reac, david.path.wiki)

############################################################################
#https://www.disgenet.org/static/disgenet2r/disgenet2r.html#performing-a-disease-enrichment
#DisGeNET enrichment
#NOTE: enter email and password for package access

disgenet_api_key <- get_disgenet_api_key(
  email = "", 
  password = "" )

Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

path.list<- unique(path.enrich$PathwayID)
dat.out<- NULL

path.list<- unique(path.enrich$PathwayID)
path.seq<- path.list[seq(100, length(path.list), 100)]
dat.out<- NULL

#for each pathway, ID which disease overlaps
#Requiring a minimum of 3 genes overlap and FDR < 0.05
#Only keeping data on nervous system diseases

for(path.in in path.list){
  
  in.dat<- path.enrich[PathwayID %in% path.in] #collect relevant pathway data
  
  tryCatch(
    {
      #Determine disease enrichment
      val<- disease_enrichment(entities =unique(in.dat$GeneID), warnings=FALSE,
                               vocabulary = "ENTREZ", database = "ALL")
      val<- val@qresult
      val<- as.data.table(separate_rows(val, c(shared_geneid, shared_symbol), sep=";"))
      val<- val[Description %in% c("Alzheimer's Disease", "Parkinson Disease", "Brain Neoplasms", 
                                   "Epilepsy", "Migraine Disorders", "Multiple Sclerosis")]
      
      val$PathwayID<- unique(in.dat$PathwayID)
      val$PathwayName<- unique(in.dat$PathwayName)
      
      #Reduce to associations meeting the minimum significance requirements
      val<- val[FDR< 0.05]
      val<- val[Count >=3]
      
      suppressWarnings()
    },
    error = function(e) {
      #print(paste("error", path.in, which(path.list %in% path.in)))
      val<- NULL
    }
  )
  
  
  dat.out<- rbind(dat.out, val)
  
  if(path.in %in% path.seq){
    print(which(path.list %in% path.in))
  }
  
}

###NOTE: SAVE FILE FOR FUTURE STEPS

setwd("")
fwrite(path.out, file="Generated_Data/Full_Enriched_Path_Disease_Links.csv")

############################################################################
###SI Figure 1

setwd("")
path.out<- fread("Generated_Data/Full_Enriched_Path_Disease_Links.csv")

#Plot differences in number of data values based on different restrictions
#p values .05, .01, .001, .0001
#gene overlap = 3, 5, 10, 20, 35, 50
#Prepare different datasets and put together

#for each combination of p-value and gene count, determine relevant data values
p.count<- function(p.value, g.value, p.print){
  
  #restrict dataset
  int.val<- path.out[FDR < p.value & Count >=g.value, list(Disease=Description, PathwayID, PathwayName, GeneID=shared_geneid, 
                                                           GeneSymbol=shared_symbol, FDR, GeneCount=Count)]
  #counts of unique items
  int.un.count<- int.val[, .(genes=length(unique(GeneID)), pathways=length(unique(PathwayID)))]
  
  #linkage counts
  pd.link<- int.val[, .(.N), by=.(PathwayID, Disease)]
  pgd.link<- int.val[, .(.N), by=.(PathwayID, GeneID, Disease)]
  
  int.un.count<- cbind(int.un.count, data.table(`path-dis`=nrow(pd.link)), 
                       data.table(`path-gene-dis`=nrow(pgd.link)))
  
  #Put values together
  dat.out<- melt(int.un.count)
  dat.out$variable<- paste("Unique", dat.out$variable)
  
  dat.out$Pvalue<- p.print
  dat.out$GeneCount<- g.value
  
  dat.out
}

range(path.out$Count)
range(path.out$FDR)

path.plot<- NULL
p.list<- c(0.0001, 0.001, 0.01, 0.05)
g.list<- c(3, 5, 10, 20, 35, 50)

for(p.in in p.list){
  for(g.in in g.list){
    
    #For each pvalue - gene overlap combination, compute unique paths, genes, and linkages
    path.plot.int<- p.count(p.in, g.in, paste("ORA p <", p.in))
    path.plot<- rbind(path.plot, path.plot.int)
    
  }
  
}

path.plot$Pvalue<- factor(path.plot$Pvalue, levels= unique(path.plot$Pvalue), ordered=TRUE)


setwd("")
jpeg("SITables_Figures/FigureS1.jpeg", width=6.5, height=6, units="in", res=300)
ggplot(path.plot) + geom_point(aes(x=GeneCount, y=value, color=GeneCount)) + theme_bw() +
  facet_grid(rows=vars(variable), cols=vars(Pvalue), scales="free_y") + theme(panel.grid = element_blank())
graphics.off()

setwd("")
fwrite(path.plot, "SITables_Figures/SI_FigureS1.csv")
