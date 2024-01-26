############################################################################

### Determine pesticides used per state and disease rates per state
###
### Author: Marissa Kosnik

############################################################################

rm(list = ls())

require(data.table)
require(tidyverse)
require (plyr)
require(readxl)

############################################################################
###USGS pesticide application data over time

###NOTE: DOWNLOAD USGS DATA AS DESCRIBED

#From a single folder of USGS pesticide application data per year 1992-2019
#read in files 

setwd("/USGS")
temp<- list.files()

#Read in all results files and join 

dat.out<- NULL
for(file.in in temp){
  
  val<- fread(file.in)
  val$STATE_FIPS_CODE <- as.character(val$STATE_FIPS_CODE)
  
  #Prep annual estimate as median low/high value per county
  val<- val[, .(KG_APPLIED=median(na.omit(c(EPEST_HIGH_KG, EPEST_LOW_KG)))), 
            by=.(YEAR, COMPOUND, STATE_FIPS_CODE, COUNTY_FIPS_CODE)]
  
  #Join by state - sum KG applied over whole state
  val<- val[, .(KG_APPLIED=sum(KG_APPLIED)), by=.(YEAR, COMPOUND, STATE_FIPS_CODE)]
  
  dat.out<- rbind(dat.out, val)
  print(file.in)
}


############################################################################
#Join to state names

#Make sure state codes match typical FIPS format
dat.out$STATE_FIPS <- as.character(dat.out$STATE_FIPS_CODE)
dat.out$STATE_FIPS<- vapply(dat.out$STATE_FIPS, function(x) ifelse(nchar(x) == 1, paste0("0", x), x), 
                            FUN.VALUE = character(1))
length(unique(dat.out$STATE_FIPS))

#Add state names
require(usa)

state.dat<- as.data.table(usa::states)
dat.out$State<- state.dat[match(dat.out$STATE_FIPS, state.dat$fips)]$name

#Check all years
test<- dat.out[, .(Chems=length(unique(COMPOUND)), KG=sum(KG_APPLIED)), by=.(YEAR)]
test
#2019 has too few chems (>200 fewer other years), remove
dat.out<- dat.out[YEAR != 2019]

#Now export list of chemical names to get matching CAS data from CompTox dashboard

setwd("")
fwrite(as.data.table(unique(dat.out$COMPOUND)), file="USGS_Pesticides.txt", col.names=FALSE)
#Run on CompTox Dashboard to determine CASRNs
#https://comptox.epa.gov/dashboard/

###NOTE: DOWNLOAD COMPTOX DASHBOARD DATA AS DESCRIBED

setwd("")
pest.cas<- fread("CCD-Batch-Search_....csv")
#Input (pesticide names in USGS) used to determine CAS

dat.out$CASRN<- pest.cas[match(dat.out$COMPOUND, pest.cas$INPUT)]$CASRN
length(unique(dat.out$CASRN))
length(unique(dat.out$COMPOUND))

test<- dat.out[, .(.N), by=.(CASRN, COMPOUND)]
length(unique(dat.out$COMPOUND)) == (length(unique(dat.out$CASRN)) + nrow(test[CASRN=="N/A"]))
#Where available, all CASRNs are unique
dat.out[CASRN=="N/A"]$CASRN<- NA

usgs.year.pest<- dat.out[, list(YEAR, COMPOUND, CASRN, STATE_FIPS, State, KG_APPLIED)]

############################################################################
#Global burden of disease data over time

###NOTE: DOWNLOAD GBD DATA AS DESCRIBED

#https://www.healthdata.org/research-analysis/gbd
setwd("")
gbd.dat<- fread("IHME-GBD_2019_DATA-....csv")
#From a file with GBD data for 1992 - 2019, 6 cause_names (listed below)
#incidence and prevalence rates per US state
#And ages = all ages, 55+ years, and age-standardized

#read in file and prep 

#Convert GBD disease names to names that will match in DisGeNET
gbd.dat[cause_name %in% "Parkinson's disease"]$cause_name<- "Parkinson Disease"
gbd.dat[cause_name %in% "Idiopathic epilepsy"]$cause_name<- "Epilepsy"
gbd.dat[cause_name %in% "Multiple sclerosis"]$cause_name<- "Multiple Sclerosis"
gbd.dat[cause_name %in% "Alzheimer's disease and other dementias"]$cause_name<- "Alzheimer's Disease"
gbd.dat[cause_name %in% "Brain and central nervous system cancer"]$cause_name<- "Brain Neoplasms"
gbd.dat[cause_name %in% "Migraine"]$cause_name<- "Migraine Disorders"
gbd.dat$Disease<- gbd.dat$cause_name

length(unique(gbd.dat$location_name))
gbd.dat<- gbd.dat[!measure_name %in% "Deaths" & metric_name %in% "Rate"] #Not using these data
gbd.dat<- gbd.dat[location_name %in% usgs.year.pest$State & age_name %in% "All ages"]
length(unique(gbd.dat$location_name))

rm(list=setdiff(ls(), c("usgs.year.pest", "gbd.dat")))

###NOTE: SAVE WORKSPACE FOR FUTURE STEPS

setwd("")
save.image(file="Workspaces/Step1_USGS_GBD_Data.RData")
