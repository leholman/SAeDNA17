#############################################
####==== South African eDNA Analysis ====####
####==== Luke E. Holman====30.03.2020====####
#############################################

###Script 0 - QC & Prep Script###

####====0.0 Packages====####

library("dplyr")
library("vegan")
library("reshape")
library("reshape2")
library("ggplot2")
library("stringr")
library("tidyr")
library("RColorBrewer")
library("readxl")
library("maps")
library("mapdata")
require(ggplot2)
require(ggmap)
require(mapproj)
require(rgeos)
require(maptools)
require(sp)
require(raster)
require(rgdal)
require(dismo)
require(Biostrings)
library("RgoogleMaps")
library("reshape")
library("metabarTOAD")
library("dada2")
library("iNEXT")
library("pairwiseAdonis")
library("betapart")
library("Hmisc")
library(UpSetR)


####====0.1 Parameters & Reading In====####

#Set some variables 
minreads <- 3
items <- NULL

#Set the seed 
set.seed("123456")

#Set wd and load up sediment data 
metadat <- read.csv(file="metadata/locations.csv")
latlongdata <- read.csv("metadata/mapdata.csv")
latlongdata$PERMori <- factor(latlongdata$PERMori,levels = c("W", "S", "E"))
tempdat <- read.csv("metadata/temp.csv")
envdat <- read.csv("metadata/EnvDat.csv")


#Set A plalette
palette(c('#51B9E0','#4bad84','#E8A016'))
#platette two 
#palette(c('#D36526','#2671B4'))

####====0.2 Data QC====####

files <- list.files("rawdata",full.names = T)
for (file in  files){
  
  rawdat <-read.csv(file=file,row.names = 1)
  
  
  #Seperate controls and samples
  samples <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="sample",gsub("(.*)_.*","\\1",basename(file))]]
  controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control",gsub("(.*)_.*","\\1",basename(file))]]
  
  
  #Filter 1 - minimum number of reads for any ID
  samples[samples< minreads] <- 0
  samples <- samples[rowSums(samples) > 0,]
  
  #Filter 2 - within samples OTU must appear in more than one sample (this works becuase there are lots of reps per site and sample)
  filtersam <- samples
  filtersam[filtersam>0 ] <- 1
  filtersam <-filtersam[rowSums(filtersam) > 1,]
  samples <- samples[rownames(samples) %in% rownames(filtersam),]
  
  #Filter 3 -Make the maximum umber of reads for each OTU in the contam the zero value in the main data
  controlsCONTAM <- controls[rowSums(controls) > 0,]
  for (contamOTU in 1:length(controlsCONTAM[,1])){
    loopOTU <- row.names(controlsCONTAM[contamOTU,])
    loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
    if (any(is.na(samples[loopOTU,]))){next}
    samples[loopOTU,samples[loopOTU,]<loopMax] <- 0
    print(paste("Cleaning contaminants",contamOTU))
  }
  
  
  
  
  #Now we write out the cleaned data
  newname <- paste("cleaned/","Cleaned.",basename(file),sep="")
  write.csv(samples,file=newname)
  
}

# Next we create a dataset that only contain OTUs represented in all 3 technical replicates #
##COI
COI <- read.csv(file = "cleaned/Cleaned.COI_DADA.csv")
rownames(COI) <- COI$X
COI <- COI[,-1]


#This little loop gets rid of any records that do not occur in every technical replicate 
for (site in unique(substr(colnames(COI),3,4))){
  loopdat <- COI[,substr(colnames(COI),3,4) == site]
  loopdat[loopdat>0] <- 1
  COI[!rowSums(loopdat)==3,substr(colnames(COI),3,4) == site] <- 0
}

## we rarefy the data
rCOI<-t(rrarefy(t(COI[rowSums(COI)>0,]),min(colSums(COI[rowSums(COI)>0,]))))
rCOI <- as.data.frame(rCOI)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
cCOI <- matrix(ncol=length(unique(substr(colnames(rCOI),3,4))),nrow = length(rCOI[,1]))
colnames(cCOI) <- unique(substr(colnames(rCOI),3,4))
rownames(cCOI) <- rownames(rCOI)
for (site in unique(substr(colnames(rCOI),3,4))){
  cCOI[,site] <- round(rowMeans(rCOI[,substr(colnames(rCOI),3,4) == site]))
}
rCOI <- as.data.frame(cCOI)
write.csv(rCOI,"cleaned/rarefied.COI.csv")


##18S
zhan <- read.csv(file = "cleaned/Cleaned.ZHAN_DADA.csv")
rownames(zhan) <- zhan$X
zhan <- zhan[,-1]


#This little loop gets rid of any records that do not occur in every technical replicate 
for (site in unique(substr(colnames(zhan),3,4))){
  loopdat <- zhan[,substr(colnames(zhan),3,4) == site]
  loopdat[loopdat>0] <- 1
  zhan[!rowSums(loopdat)==3,substr(colnames(zhan),3,4) == site] <- 0
}
## we rarefy the data
r18S<-t(rrarefy(t(zhan[rowSums(zhan)>0,]),min(colSums(zhan[rowSums(zhan)>0,]))))
r18S <- as.data.frame(r18S)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
c18S <- matrix(ncol=length(unique(substr(colnames(r18S),3,4))),nrow = length(r18S[,1]))
colnames(c18S) <- unique(substr(colnames(r18S),3,4))
rownames(c18S) <- rownames(r18S)
for (site in unique(substr(colnames(r18S),3,4))){
  c18S[,site] <- round(rowMeans(r18S[,substr(colnames(r18S),3,4) == site]))
}
r18S <- as.data.frame(c18S)
write.csv(r18S,"cleaned/rarefied.18S.csv")

##ProK
ProK <- read.csv(file = "cleaned/Cleaned.PROK_DADA.csv")
rownames(ProK) <- ProK$X
ProK <- ProK[,-1]


#This little loop gets rid of any records that do not occur in every technical replicate 
for (site in unique(substr(colnames(ProK),1,2))){
  loopdat <- ProK[,substr(colnames(ProK),1,2) == site]
  loopdat[loopdat>0] <- 1
  ProK[!rowSums(loopdat)==3,substr(colnames(ProK),1,2) == site] <- 0
}

## we rarefy the data
rProK<-t(rrarefy(t(ProK[rowSums(ProK)>0,]),min(colSums(ProK[rowSums(ProK)>0,]))))
rProK <- as.data.frame(rProK)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
cProK <- matrix(ncol=length(unique(substr(colnames(rProK),1,2))),nrow = length(rProK[,1]))
colnames(cProK) <- unique(substr(colnames(rProK),1,2))
rownames(cProK) <- rownames(rProK)
for (site in unique(substr(colnames(rProK),1,2))){
  cProK[,site] <- round(rowMeans(rProK[,substr(colnames(rProK),1,2) == site]))
}
rProK <- as.data.frame(cProK)

write.csv(rProK,"cleaned/rarefied.ProK.csv")


####====0.3. Taxonomy & Final QC====####

# Process the taxonomy 
#COI
COIassignments <- ParseTaxonomy(blastoutput = "Taxonomy/COI.dada2.blast.270819.txt.gz",lineages = "Taxonomy/lineages-2019-08-14.csv.gz")
#these are then put in the same order as the ecological data for later analyses
COIassignments <- COIassignments[match(rownames(rCOI),as.character(COIassignments$OTU)),]

#18S
z18Sassignments <- ParseTaxonomy(blastoutput = "Taxonomy/18S.blast.results.txt.gz",lineages = "Taxonomy/lineages-2019-08-14.csv.gz", pctThreshold = 99)
#as previous
z18Sassignments <- z18Sassignments[match(rownames(r18S),as.character(z18Sassignments$OTU)),]

#16S -testing the dada2 assignment RDP
ProKOTUS <- readDNAStringSet("/Volumes/BackUp2/Ext.HT.Sequence.data/2018.SAeDNA/SAeDNA3.1/5.OTUs/dada2.fasta")
ProKtaxaALL <- assignTaxonomy(ProKOTUS, "~/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE,verbose = TRUE)
ProKtaxaALL$OTU <- names(ProKOTUS)
ProKOTUS <- ProKOTUS[ProKOTUS@ranges@NAMES %in% rownames(rProK)]
ProKtaxa <- assignTaxonomy(ProKOTUS, "~/Downloads/silva_nr_v132_train_set.fa.gz", multithread=TRUE,verbose = TRUE)
ProKspp <- assignSpecies(ProKOTUS, "~/Downloads/silva_species_assignment_v132.fa.gz",verbose = TRUE)
ProKassignments <- ParseTaxonomy(blastoutput = "Taxonomy/ProK.dada2.22082019.txt.gz",lineages = "Taxonomy/lineages-2019-08-14.csv.gz", pctThreshold = 99)
ProKassignments <- ProKassignments[match(rownames(rProK),as.character(ProKassignments$OTU)),]

ProKtaxa <- as.data.frame(ProKtaxa)
ProKtaxa$OTUname <- ProKOTUS@ranges@NAMES

table(as.character(ProKtaxa[match(ProKtaxa$OTUname,row.names(rProK)),4]))

write.csv(COIassignments,"Taxonomy/CleanedTaxonomy/COI.blast.csv")
write.csv(z18Sassignments,"Taxonomy/CleanedTaxonomy/18S.blast.csv")
write.csv(ProKassignments,"Taxonomy/CleanedTaxonomy/ProK.blast.csv")
write.csv(ProKtaxa,"Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.csv")
write.csv(ProKtaxaALL,"Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.ALL.csv")

#Use these to avoid processing data again if the above has already been run
#COIassignments <- read.csv("Taxonomy/CleanedTaxonomy/COI.blast.csv")
#z18Sassignments <- read.csv("Taxonomy/CleanedTaxonomy/18S.blast.csv")
#ProKtaxa <- read.csv("Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.csv")
#ProKtaxaALL <- read.csv("Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.ALL.csv")

##Now we collapse OTUs assigned to the same taxa
#COI
sCOIassignments<-COIassignments[as.character(COIassignments$OTU) %in% rownames(rCOI),]
sCOIassignments <- sCOIassignments[sCOIassignments$assignmentQual=="High", ]

#Write out the COI data for checking against WRIMS
write.csv(sCOIassignments,file="Taxonomy/NNS/COItobechecked.csv")

for (name in as.character(names(table(sCOIassignments$species)[table(sCOIassignments$species)>1]))){
  collapseOTUs <- as.character(sCOIassignments$OTU[sCOIassignments$species==name]) 
  MotherOTU <- names(sort(rowSums(rCOI[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  rCOI[MotherOTU,] <- rCOI[MotherOTU,] + colSums(rCOI[collapseOTUs,])
  rCOI <- rCOI[-match(collapseOTUs,rownames(rCOI)),]
}

#18S
sz18Sassignments <- z18Sassignments[as.character(z18Sassignments$OTU) %in% rownames(r18S),]
sz18Sassignments <- sz18Sassignments[sz18Sassignments$assignmentQual=="High", ]

for (name in as.character(names(table(sz18Sassignments$species)[table(sz18Sassignments$species)>1]))){
  collapseOTUs <- as.character(sz18Sassignments$OTU[sz18Sassignments$species==name]) 
  MotherOTU <- names(sort(rowSums(r18S[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  r18S[MotherOTU,] <- r18S[MotherOTU,] + colSums(r18S[collapseOTUs,])
  r18S <- r18S[-match(collapseOTUs,rownames(r18S)),]
}

####====0.3 Control Data====####
#COI
file <- "rawdata/COI_DADA.csv"
rawdat <- read.csv(file=file,row.names = 1)
controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control",gsub("(.*)_.*","\\1",basename(file))]]
controls <- controls[rowSums(controls)>0,]
COIcontrols <- data.frame(matrix(nrow= dim(controls)[1],ncol=11))
rownames(COIcontrols) <- rownames(controls)
colnames(COIcontrols) <- colnames(COIassignments)
temp <- as.data.frame(apply(COIassignments,2, as.character),stringsAsFactors=FALSE)
for (row in 1:dim(controls)[1]){
  if(!rownames(controls)[row] %in% COIassignments$OTU){
    next()}
  COIcontrols[rownames(COIcontrols)[row],] <-  temp[na.omit(match(rownames(COIcontrols)[row],as.character(temp$OTU))),]
}

COIcontrols <- cbind(controls,COIcontrols)
#Stats for COI
mean(colSums(COIcontrols[1: 1:(dim(controls)[2])]))
sd(colSums(COIcontrols[1: 1:(dim(controls)[2])]))
table(COIcontrols$assignmentQual)
temp <- COIcontrols[1: 1:(dim(controls)[2])]
temp[temp>1] <- 1
mean(colSums(temp))
sd(colSums(temp))



write.csv(COIcontrols,file="controls/COIcontrol.csv")

#18S
file <- "rawdata/ZHAN_DADA.csv"

rawdat <- read.csv(file=file,row.names = 1)
controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control",gsub("(.*)_.*","\\1",basename(file))]]
controls <- controls[rowSums(controls)>0,]
z18Scontrols <- data.frame(matrix(nrow= dim(controls)[1],ncol=11))
rownames(z18Scontrols) <- rownames(controls)
colnames(z18Scontrols) <- colnames(z18Sassignments)
temp <- as.data.frame(apply(z18Sassignments,2, as.character),stringsAsFactors=FALSE)
for (row in 1:dim(controls)[1]){
  if(!rownames(controls)[row] %in% z18Sassignments$OTU){
    next()}
  z18Scontrols[rownames(z18Scontrols)[row],] <-  temp[na.omit(match(rownames(z18Scontrols)[row],as.character(temp$OTU))),]
}

z18Scontrols <- cbind(controls,z18Scontrols)
#Stats for 18S
mean(colSums(z18Scontrols[1: 1:(dim(controls)[2])]))
sd(colSums(z18Scontrols[1: 1:(dim(controls)[2])]))
temp <- z18Scontrols[1: 1:(dim(controls)[2])]
temp[temp>1] <- 1
mean(colSums(temp))
sd(colSums(temp))


table(z18Scontrols$assignmentQual)

write.csv(z18Scontrols,file="controls/18Scontrol.csv")



#16S
file <- "rawdata/PROK_DADA.csv"
rawdat <- read.csv(file=file,row.names=1)
controls <- rawdat[colnames(rawdat) %in% metadat[metadat$Type=="control",gsub("(.*)_.*","\\1",basename(file))]]
controls <- controls[rowSums(controls)>0,]
ProKcontrols <- data.frame(matrix(nrow= dim(controls)[1],ncol=8))
rownames(ProKcontrols) <- rownames(controls)
colnames(ProKcontrols) <- colnames(ProKtaxaALL)
temp <- as.data.frame(apply(ProKtaxaALL,2, as.character),stringsAsFactors=FALSE)
for (row in 1:dim(controls)[1]){
  if(!rownames(controls)[row] %in% ProKtaxaALL$OTU){
    next()}
  ProKcontrols[rownames(ProKcontrols)[row],] <-  temp[na.omit(match(rownames(ProKcontrols)[row],as.character(temp$OTU))),]
}

ProKcontrols <- cbind(controls,ProKcontrols)
#Stats for 16S
mean(colSums(ProKcontrols[1: 1:(dim(controls)[2])]))
sd(colSums(ProKcontrols[1: 1:(dim(controls)[2])]))
table(ProKcontrols$assignmentQual)
temp <- ProKcontrols[1: 1:(dim(controls)[2])]
temp[temp>1] <- 1
mean(colSums(temp))
sd(colSums(temp))

#OTU19 contribution
mean(as.numeric(ProKcontrols[match("OTU_19",row.names(ProKcontrols)),2:6]))

#Stats without OTU 19 

ProKcontrolsTrunc <- ProKcontrols[-match("OTU_19",row.names(ProKcontrols)),]
#Stats for 16S
mean(colSums(ProKcontrolsTrunc[1: 1:(dim(controls)[2])]))
sd(colSums(ProKcontrolsTrunc[1: 1:(dim(controls)[2])]))

temp <- ProKcontrolsTrunc[1: 1:(dim(controls)[2])]
temp[temp>1] <- 1
mean(colSums(temp))
sd(colSums(temp))



write.csv(ProKcontrols,file="controls/ProKcontrol.csv")

##How many OTUs assigned taxonomy?
table(COIassignments$assignmentQual[match(rownames(rCOI),COIassignments$OTU)],useNA="always")

table(z18Sassignments$assignmentQual[match(rownames(r18S),z18Sassignments$OTU)],useNA="always")

table(ProKtaxa$Genus,useNA="always")
