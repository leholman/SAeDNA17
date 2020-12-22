#############################################
####==== South African eDNA Analysis ====####
####==== Luke E. Holman====16.10.2020====####
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
library(seqinr)


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
files <- files[grep(".csv",files)]

for (file in  files){
  
  rawdat <-read.csv(file=file,row.names = 1)
  
  
  #Separate controls and samples
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

# Process the blast taxonomy 
#COI
#NOTE - for upload the COI file has been split, unzip concatenate and rezip beforre running the below line
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
write.csv(rCOI,"cleaned/rarefied.COI.csv")


#18S we dont do this becuase 18S gives genus at best for this gene fragment
#sz18Sassignments <- z18Sassignments[as.character(z18Sassignments$OTU) %in% rownames(r18S),]
#sz18Sassignments <- sz18Sassignments[sz18Sassignments$assignmentQual=="High", ]

#for (name in as.character(names(table(sz18Sassignments$species)[table(sz18Sassignments$species)>1]))){
#  collapseOTUs <- as.character(sz18Sassignments$OTU[sz18Sassignments$species==name]) 
#  MotherOTU <- names(sort(rowSums(r18S[collapseOTUs,]),decreasing = TRUE))[1]
#  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
#  r18S[MotherOTU,] <- r18S[MotherOTU,] + colSums(r18S[collapseOTUs,])
#  r18S <- r18S[-match(collapseOTUs,rownames(r18S)),]
#}


#RDP data
RDP.COI.assign <- read.table("Taxonomy/RDP.assigns/RDP.class.COI.v4.txt",sep="\t")
RDP.18S.assign <- read.delim("Taxonomy/RDP.assigns/RDP.class.18S.v3.2.txt",sep="\t",header = F)

####====0.4 Creating Data Subsets====####

####Initial reviews suggested reworking the data into taxonomic rather than marker based subsets 
## First let's try doing that at domain level
#Eukaryote data
#COI

##COI
COI <- read.csv(file = "cleaned/Cleaned.COI_DADA.csv",row.names = 1)

#Truncate to include records across all three tech reps
for (site in unique(substr(colnames(COI),3,4))){
  loopdat <- COI[,substr(colnames(COI),3,4) == site]
  loopdat[loopdat>0] <- 1
  COI[!rowSums(loopdat)==3,substr(colnames(COI),3,4) == site] <- 0
}

##Subset by domain 

COI.euk <- COI[rownames(COI) %in% RDP.COI.assign$V1[RDP.COI.assign$V6=="Eukaryota" & RDP.COI.assign$V8>0.3],]


## we rarefy the data
rCOI.euk<-t(rrarefy(t(COI.euk[rowSums(COI.euk)>0,]),min(colSums(COI.euk[rowSums(COI.euk)>0,]))))
rCOI.euk <- as.data.frame(rCOI.euk)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
cCOI.euk <- matrix(ncol=length(unique(substr(colnames(rCOI.euk),3,4))),nrow = length(rCOI.euk[,1]))
colnames(cCOI.euk) <- unique(substr(colnames(rCOI.euk),3,4))
rownames(cCOI.euk) <- rownames(rCOI.euk)
for (site in unique(substr(colnames(rCOI.euk),3,4))){
  cCOI.euk[,site] <- round(rowMeans(rCOI.euk[,substr(colnames(rCOI.euk),3,4) == site]))
}
rCOI.euk <- as.data.frame(cCOI.euk)

#COI
sCOIassignments<-COIassignments[as.character(COIassignments$OTU) %in% rownames(rCOI.euk),]
sCOIassignments <- sCOIassignments[sCOIassignments$assignmentQual=="High", ]

for (name in as.character(names(table(sCOIassignments$species)[table(sCOIassignments$species)>1]))){
  collapseOTUs <- as.character(sCOIassignments$OTU[sCOIassignments$species==name]) 
  MotherOTU <- names(sort(rowSums(rCOI.euk[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  rCOI.euk[MotherOTU,] <- rCOI.euk[MotherOTU,] + colSums(rCOI.euk[collapseOTUs,])
  rCOI.euk <- rCOI.euk[-match(collapseOTUs,rownames(rCOI.euk)),]
}

write.csv(rCOI.euk,"cleaned/rarefied.COI.euk.csv")



#18S
zhan <- read.csv(file = "cleaned/Cleaned.ZHAN_DADA.csv",row.names = 1)

#Truncate to include records across all three tech reps
for (site in unique(substr(colnames(zhan),3,4))){
  loopdat <- zhan[,substr(colnames(zhan),3,4) == site]
  loopdat[loopdat>0] <- 1
  zhan[!rowSums(loopdat)==3,substr(colnames(zhan),3,4) == site] <- 0
}

##Subset by domain 

zhan.euk <- zhan[rownames(zhan) %in% RDP.18S.assign$V1[RDP.18S.assign$V6=="Eukaryota" & RDP.18S.assign$V8>0.3],]

## we rarefy the data
r18S.euk<-t(rrarefy(t(zhan.euk[rowSums(zhan.euk)>0,]),min(colSums(zhan.euk[rowSums(zhan.euk)>0,]))))
r18S.euk <- as.data.frame(r18S.euk)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
c18S.euk <- matrix(ncol=length(unique(substr(colnames(r18S.euk),3,4))),nrow = length(r18S.euk[,1]))
colnames(c18S.euk) <- unique(substr(colnames(r18S.euk),3,4))
rownames(c18S.euk) <- rownames(r18S.euk)
for (site in unique(substr(colnames(r18S.euk),3,4))){
  c18S.euk[,site] <- round(rowMeans(r18S.euk[,substr(colnames(r18S.euk),3,4) == site]))
}
r18S.euk <- as.data.frame(c18S.euk)
write.csv(r18S.euk,"cleaned/rarefied.18S.euk.csv")


#Prokaryote data
#16S is all prokaryotic apart form the below ASV

#ProKtaxa[2723,]

#this BLASTs as a prokaryotic spp so it is an error

##Now let's try some sub domain level subsets 

#Metazoans

##COI
COI <- read.csv(file = "cleaned/Cleaned.COI_DADA.csv",row.names = 1)

#Truncate to include records across all three tech reps
for (site in unique(substr(colnames(COI),3,4))){
  loopdat <- COI[,substr(colnames(COI),3,4) == site]
  loopdat[loopdat>0] <- 1
  COI[!rowSums(loopdat)==3,substr(colnames(COI),3,4) == site] <- 0
}

##Subset by kingdom

COI.met <- COI[rownames(COI) %in% RDP.COI.assign$V1[RDP.COI.assign$V9=="Metazoa" & RDP.COI.assign$V11>0.3],]

## we rarefy the data
rCOI.met<-t(rrarefy(t(COI.met[rowSums(COI.met)>0,]),min(colSums(COI.met[rowSums(COI.met)>0,]))))
rCOI.met <- as.data.frame(rCOI.met)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
cCOI.met <- matrix(ncol=length(unique(substr(colnames(rCOI.met),3,4))),nrow = length(rCOI.met[,1]))
colnames(cCOI.met) <- unique(substr(colnames(rCOI.met),3,4))
rownames(cCOI.met) <- rownames(rCOI.met)
for (site in unique(substr(colnames(rCOI.met),3,4))){
  cCOI.met[,site] <- round(rowMeans(rCOI.met[,substr(colnames(rCOI.met),3,4) == site]))
}
rCOI.met <- as.data.frame(cCOI.met)

#COI
sCOIassignments<-COIassignments[as.character(COIassignments$OTU) %in% rownames(rCOI.met),]
sCOIassignments <- sCOIassignments[sCOIassignments$assignmentQual=="High", ]

for (name in as.character(names(table(sCOIassignments$species)[table(sCOIassignments$species)>1]))){
  collapseOTUs <- as.character(sCOIassignments$OTU[sCOIassignments$species==name]) 
  MotherOTU <- names(sort(rowSums(rCOI.met[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  rCOI.met[MotherOTU,] <- rCOI.met[MotherOTU,] + colSums(rCOI.met[collapseOTUs,])
  rCOI.met <- rCOI.met[-match(collapseOTUs,rownames(rCOI.met)),]
}

write.csv(rCOI.met,"cleaned/rarefied.COI.met.csv")

##Now the metazoans in 18S 

z18S <- read.csv(file = "cleaned/Cleaned.ZHAN_DADA.csv",row.names = 1)

#Truncate to include records across all three tech reps
for (site in unique(substr(colnames(z18S),3,4))){
  loopdat <- z18S[,substr(colnames(z18S),3,4) == site]
  loopdat[loopdat>0] <- 1
  z18S[!rowSums(loopdat)==3,substr(colnames(z18S),3,4) == site] <- 0
}

##Subset by kingdom

z18S.met <- z18S[rownames(z18S) %in% RDP.18S.assign$V1[RDP.18S.assign$V9=="Metazoa_Animalia" & RDP.18S.assign$V11>0.3],]

## we rarefy the data
r18S.met<-t(rrarefy(t(z18S.met[rowSums(z18S.met)>0,]),min(colSums(z18S.met[rowSums(z18S.met)>0,]))))
r18S.met <- as.data.frame(r18S.met)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
c18S.met <- matrix(ncol=length(unique(substr(colnames(r18S.met),3,4))),nrow = length(r18S.met[,1]))
colnames(c18S.met) <- unique(substr(colnames(r18S.met),3,4))
rownames(c18S.met) <- rownames(r18S.met)
for (site in unique(substr(colnames(r18S.met),3,4))){
  c18S.met[,site] <- round(rowMeans(r18S.met[,substr(colnames(r18S.met),3,4) == site]))
}
r18S.met <- as.data.frame(c18S.met)

write.csv(r18S.met,"cleaned/rarefied.18S.met.csv")


##Defining protista COI

taxaData <- RDP.COI.assign[match(rownames(rCOI.euk),RDP.COI.assign$V1),]
COIprotists <- read.csv("Taxonomy/ProtistDesignations/COI.protist.classificiations.csv")

taxaData2 <- taxaData[taxaData$V12 %in% COIprotists$ID[COIprotists$Level=="Phylum"] & taxaData$V14>0.3,]
taxaData3 <- taxaData[taxaData$V15 %in% COIprotists$ID[COIprotists$Level=="Class"] & taxaData$V17>0.3,]
taxaData4 <- taxaData[taxaData$V18 %in% COIprotists$ID[COIprotists$Level=="Order"] & taxaData$V20>0.3,]

COI.protist.taxa <- rbind(taxaData2,taxaData3,taxaData4)

COI <- read.csv(file = "cleaned/Cleaned.COI_DADA.csv",row.names = 1)

#Truncate to include records across all three tech reps
for (site in unique(substr(colnames(COI),3,4))){
  loopdat <- COI[,substr(colnames(COI),3,4) == site]
  loopdat[loopdat>0] <- 1
  COI[!rowSums(loopdat)==3,substr(colnames(COI),3,4) == site] <- 0
}

##Subset by protists 

COI.pts <- COI[rownames(COI) %in% COI.protist.taxa$V1,]

## we rarefy the data
rCOI.pts<-t(rrarefy(t(COI.pts[rowSums(COI.pts)>0,]),min(colSums(COI.pts[rowSums(COI.pts)>0,]))))
rCOI.pts <- as.data.frame(rCOI.pts)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
cCOI.pts <- matrix(ncol=length(unique(substr(colnames(rCOI.pts),3,4))),nrow = length(rCOI.pts[,1]))
colnames(cCOI.pts) <- unique(substr(colnames(rCOI.pts),3,4))
rownames(cCOI.pts) <- rownames(rCOI.pts)
for (site in unique(substr(colnames(rCOI.pts),3,4))){
  cCOI.pts[,site] <- round(rowMeans(rCOI.pts[,substr(colnames(rCOI.pts),3,4) == site]))
}
rCOI.pts <- as.data.frame(cCOI.pts)

#COI
sCOIassignments<-COIassignments[as.character(COIassignments$OTU) %in% rownames(rCOI.pts),]
sCOIassignments <- sCOIassignments[sCOIassignments$assignmentQual=="High", ]

for (name in as.character(names(table(sCOIassignments$species)[table(sCOIassignments$species)>1]))){
  collapseOTUs <- as.character(sCOIassignments$OTU[sCOIassignments$species==name]) 
  MotherOTU <- names(sort(rowSums(rCOI.pts[collapseOTUs,]),decreasing = TRUE))[1]
  collapseOTUs <- collapseOTUs[-grep(MotherOTU,collapseOTUs)]
  rCOI.pts[MotherOTU,] <- rCOI.pts[MotherOTU,] + colSums(rCOI.pts[collapseOTUs,])
  rCOI.pts <- rCOI.pts[-match(collapseOTUs,rownames(rCOI.pts)),]
}

write.csv(rCOI.pts,"cleaned/rarefied.COI.pts.csv")


##Defining protista 18S

r18S.euk <- read.csv("cleaned/rarefied.18S.euk.csv",row.names = 1)
taxaData <- RDP.18S.assign[match(rownames(r18S.euk),RDP.18S.assign$V1),]

z18Sprotists <- read.csv("Taxonomy/ProtistDesignations/18S.protist.classificiations.csv")

taxaData2 <- taxaData[taxaData$V9 %in% z18Sprotists$ID[z18Sprotists$Level=="Kingdom.ish"] & taxaData$V14>0.3,]
taxaData3 <- taxaData[taxaData$V12 %in% z18Sprotists$ID[z18Sprotists$Level=="Phylum"] & taxaData$V17>0.3,]
taxaData4 <- taxaData[taxaData$V15 %in% z18Sprotists$ID[z18Sprotists$Level=="Class"] & taxaData$V20>0.3,]

z18S.protist.taxa <- rbind(taxaData2,taxaData3,taxaData4)


z18S <- read.csv(file = "cleaned/Cleaned.ZHAN_DADA.csv",row.names = 1)

#Truncate to include records across all three tech reps
for (site in unique(substr(colnames(z18S),3,4))){
  loopdat <- z18S[,substr(colnames(z18S),3,4) == site]
  loopdat[loopdat>0] <- 1
  z18S[!rowSums(loopdat)==3,substr(colnames(z18S),3,4) == site] <- 0
}

##Subset by protists

z18S.pts <- z18S[rownames(z18S) %in% z18S.protist.taxa$V1,]

## we rarefy the data
r18S.pts<-t(rrarefy(t(z18S.pts[rowSums(z18S.pts)>0,]),min(colSums(z18S.pts[rowSums(z18S.pts)>0,]))))
r18S.pts <- as.data.frame(r18S.pts)

#Finally we collapse the technical replicates (by averaging them) into a matrix 
c18S.pts <- matrix(ncol=length(unique(substr(colnames(r18S.pts),3,4))),nrow = length(r18S.pts[,1]))
colnames(c18S.pts) <- unique(substr(colnames(r18S.pts),3,4))
rownames(c18S.pts) <- rownames(r18S.pts)
for (site in unique(substr(colnames(r18S.pts),3,4))){
  c18S.pts[,site] <- round(rowMeans(r18S.pts[,substr(colnames(r18S.pts),3,4) == site]))
}
r18S.pts <- as.data.frame(c18S.pts)

write.csv(r18S.pts,"cleaned/rarefied.18S.pts.csv")




####====0.5 Control Data====####
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


####====0.6 Unused Code====####

#first we pull in the RDP data
RDP.COI.assign <- read.table("Taxonomy/RDP.assigns/RDP.class.COI.v4.txt",sep="\t")

#now lets compare the phyla level assignments to BLAST
test <- as.data.frame(table(RDP.COI.assign$V12[match(COIassignments$OTU,RDP.COI.assign$V1)],COIassignments$phylum))
test$match <- ifelse(as.character(test$Var1)==as.character(test$Var2),"match","no-match")
#get rid of blank assignments
test <- test[test$Var2!="",]
test <- test[test$Var1!="",]
#get rid of broadly unassigned
test <- test[test$Var1!="undef_undef_Eukaryota",]

#sum mismatches
correct <- sum(test$Freq[test$match=="match"])
incorrect <- sum(test$Freq[test$match=="no-match"])

correct/(correct+incorrect)*100
#~87% correct at phyla level

#now let's evaluate the 18S data
RDP.18S.assign <- read.delim("Taxonomy/RDP.assigns/RDP.class.18S.v3.2.txt",sep="\t",header = F)

#another BLAST comparison
phylaRDP <- RDP.18S.assign$V12[match(z18Sassignments$OTU,RDP.18S.assign$V1)]
#subset by cutoff
subset <- RDP.18S.assign$V11[match(z18Sassignments$OTU,RDP.18S.assign$V1)]>0.5
phylaRDP <- phylaRDP[subset]

phylaBLAST <- z18Sassignments$phylum[subset]

test <- as.data.frame(table(phylaRDP,phylaBLAST))
#get rid of blank assignments
test <- test[test$phylaRDP!="",]
test <- test[test$phylaBLAST!="",]
#get rid of broadly unassigned
test <- test[test$phylaRDP!="Eukaryota_undef_Eukaryota_undef_undef",]
test$match <- ifelse(as.character(test$phylaRDP)==as.character(test$phylaBLAST),"match","no-match")

#sum mismatches
correct <- sum(test$Freq[test$match=="match"])
incorrect <- sum(test$Freq[test$match=="no-match"])

correct/(correct+incorrect)*100
#73% - but most issues down to nomenclature differences rather than error
