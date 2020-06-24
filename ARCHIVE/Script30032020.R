##### South African eDNA Analysis #####
#Luke E. Holman 20.08.2019
####Packages####

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

#####


#### Settings and Setup####

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

#Set A plalette
palette(c('#51B9E0','#4bad84','#E8A016'))
#platette two 
#palette(c('#D36526','#2671B4'))
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

#### Next we create a dataset that only contain OTUs represented in all 3 technical replicates ####
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

#### Process the taxonomy 
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

#####OutputControlData###
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



#### Map #### suggestion - https://stackoverflow.com/questions/26114962/include-scale-and-coordinates-in-r-map

#First lets build a map with our samples 
pdf("figures/mapv1.pdf",width=7,height=5)
map("world", xlim=c(14,34), ylim=c(-36,-26), col="gray", fill=TRUE)
points(latlongdata$Long[latlongdata$type=="m"],latlongdata$Lat[latlongdata$type=="m"],pch=19,cex=1,col="darkred")
points(latlongdata$Long[latlongdata$type=="n"],latlongdata$Lat[latlongdata$type=="n"],pch=19,cex=0.5,col="darkblue")
map.scale(ratio=FALSE,cex=0.5)
map.axes()
text(latlongdata$Long[latlongdata$ori=="W" & latlongdata$type=="m"]-0.8,latlongdata$Lat[latlongdata$ori=="W" & latlongdata$type=="m"],labels=latlongdata$sitecode[ latlongdata$type=="m"& latlongdata$ori=="W"],cex=0.5)
text(latlongdata$Long[latlongdata$ori=="S"& latlongdata$type=="m"],latlongdata$Lat[latlongdata$ori=="S"& latlongdata$type=="m"]-0.5,labels=latlongdata$sitecode[latlongdata$ori=="S"& latlongdata$type=="m"],cex=0.5)
text(latlongdata$Long[latlongdata$ori=="S"& latlongdata$type=="n"],latlongdata$Lat[latlongdata$ori=="S"& latlongdata$type=="n"]-0.5,labels=latlongdata$sitecode[latlongdata$ori=="S"& latlongdata$type=="n"],cex=0.5)
text(latlongdata$Long[latlongdata$ori=="E"& latlongdata$type=="m"]+0.8,latlongdata$Lat[latlongdata$ori=="E"& latlongdata$type=="m"],labels=latlongdata$sitecode[latlongdata$ori=="E"& latlongdata$type=="m"],cex=0.5)
text(latlongdata$Long[latlongdata$ori=="E"& latlongdata$type=="n"]+0.8,latlongdata$Lat[latlongdata$ori=="E"& latlongdata$type=="n"],labels=latlongdata$sitecode[latlongdata$ori=="E"& latlongdata$type=="n"],cex=0.5)
#some need to be done by hand
#SM
text(latlongdata$Long[latlongdata$sitecode=="SM"]-0.8,latlongdata$Lat[latlongdata$sitecode=="SM"]-0.4,labels="SM",cex=0.5)
#HB
text(latlongdata$Long[latlongdata$sitecode=="HB"]-0.6,latlongdata$Lat[latlongdata$sitecode=="HB"]-0.3,labels="HB",cex=0.5)
#MN
text(latlongdata$Long[latlongdata$sitecode=="MN"]-0.4,latlongdata$Lat[latlongdata$sitecode=="MN"]-0.5,labels="MN",cex=0.5)
#NK
text(latlongdata$Long[latlongdata$sitecode=="NK"]-0.5,latlongdata$Lat[latlongdata$sitecode=="NK"]-0.3,labels="NK",cex=0.5)

dev.off()
map.scale(ratio=FALSE,cex=0.5)

#Now we build a map only with the samples in the dataset
smallLatLong <- latlongdata[latlongdata$sitecode %in% names(rCOI),]

pdf("figures/Figure1/mapv2.pdf",width=9,height=6.5)
m<- map("world", xlim=c(16,33), ylim=c(-36,-28), col="gray", fill=TRUE,mar=c(4, 4, 4, 4))

#
rect(15,-38,20,-27,border = NA,col="#78C1E1")
rect(20,-38,26.5,-27,border = NA,col="#139F73")
rect(26.5,-38,34,-27,border = NA,col="#E69F03")

m<- map("world", xlim=c(16,33), ylim=c(-36,-28), col="gray", fill=TRUE,mar=c(4, 4, 4, 4),add=TRUE)


map.scale(ratio=FALSE,cex=0.8,lwd=6,relwidth = 0.16)
addCompass(31,-34.5,cex=0.7)

xat <- pretty(m$range[1:2],n = 10)
xlab <- parse(text=degreeLabelsEW(xat))

yat <- pretty(m$range[3:4])
ylab <- parse(text=degreeLabelsNS(yat))


box()
axis(1, at=xat, labels=xlab)
axis(2, las=TRUE, at=yat, labels=ylab)
axis(3, at=xat, labels=xlab)
axis(4, las=TRUE, at=yat, labels=ylab)




points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=1,cex=2,lwd=4,col="white")
points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=1,cex=2,lwd=2,col="darkblue")


points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=4,cex=2,lwd=4,col="white")
points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=4,cex=2,lwd=2,col="darkred")

#Greyscale Os
#points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=16,cex=2.3,col="black")
#points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=16,cex=2,col="white")

#points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=1,cex=2,lwd=4,col="white")
#points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=1,cex=2,lwd=2,col="black")

dev.off()



##WestCoast
pdf("figures/map.west.pdf",width=4,height=7)
map('worldHires', xlim=c(17.5,19), ylim=c(-34.5,-32.5), col="gray", fill=TRUE,resolution=0,mar=c(4, 4, 4, 4))
map.axes()
dev.off()

##South
pdf("figures/map.south.pdf",width=7,height=4)
map('worldHires', xlim=c(21.5,28.5), ylim=c(-34.5,-32.5), col="gray", fill=TRUE,resolution=0,mar=c(4, 4, 4, 4))
map.axes()
dev.off()

##East
pdf("figures/map.east.pdf",width=7,height=4)
map('worldHires', xlim=c(30.5,32.5), ylim=c(-30.5,-28.5), col="gray", fill=TRUE,resolution=0,mar=c(4, 4, 4, 4))
map.axes()
dev.off()



###  Patterns of alpha  taxonomic diversity across the coast ####
#SppRichness
temp <- rCOI
temp[temp > 1] <- 1
Coastal.alpha <- data.frame("ID"=names(temp),"COI"= colSums(temp))
temp <- r18S
temp[temp > 1] <- 1
Coastal.alpha$z18S <- colSums(temp)
temp <- rProK
temp[temp > 1] <- 1
Coastal.alpha$ProK <- colSums(temp)
#reorder for plotting
Coastal.alpha <- Coastal.alpha[na.omit(match(latlongdata$sitecode,Coastal.alpha$ID)),]

#richness plot
pdf("figures/OTU.Richness.pdf",width = 9,height=5)
plot(1:18,Coastal.alpha$COI,xaxt='n',xlab="",ylab="OTUrichness",pch=16, col="Darkred",ylim=c(200,950))
points(1:18,Coastal.alpha$z18S,col="Darkblue",pch=16)
points(1:18,Coastal.alpha$ProK,col="Darkgreen",pch=16)
axis(1,at=1:18,labels=Coastal.alpha$ID,cex=0.7)
legend("topright",legend=c("COI","18S","16S"),col=c("Darkred","Darkblue","Darkgreen"),pch=16)
dev.off()

#richness by coast
Coastal.alpha.tall <- melt(Coastal.alpha)
Coastal.alpha.tall$coast <- as.character(latlongdata$PERMori[match(as.character(Coastal.alpha.tall$ID),as.character(latlongdata$sitecode))])


Coastal.alpha.tall$comb <- factor(paste(Coastal.alpha.tall$variable,Coastal.alpha.tall$coast,sep="."),
                                     levels=unique(paste(Coastal.alpha.tall$variable,Coastal.alpha.tall$coast,sep="."))[c(1,4,7,2,5,8,3,6,9)])


#Test for sig diff by coast 
##ANOVA assump - Are the residuals normally distributed? 

resCOI <- residuals(lm(log10(value)~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="COI",]))
shapiro.test(resCOI)
#yes if log10 transformed
res18S <- residuals(lm(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="z18S",]))
shapiro.test(res18S)
#yes
resProK <- residuals(lm(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="ProK",]))
shapiro.test(resProK)
#yes

## Do coasts have equal variance?
bartlett.test(log10(value)~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="COI",])
#yes
bartlett.test(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="z18S",])
#yes
bartlett.test(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="COI",])
#yes

ANOVA.COI <- aov(log10(value)~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="COI",])
ANOVA.18S <- aov(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="z18S",])
ANOVA.ProK <- aov(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="ProK",])
summary(ANOVA.COI) #Global P =0.051
summary(ANOVA.18S) #Global P =0.195
summary(ANOVA.ProK) #Global P =0.0065

TukeyHSD(ANOVA.COI)
TukeyHSD(ANOVA.18S)
TukeyHSD(ANOVA.ProK)

pdf("figures/OTU.coast.Richness2.pdf",width = 9,height=4)
palette(c('#54B4E9','#F0E442','#CC79A7'))
plot(jitter(as.numeric(Coastal.alpha.tall$comb),0.8),
     Coastal.alpha.tall$value,pch=16,
     xaxt="n",
     xlab="Coast",
     ylab="OTUs")
axis(1,at=1:9,labels=rep(c("COI","18S","16S"),3),cex=0.8)

##Lets draw in those mean values
for (num in 1:9){
  run.mean <- mean(Coastal.alpha.tall$value[as.numeric(Coastal.alpha.tall$comb)==num])
  lines(c(num-0.2,num+0.2),c(run.mean,run.mean),lwd=2)
}
dev.off()



#Taxonomic profile
#set NA as "" in taxonomy file

##We first make a function that calculates the count table for us in two ways
CountTable <- function(in.taxonomy,in.data,output="Count"){
  if(length(in.taxonomy)!=length(in.data[,1])){stop("Dataframe and corresponding taxonomy are not the same length")}
  in.taxonomy[is.na(in.taxonomy)] <- ""
  out.dat <- as.data.frame(matrix(ncol=length(in.data[1,]),nrow=length(unique(in.taxonomy))))
  rownames(out.dat) <- sort(unique(in.taxonomy))
  colnames(out.dat) <- colnames(in.data)    
  out.dat.abundance <- out.dat
  for (sample in 1:length(in.data[1,])){
    out.dat[,sample] <- table(in.taxonomy[in.data[,sample]>0])[match(sort(unique(in.taxonomy)),names(table(in.taxonomy[in.data[,sample]>0])))]
    out.dat.abundance[,sample] <- aggregate(in.data[,sample], by=list(Category=in.taxonomy), FUN=sum)[,2]
  }
  out.dat[is.na(out.dat)] <- 0
  rownames(out.dat)[1] <- "Unassigned"
  if(output=="Count"){return(out.dat)}else if(
    output=="Abundance"){return(out.dat.abundance)}
}
#Then we write a function for concatenating small abundance groups per sample
minAbundance <- function(inputtable=NA,minAbun= 0.01){
  inputtable <- rbind(inputtable,rep(0,dim(inputtable)[1]))
  rownames(inputtable)[dim(inputtable)[1]] <- "Others"
  for (row in 1:dim(inputtable)[2]){
    min <- sum(inputtable[,row])*minAbun
    others <- sum(inputtable[inputtable[,row]<min,row])
    inputtable[inputtable[,row]<min,row] <- 0
    inputtable["Others",row] <- others
    inputtable <- inputtable[rowSums(inputtable)>1,]
  }
  return(inputtable)
}

#Lets make those tables
COI.taxa.abun <- minAbundance(CountTable(COIassignments[match(rownames(rCOI),COIassignments$OTU),"phylum"],rCOI,output="Abundance"))
z18S.taxa.abun <- minAbundance(CountTable(z18Sassignments[match(rownames(r18S),z18Sassignments$OTU),"phylum"],r18S,output="Abundance"))
ProK.taxa.abun <- minAbundance(CountTable(as.character(ProKtaxa$Phylum),rProK,output="Abundance"))

COI.taxa.count <- minAbundance(CountTable(COIassignments[match(rownames(rCOI),COIassignments$OTU),"phylum"],rCOI,output="Count"))
z18S.taxa.count <- minAbundance(CountTable(z18Sassignments[match(rownames(r18S),z18Sassignments$OTU),"phylum"],r18S,output="Count"))
ProK.taxa.count <- minAbundance(CountTable(as.character(ProKtaxa$Phylum),rProK,output="Count"))

#Number of phyla 
CountTable(COIassignments[match(rownames(rCOI),COIassignments$OTU),"phylum"],rCOI,output="Count")

table(COIassignments[match(rownames(r18S),z18Sassignments$OTU),"phylum"])


#Add a little thing for more colours 
getPalette = colorRampPalette(brewer.pal(9, "Pastel1"))

pdf("figures/Figure1/CountTab.pdf",width = 11, height = 7)
par(mfrow=c(3,1),mai = c(0.3, 0.5, 0.1,1),xpd=TRUE)
barplot(prop.table(as.matrix(COI.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(COI.taxa.count)))]),2),col=getPalette(dim(COI.taxa.count)[1]),axisnames=FALSE)
legend(22,1,rownames(COI.taxa.count),fill=getPalette(dim(COI.taxa.count)[1]),cex=0.6)
barplot(prop.table(as.matrix(z18S.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(z18S.taxa.count)))]),2),col=getPalette(dim(z18S.taxa.count)[1]),axisnames=FALSE)
legend(22,1,rownames(z18S.taxa.count),fill =getPalette(dim(z18S.taxa.count)[1]),cex=0.6)
barplot(prop.table(as.matrix(ProK.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.count)))]),2),col=getPalette(dim(ProK.taxa.count)[1]))
legend(22,1,rownames(ProK.taxa.count),fill=getPalette(dim(ProK.taxa.count)[1]),cex=0.6)
dev.off()

pdf("figures/Figure1/CountTab.noUnassign.pdf",width = 11, height = 7)
par(mfrow=c(3,1),mai = c(0.3, 0.5, 0.2,1),xpd=TRUE)
barplot(prop.table(as.matrix(COI.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(COI.taxa.count)))]),2)[-1,],col=getPalette(dim(COI.taxa.count)[1]),axisnames=FALSE,main="COI")
legend(22,0.4,rownames(COI.taxa.count)[-1],fill=getPalette(dim(COI.taxa.count)[1]-1),cex=0.6)
barplot(prop.table(as.matrix(z18S.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(z18S.taxa.count)))]),2)[-1,],col=getPalette(dim(z18S.taxa.count)[1]),axisnames=FALSE,main="18S")
legend(22,0.75,rownames(z18S.taxa.count)[-1],fill=getPalette(dim(z18S.taxa.count)[1]-1),cex=0.6)
barplot(prop.table(as.matrix(ProK.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.count)))]),2)[-1,],col=getPalette(dim(ProK.taxa.count)[1]),main="16S")
legend(22,1,rownames(ProK.taxa.count)[-1],fill=getPalette(dim(ProK.taxa.count)[1]-1),cex=0.6)
dev.off()


pdf("figures/Figure1/AbundanceTab.pdf",width = 11, height = 7)
par(mfrow=c(3,1),mai = c(0.3, 0.5, 0.1,1),xpd=TRUE)
barplot(prop.table(as.matrix(COI.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(COI.taxa.abun)))]),2),col=getPalette(dim(COI.taxa.count)[1]),axisnames=FALSE)
legend(22,1,rownames(COI.taxa.abun),fill=getPalette(dim(COI.taxa.abun)[1]),cex=0.6)
barplot(prop.table(as.matrix(z18S.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(z18S.taxa.abun)))]),2),col=getPalette(dim(z18S.taxa.count)[1]),axisnames=FALSE)
legend(22,1,rownames(z18S.taxa.abun),fill=getPalette(dim(z18S.taxa.abun)[1]),cex=0.6)
barplot(prop.table(as.matrix(ProK.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.abun)))]),2),col=getPalette(dim(ProK.taxa.count)[1]))
legend(22,1,rownames(ProK.taxa.abun),fill=getPalette(dim(ProK.taxa.abun)[1]),cex=0.6)
dev.off()

pdf("figures/Figure1/AbundanceTab.noUnassign.pdf",width = 11, height = 7)
par(mfrow=c(3,1),mai = c(0.3, 0.6, 0.1,1.3),xpd=TRUE)
barplot(prop.table(as.matrix(COI.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(COI.taxa.abun)))]),2)[-1,],col=getPalette(dim(COI.taxa.count)[1]),axisnames=FALSE,main="COI")
legend(22,0.8,rownames(COI.taxa.abun)[-1],fill=getPalette(dim(COI.taxa.abun)[1]-1),cex=0.6)
barplot(prop.table(as.matrix(z18S.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(z18S.taxa.abun)))]),2)[-1,],col=getPalette(dim(z18S.taxa.count)[1]),axisnames=FALSE,main="18S")
legend(22,1,rownames(z18S.taxa.abun)[-1],fill=getPalette(dim(z18S.taxa.abun)[1]-1),cex=0.6)
barplot(prop.table(as.matrix(ProK.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.abun)))]),2)[-1,],col=getPalette(dim(ProK.taxa.count)[1]),main="16S")
legend(22,1,rownames(ProK.taxa.abun)[-1],fill=getPalette(dim(ProK.taxa.abun)[1]-1),cex=0.6)
dev.off()


#Lets compare to Awad2002
awad <- read.csv("metadata/awad.richness.csv")
awad.all <- read.csv("metadata/awad.total.richness.csv")
plot(jitter(as.numeric(awad.all$Coast),0.8),
     awad.all$SpeciesRichness,pch=16,col=as.numeric(awad.all$Coast),
     xaxt="n",xlab="Coast",
     ylab="OTUs")
axis(1,at=c(1,2,3),labels=c("West","South","East"),cex=0.8)
plot(jitter(as.numeric(awad$Coast),0.8),
     awad$Richness,pch=16,col=as.numeric(awad$Coast),
     xaxt="n",xlab="Coast",
     ylab="OTUs")
axis(1,at=c(1,2,3),labels=c("West","South","East"),cex=0.8)



##Lets pull in the whole awad dataset and look at the patterns. 
palette(c('#51B9E0','#2C7A59','#E8A016'))
awad.all <- read.csv("rawdata/Awad.total.csv")

#PCoA
awad.dist <- vegdist(t(awad.all), "jaccard")
pcoa <- cmdscale(awad.dist)

#nMDS
nMDS <- metaMDS(t(awad.all),"jaccard")

#Plot
pdf("figures/awad.beta.pdf",width = 9,height=5.5)
par(mfrow=c(1,2))
plot(nMDS$points,col="lightblue",pch=16,main="Awad2002 - nMDS - Jaccard",cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points),cex=0.7)

plot(pcoa,col="lightblue",pch=16,cex=3,main="Awad2002 - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa),cex=0.7)
dev.off()
#Dendrogram
par(mfrow=c(1,1))
#the below line makes some weights to reorder the leaves as much as possible to match coasts
test <- hclust(vegdist(t(awad.all), "jaccard"),method="average")
weights <- match(test$labels,latlongdata$sitecode)
pdf("figures/dend.awad.pdf",width = 9,height=5.5)
plot(reorder(hclust(vegdist(t(awad.all), "jaccard"),method="average"),as.numeric(gsub("[^0-9\\.]", "", colnames(awad.all)) )),main="Awad")
dev.off()

#For fun lets put the whole dataset toghtehr and see what it looks like!

rCombined <- rbind(r18S,rCOI,rProK)

#PCoA
combined <- vegdist(t(rCombined), "jaccard")
pcoa <- cmdscale(combined)

#nMDS
nMDS <- metaMDS(t(rCombined),"jaccard")

#Plot
pdf("figures/combined.beta.pdf",width = 9,height=5.5)
par(mfrow=c(1,2))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="Combined - nMDS - Jaccard",cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points))

plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="Combined - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
dev.off()
#Dendrogram
par(mfrow=c(1,1))
#the below line makes some weights to reorder the leaves as much as possible to match coasts
test <- hclust(vegdist(t(rCombined), "jaccard"))
weights <- match(test$labels,latlongdata$sitecode)
pdf("figures/dend.awad.pdf",width = 9,height=5.5)
plot(reorder(hclust(vegdist(t(rCombined), "jaccard"),method="average"),weights),main="Combined")
dev.off()


### Nows let's look at beta diversity 

##First COI
#PCoA
rCOI.dist <- vegdist(t(rCOI), "jaccard")
pcoa <- cmdscale(rCOI.dist)

#nMDS
nMDS <- metaMDS(t(rCOI),"jaccard")

##PERMANOVA COI
PermDispCOI <- betadisper(vegdist(t(rCOI), "jaccard"),latlongdata$PERMori[match(colnames(rCOI),latlongdata$sitecode)])
anova(PermDispCOI)
#No sig difference in multivariance homogenity
adonis(vegdist(t(rCOI), "jaccard")~latlongdata$PERMori[match(colnames(rCOI),latlongdata$sitecode)])
pairwise.adonis(vegdist(t(rCOI), "jaccard"),latlongdata$PERMori[match(colnames(rCOI),latlongdata$sitecode)])
#p<0.001


#Plot
pdf("figures/COI.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordiellipse(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            kind = "sd",
            draw = "polygon",
            lty=0,
            col=1:3)
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points))
text(-0.45,0.3,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

   #What about plotting the relationship between temperature and nMDS axes?
pdf("figures/COI.temp.nMDS.pdf",width = 5.2,height=5)
par(mar=c(4,4,1,1))
plot(nMDS$points[,1],tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)],pch=16,ylab="2-year Mean Temp (C)",xlab="nMDS Axis1",cex=1.2)
abline(lm(tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]~nMDS$points[,1]),col="red",lwd=1.3)
dev.off()

summary(lm(nMDS$points[,1]~tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]))

#Lets partition the diversity into turnover and nestedness NEW

tempCOI <-rCOI
tempCOI[tempCOI > 1] <- 1
test <- beta.pair(t(tempCOI), index.family = "jaccard")
test.amount <- beta.multi(t(tempCOI), index.family = "jaccard")


nMDSall <- metaMDS(test$beta.jac)
nMDSturn <- metaMDS(test$beta.jtu)
nMDSnest<- metaMDS(test$beta.jne)
par(mfrow=c(1,3))
plot(nMDSall$points,col=latlongdata$PERMori[match(rownames(nMDSall$points),latlongdata$sitecode)],pch=16,main="All",cex=3)
text(nMDSall$points[,1],nMDSall$points[,2],labels=rownames(nMDSall$points))
plot(nMDSturn$points,col=latlongdata$PERMori[match(rownames(nMDSturn$points),latlongdata$sitecode)],pch=16,main="Turn",cex=3)
text(nMDSturn$points[,1],nMDSturn$points[,2],labels=rownames(nMDSturn$points))
plot(nMDSnest$points,col=latlongdata$PERMori[match(rownames(nMDSnest$points),latlongdata$sitecode)],pch=16,main="Nest",cex=3)
text(nMDSnest$points[,1],nMDSnest$points[,2],labels=rownames(nMDSall$points))

#Plot a dendrogram (UPGMA)
par(mfrow=c(1,1))
#the below line makes some weights to reorder the leaves as much as possible to match coasts
weights <- match(test$labels,latlongdata$sitecode)
pdf("figures/dend.COI.pdf",width = 9,height=5.5)
plot(reorder(hclust(vegdist(t(rCOI), "jaccard"),method="average"),weights),main="COI")
dev.off()

##Now 18S
#PCoA
r18S.dist <- vegdist(t(r18S), "jaccard")
pcoa <- cmdscale(r18S.dist)

#nMDS
nMDS <- metaMDS(t(r18S),"jaccard")

##PERMANOVA 18S
PermDisp18S <- betadisper(vegdist(t(r18S), "jaccard"),latlongdata$PERMori[match(colnames(r18S),latlongdata$sitecode)])
anova(PermDisp18S)
#No sig difference in multivariance homogenity
adonis(vegdist(t(r18S), "jaccard")~latlongdata$PERMori[match(colnames(r18S),latlongdata$sitecode)])
pairwise.adonis(vegdist(t(r18S), "jaccard"),latlongdata$PERMori[match(colnames(r18S),latlongdata$sitecode)])
#p<0.001


#Plot
pdf("figures/r18S.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordiellipse(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            kind = "sd",
            draw = "polygon",
            lty=0,
            col=1:3)
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points))
text(-0.375,0.475,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

#What about plotting the relationship between temperature and nMDS axes?
pdf("figures/18S.temp.nMDS.pdf",width = 5,height=5)
par(mar=c(4,4,1,1))
plot(nMDS$points[,1],tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)],pch=16,ylab="2-year Mean Temp (C)",xlab="nMDS Axis1",cex=1.2)
abline(lm(tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]~nMDS$points[,1]),col="red",lwd=1.3)
dev.off()

summary(lm(nMDS$points[,1]~tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]))


#Lets partition the diversity into turnover and nestedness NEW

temp18S <-r18S
temp18S[temp18S > 1] <- 1
test <- beta.pair(t(temp18S), index.family = "jaccard")
test.amount <- beta.multi(t(temp18S), index.family = "jaccard")


nMDSall <- metaMDS(test$beta.jac)
nMDSturn <- metaMDS(test$beta.jtu)
nMDSnest<- metaMDS(test$beta.jne)
par(mfrow=c(1,3))
plot(nMDSall$points,col=latlongdata$PERMori[match(rownames(nMDSall$points),latlongdata$sitecode)],pch=16,main="All",cex=3)
text(nMDSall$points[,1],nMDSall$points[,2],labels=rownames(nMDSall$points))
plot(nMDSturn$points,col=latlongdata$PERMori[match(rownames(nMDSturn$points),latlongdata$sitecode)],pch=16,main="Turn",cex=3)
text(nMDSturn$points[,1],nMDSturn$points[,2],labels=rownames(nMDSturn$points))
plot(nMDSnest$points,col=latlongdata$PERMori[match(rownames(nMDSnest$points),latlongdata$sitecode)],pch=16,main="Nest",cex=3)
text(nMDSnest$points[,1],nMDSnest$points[,2],labels=rownames(nMDSall$points))

#Let's exmaine the differences in a dendrogram (UPGMA)
par(mfrow=c(1,1))
#the below line makes some weights to reorder the leaves as much as possible to match coasts
weights <- match(test$labels,latlongdata$sitecode)
pdf("figures/dend.18S.pdf",width = 9,height=5.5)
plot(reorder(hclust(vegdist(t(r18S), "jaccard"),method="average"),weights),main="18S")
dev.off()

##Now ProK
#PCoA
ProK.dist <- vegdist(t(rProK), "jaccard")
pcoa <- cmdscale(ProK.dist)

#nMDS
nMDS <- metaMDS(t(rProK),"jaccard")

##PERMANOVA ProK
PermDispProK <- betadisper(vegdist(t(rProK), "jaccard"),latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
anova(PermDispProK)
TukeyHSD(PermDispProK)
plot(TukeyHSD(PermDispProK))
#Difference in variance between East and West 
adonis(vegdist(t(rProK), "jaccard")~latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
pairwise.adonis(vegdist(t(rProK), "jaccard"),latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
#p<0.001


#Plot
pdf("figures/ProK.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordiellipse(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            kind = "sd",
            draw = "polygon",
            lty=0,
            col=1:3)
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points))
text(-0.475,0.25,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

#What about plotting the relationship between temperature and nMDS axes?
pdf("figures/16S.temp.nMDS.pdf",width = 5,height=5)
par(mar=c(4,4,1,1))
plot(nMDS$points[,1],tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)],pch=16,ylab="2-year Mean Temp (C)",xlab="nMDS Axis1",cex=1.2)
abline(lm(tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]~nMDS$points[,1]),col="red",lwd=1.3)
dev.off()

summary(lm(nMDS$points[,1]~tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]))



#Lets partition the diversity into turnover and nestedness NEW

tempProK <-rProK
tempProK[tempProK > 1] <- 1
test <- beta.pair(t(tempProK), index.family = "jaccard")
test.amount <- beta.multi(t(tempProK), index.family = "jaccard")


nMDSall <- metaMDS(test$beta.jac)
nMDSturn <- metaMDS(test$beta.jtu)
nMDSnest<- metaMDS(test$beta.jne)
par(mfrow=c(1,3))
plot(nMDSall$points,col=latlongdata$PERMori[match(rownames(nMDSall$points),latlongdata$sitecode)],pch=16,main="All",cex=3)
text(nMDSall$points[,1],nMDSall$points[,2],labels=rownames(nMDSall$points))
plot(nMDSturn$points,col=latlongdata$PERMori[match(rownames(nMDSturn$points),latlongdata$sitecode)],pch=16,main="Turn",cex=3)
text(nMDSturn$points[,1],nMDSturn$points[,2],labels=rownames(nMDSturn$points))
plot(nMDSnest$points,col=latlongdata$PERMori[match(rownames(nMDSnest$points),latlongdata$sitecode)],pch=16,main="Nest",cex=3)
text(nMDSnest$points[,1],nMDSnest$points[,2],labels=rownames(nMDSall$points))

#Let's exmaine the differences in a dendrogram (UPGMA)
#the below line makes some weights to reorder the leaves as much as possible to match coasts
weights <- match(test$labels,latlongdata$sitecode)
pdf("figures/dend.16S.pdf",width = 9,height=5.5)
plot(reorder(hclust(vegdist(t(rProK), "jaccard"),method="average"),weights),main="16S")
dev.off()

###Now lets explore the dissimilarity x distance relationship 

palette(c('#D36526','#2671B4'))

##First COI
data <- rCOI

data.dist <- vegdist(t(data), "jaccard",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA

pairwise <- melt(data.dist2)
pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
pairwise$comp <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
}

#theDistanceData
Kmdata <-read.csv("metadata/distance.csv")
rownames(Kmdata) <- Kmdata$sitecode
Kmdata <- Kmdata[,-1]
Kmdatapair <- melt(as.matrix(Kmdata))
Kmdatapair <- Kmdatapair[Kmdatapair$value != 0 & !is.na(Kmdatapair$value),]
Kmdatapair$comp <- rep(NA,length(Kmdatapair[,1]))
#silly loop to get names matching 
for (row in 1:length(Kmdatapair[,1])){
  Kmdatapair$comp[row] <- paste(sort(c(as.character(Kmdatapair$Var1[row]),as.character(Kmdatapair$Var2[row]))),collapse=".")
}
hist(Kmdatapair$value,breaks=20,col="lightblue")

##Plot time
pairwise$comp %in% Kmdatapair$comp

hist(pairwise$value,breaks=20) 
pdf("figures/Overall.COI.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="COI")
dev.off()

#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
}

#We get rid of mixed samples
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]



#Lets do some stats! 
##First lets rename everything to make it easier to understand
site.similarity <- 1-pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type


#Lets exmaine the diagnostic plots 
plot(lm(log10(site.similarity)~distance*site.type))
##They look good
model <- summary(lm(log10(site.similarity)~distance*site.type))
##Stats all signficant

sink("model.output/Dist.Decay.lm.COI.txt")
print(summary(lm(log10(site.similarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.similarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(site.similarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.COI.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/log10Type.COI.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)

lines(indata$distance[indata$site.type=="Artificial"],indata$fit[indata$site.type=="Artificial"],col="darkred",lwd=2)
polygon(x = c(indata$distance[indata$site.type=="Artificial"], rev(indata$distance[indata$site.type=="Artificial"])),
        y = c(indata$lwr[indata$site.type=="Artificial"], 
              rev(indata$upr[indata$site.type=="Artificial"])),
        col =  adjustcolor("orangered", alpha.f = 0.10), border = NA)


lines(indata$distance[indata$site.type=="Natural"],indata$fit[indata$site.type=="Natural"],col="darkblue",lwd=2)
polygon(x = c(indata$distance[indata$site.type=="Natural"], rev(indata$distance[indata$site.type=="Natural"])),
        y = c(indata$lwr[indata$site.type=="Natural"], 
              rev(indata$upr[indata$site.type=="Natural"])),
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)


dev.off()






##Now 18S
data <- r18S

data.dist <- vegdist(t(data), "jaccard",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA

pairwise <- melt(data.dist2)
pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
pairwise$comp <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
}

#theDistanceData
Kmdata <-read.csv("metadata/distance.csv")
rownames(Kmdata) <- Kmdata$sitecode
Kmdata <- Kmdata[,-1]
Kmdatapair <- melt(as.matrix(Kmdata))
Kmdatapair <- Kmdatapair[Kmdatapair$value != 0 & !is.na(Kmdatapair$value),]
Kmdatapair$comp <- rep(NA,length(Kmdatapair[,1]))
#silly loop to get names matching 
for (row in 1:length(Kmdatapair[,1])){
  Kmdatapair$comp[row] <- paste(sort(c(as.character(Kmdatapair$Var1[row]),as.character(Kmdatapair$Var2[row]))),collapse=".")
}
hist(Kmdatapair$value,breaks=20,col="lightblue")

##Plot time
pairwise$comp %in% Kmdatapair$comp

hist(pairwise$value,breaks=20) 
pdf("figures/Overall.18S.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "lightblue",main="18S")
dev.off()

#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
}

#We get rid of mixed samples
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]


#Lets do some stats! 
##First lets rename everything to make it easier to understand
site.similarity <- 1-pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type


#Lets exmaine the diagnostic plots 
plot(lm(log10(site.similarity)~distance*site.type))
##They look good
summary(lm(log10(site.similarity)~distance*site.type))
##Sig effect but no difference between site types. 

sink("model.output/Dist.Decay.lm.18S.txt")
print(summary(lm(log10(site.similarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.similarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(site.similarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.18S.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="18S")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/log10Type.18S.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="18S")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)

lines(indata$distance[indata$site.type=="Artificial"],indata$fit[indata$site.type=="Artificial"],col="darkred",lwd=2)
polygon(x = c(indata$distance[indata$site.type=="Artificial"], rev(indata$distance[indata$site.type=="Artificial"])),
        y = c(indata$lwr[indata$site.type=="Artificial"], 
              rev(indata$upr[indata$site.type=="Artificial"])),
        col =  adjustcolor("orangered", alpha.f = 0.10), border = NA)


lines(indata$distance[indata$site.type=="Natural"],indata$fit[indata$site.type=="Natural"],col="darkblue",lwd=2)
polygon(x = c(indata$distance[indata$site.type=="Natural"], rev(indata$distance[indata$site.type=="Natural"])),
        y = c(indata$lwr[indata$site.type=="Natural"], 
              rev(indata$upr[indata$site.type=="Natural"])),
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)


dev.off()


##Now ProK
data <- rProK

data.dist <- vegdist(t(data), "jaccard",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA

pairwise <- melt(data.dist2)
pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
pairwise$comp <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
}

#theDistanceData
Kmdata <-read.csv("metadata/distance.csv")
rownames(Kmdata) <- Kmdata$sitecode
Kmdata <- Kmdata[,-1]
Kmdatapair <- melt(as.matrix(Kmdata))
Kmdatapair <- Kmdatapair[Kmdatapair$value != 0 & !is.na(Kmdatapair$value),]
Kmdatapair$comp <- rep(NA,length(Kmdatapair[,1]))
#silly loop to get names matching 
for (row in 1:length(Kmdatapair[,1])){
  Kmdatapair$comp[row] <- paste(sort(c(as.character(Kmdatapair$Var1[row]),as.character(Kmdatapair$Var2[row]))),collapse=".")
}
hist(Kmdatapair$value,breaks=20,col="lightblue")

##Plot time
pairwise$comp %in% Kmdatapair$comp

hist(pairwise$value,breaks=20) 
pdf("figures/Overall.ProK.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "lightblue",main="16S")
dev.off()

#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
}

#We get rid of mixed samples
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]


#Lets do some stats! 
##First lets rename everything to make it easier to understand
site.similarity <- 1-pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type


#Lets exmaine the diagnostic plots 
plot(lm(log10(site.similarity)~distance*site.type))
##They look good
summary(lm(log10(site.similarity)~distance*site.type))
##Sig effect but no difference between site types. 

sink("model.output/Dist.Decay.lm.16S.txt")
print(summary(lm(log10(site.similarity)~distance*site.type)))
sink()


#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.similarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(site.similarity)~distance*site.type),indata,interval="confidence"))



#NowPlot to visualise effect
pdf("figures/Type.ProK.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="16S")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/log10Type.ProK.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="16S")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)

lines(indata$distance[indata$site.type=="Artificial"],indata$fit[indata$site.type=="Artificial"],col="darkred",lwd=2)
polygon(x = c(indata$distance[indata$site.type=="Artificial"], rev(indata$distance[indata$site.type=="Artificial"])),
        y = c(indata$lwr[indata$site.type=="Artificial"], 
              rev(indata$upr[indata$site.type=="Artificial"])),
        col =  adjustcolor("orangered", alpha.f = 0.10), border = NA)


lines(indata$distance[indata$site.type=="Natural"],indata$fit[indata$site.type=="Natural"],col="darkblue",lwd=2)
polygon(x = c(indata$distance[indata$site.type=="Natural"], rev(indata$distance[indata$site.type=="Natural"])),
        y = c(indata$lwr[indata$site.type=="Natural"], 
              rev(indata$upr[indata$site.type=="Natural"])),
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)

dev.off()



#### Now lets seperate out the phyla and check for the same patterns $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


dataset <- rCOI[]
loopdatasetname <- "Bacillariophyta"

palette(c('#3473B1','#2C7A59','#E8A016'))

megaplot <- function(dataset,loopdatasetname,location){
palette(c('#3473B1','#2C7A59','#E8A016'))

title <- paste0(loopdatasetname," n=",dim(dataset)[1])
pdf(paste0(location,loopdatasetname,".plot.pdf"),width = 8,height=11)
par(mfrow=c(3,2))
temp <- dataset[colSums(dataset)>0]
dataset <- dataset[colSums(dataset)>0]

#Alpha richness per coast (sig for relationship) 
Loop.Coastal.alpha <- data.frame("ID"=names(temp),"COI"= colSums(temp))
temp[temp > 1] <- 1
Loop.Coastal.alpha$Data <- colSums(temp)
#reorder for plotting
Loop.Coastal.alpha <- Loop.Coastal.alpha[na.omit(match(latlongdata$sitecode,Loop.Coastal.alpha$ID)),]
Loop.Coastal.alpha$coast <- as.character(latlongdata$PERMori[match(as.character(Loop.Coastal.alpha$ID),as.character(latlongdata$sitecode))])


Loop.ANOVA <- aov(Data~coast,data=Loop.Coastal.alpha)

Loop.ANOVA <- lm(Data~coast,data=Loop.Coastal.alpha)
temp <- summary(Loop.ANOVA)
temp$coefficients[,4]


plot(jitter(as.numeric(as.factor(Loop.Coastal.alpha$coast)),0.6),
     Loop.Coastal.alpha$Data,pch=16,col=as.numeric(as.factor(Loop.Coastal.alpha$coast)),
     xaxt="n",xlab="Coast",
     ylab="OTUs",
     main=title)
axis(1,at=c(1,2,3),labels=c("West","South","East"),cex=0.8)

##Lets draw in those mean values
for (num in 1:3){
  run.mean <- mean(Loop.Coastal.alpha$Data[as.numeric(as.factor(Loop.Coastal.alpha$coast))==num])
  lines(c(num-0.1,num+0.1),c(run.mean,run.mean),lwd=2)
}
legend("topright",legend=c("COI","18S","16S"),col=1:3,pch=16)

#Dendrogram and nMDS for jaccard (sig for PERMANOVA) 
test <- hclust(vegdist(t(dataset), "jaccard"))
weights <- match(test$labels,latlongdata$sitecode)
plot(reorder(hclust(vegdist(t(dataset), "jaccard"),method="average"),weights),main=title)

#nMDS
nMDS <- metaMDS(t(dataset),"jaccard")
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main=paste0("nMDS ",title," stress=",round(nMDS$stress,3)),cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points))

#1st nMDS axis to temp 
temp <- summary(lm(nMDS$points[,1]~tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]))
round(temp$coefficients[2,4],5)
plot(tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)],nMDS$points[,1],pch=16,col="lightblue",xlab="2-year Mean Temp (C)",ylab="nMDS Axis1",main=paste0(title," P=",round(temp$coefficients[2,4],5)))
text(tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)],nMDS$points[,1],labels=names(nMDS$points[,1]))
abline(lm(nMDS$points[,1]~tempdat$TempAvr[match(names(nMDS$points[,1]),tempdat$Site.Code)]),col="red")

palette(c('#D36526','#2671B4'))

#Distance decay relationship (stats)
loop.data.dist <- vegdist(t(dataset), "jaccard",upper = FALSE)
loop.data.dist2 <- as.matrix(loop.data.dist)
loop.data.dist2[upper.tri(loop.data.dist2)] <- NA

pairwise <- melt(loop.data.dist2)
pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
pairwise$comp <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
}
#Plot the overall realtionship
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "lightblue",main=title)


#get rid of the mixed comparisons
pairwise$type <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
}

#We get rid of mixed samples
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]

#stats
site.similarity <- 1-pairwise.nomixed$value
index.withoutzeros <- which(!site.similarity==0)

if(!length(index.withoutzeros)==0){

site.similarity <- site.similarity[index.withoutzeros]
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
distance <- distance[index.withoutzeros]
site.type <- pairwise.nomixed$type
site.type <- site.type[index.withoutzeros]

if(length(unique(site.type))>1 & length(site.similarity)>4){

loopdecay <- summary(lm(log10(site.similarity)~distance*site.type))

paste0("Padd=",round(loopdecay$coefficients,5)[3,4]," Pint=",round(loopdecay$coefficients,5)[4,4])

plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main=paste0("Padd=",round(loopdecay$coefficients,5)[3,4]," Pint=",round(loopdecay$coefficients,5)[4,4]))
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
x <-Kmdatapair$value[match(pairwise.nomixed$comp[index.withoutzeros],Kmdatapair$comp)]
y <- log10(1-pairwise.nomixed$value[index.withoutzeros])
abline(lm(y[site.type=="Artificial"] ~ x[site.type=="Artificial"]),lwd=2,col=1)
abline(lm(y[site.type=="Natural"] ~ x[site.type=="Natural"]),lwd=2,col=2)
}}
dev.off()
}

#loop over phyla in each marker
#COI 
for (phylum in as.character(unique(COIassignments$phylum))){
  
  temp <- COIassignments[!is.na(COIassignments$phylum),]
  runningdat <- rCOI[na.omit(match(as.character(temp$OTU[which(temp$phylum==phylum)]),rownames(rCOI))),]
  if (dim(runningdat)[1]==0){next()}
  if (dim(runningdat)[1]<20){next()}
  megaplot(runningdat,phylum,"figures/phyla.data/COI/")
}

for (phylum in as.character(unique(z18Sassignments$phylum))){
  
  temp <- z18Sassignments[!is.na(z18Sassignments$phylum),]
  runningdat <- r18S[na.omit(match(as.character(temp$OTU[which(temp$phylum==phylum)]),rownames(r18S))),]
  if (dim(runningdat)[1]==0){next()}
  if (dim(runningdat)[1]<20){next()}
  megaplot(runningdat,phylum,"figures/phyla.data/18S/")
}

for (phylum in as.character(unique(ProKtaxa$Phylum))){
  
  temp <- ProKtaxa[!is.na(ProKtaxa$Phylum),]
  runningdat <- rProK[na.omit(match(as.character(temp$OTUname[which(temp$Phylum==phylum)]),rownames(rProK))),]
  if (dim(runningdat)[1]==0){next()}
  if (dim(runningdat)[1]<20){next()}
  megaplot(runningdat,phylum,"figures/phyla.data/16S/")
}

runningdat <- rCOI[na.omit(match(as.character(COIassignments$OTU[which(!COIassignments$assignmentQual=="None")]),rownames(rCOI))),]
megaplot(runningdat,"All","figures/phyla.data/COI/")
runningdat <- r18S[na.omit(match(as.character(z18Sassignments$OTU[which(!z18Sassignments$assignmentQual=="None")]),rownames(r18S))),]
megaplot(runningdat,"All","figures/phyla.data/18S/")
runningdat <- rProK[na.omit(match(as.character(ProKtaxa$OTUname[which(!COIassignments$assignmentQual=="None")]),rownames(rCOI))),]
megaplot(runningdat,"All","figures/phyla.data/COI/")

##pull in invasive data for COI
NNScoitax <- read.csv("Taxonomy/NNS/COImatched.csv")
as.character(NNScoitax$OTU[NNScoitax$Status=="nonnative"])
runningdat <- rCOI[na.omit(match(as.character(as.character(NNScoitax$OTU[NNScoitax$Status=="nonnative"])),rownames(rCOI))),]
megaplot(runningdat,"NNS","figures/phyla.data/")


####Now lets get everyhting into a table for analysis

HeaderOut <- c("ID","Nobs","Pperm.coast","Pdisp","DistDecayPsiteadd","DistDecayPsiteInt")
##COI

COIphylumstats <- matrix(ncol=length(HeaderOut),nrow=length(as.character(unique(COIassignments$phylum))))
colnames(COIphylumstats) <- HeaderOut
rownames(COIphylumstats) <- as.character(unique(COIassignments$phylum))


for (phylum in as.character(unique(COIassignments$phylum))){
  runningstats <- rep(NA,length(HeaderOut))
  temp <- COIassignments[!is.na(COIassignments$phylum),]
  runningdat <- rCOI[na.omit(match(as.character(temp$OTU[which(temp$phylum==phylum)]),rownames(rCOI))),]
  if (dim(runningdat)[1]==0){next()}
  if (dim(runningdat)[1]<20){next()}
  test <- runningdat
  test[test>0] <- 1
  rowSums(test)
  runningstats[1] <- phylum
  runningstats[2] <- dim(runningdat)[1]
  
  #PERMANOVA
  ###
  if(!any(vegdist(t(runningdat), "jaccard")=="NaN")){
  perm <- adonis(vegdist(t(runningdat), "jaccard")~latlongdata$PERMori[match(colnames(runningdat),latlongdata$sitecode)])
  runningstats[3] <- perm$aov.tab$`Pr(>F)`[1]
  
  dispersion <- betadisper(vegdist(t(runningdat), "jaccard"),latlongdata$PERMori[match(colnames(runningdat),latlongdata$sitecode)])
  sig <- anova(dispersion)
  runningstats[4] <- round(sig$`Pr(>F)`[1],4)
  
  #DistDecay
  loop.data.dist <- vegdist(t(runningdat), "jaccard",upper = FALSE)
  loop.data.dist2 <- as.matrix(loop.data.dist)
  loop.data.dist2[upper.tri(loop.data.dist2)] <- NA
  
  pairwise <- melt(loop.data.dist2)
  pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
  pairwise$comp <- rep(NA,length(pairwise[,1]))
  for (row in 1:length(pairwise[,1])){
    pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
  }
  
  pairwise$type <- rep(NA,length(pairwise[,1]))
  for (row in 1:length(pairwise[,1])){
    first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
    second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
    if(length(unique(c(first,second)))>1){
      pairwise$type[row] <- "Mixed"
    } else if (unique(c(first,second))=="m"){
      pairwise$type[row] <- "Artificial"
    } else if (unique(c(first,second))=="n"){
      pairwise$type[row] <- "Natural"}
  }
  
  #We get rid of mixed samples
  pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]
  
  #stats
  site.similarity <- 1-pairwise.nomixed$value
  index.withoutzeros <- which(!site.similarity==0)
  
  if(!length(index.withoutzeros)==0){
    
    site.similarity <- site.similarity[index.withoutzeros]
    distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
    distance <- distance[index.withoutzeros]
    site.type <- pairwise.nomixed$type
    site.type <- site.type[index.withoutzeros]
    
    if(length(unique(site.type))>1 & length(site.similarity)>4){
      
      loopdecay <- summary(lm(log10(site.similarity)~distance*site.type))
      
      runningstats[5] <- round(loopdecay$coefficients,5)[3,4]
      runningstats[6] <- round(loopdecay$coefficients,5)[4,4]
      
    }}
  }
  COIphylumstats[match(runningstats[1],rownames(COIphylumstats)),] <- runningstats
}
  

##18S

z18Sphylumstats <- matrix(ncol=length(HeaderOut),nrow=length(as.character(unique(z18Sassignments$phylum))))
colnames(z18Sphylumstats) <- HeaderOut
rownames(z18Sphylumstats) <- as.character(unique(z18Sassignments$phylum))


for (phylum in as.character(unique(z18Sassignments$phylum))){
  runningstats <- rep(NA,length(HeaderOut))
  temp <- z18Sassignments[!is.na(z18Sassignments$phylum),]
  runningdat <- r18S[na.omit(match(as.character(temp$OTU[which(temp$phylum==phylum)]),rownames(r18S))),]
  if (dim(runningdat)[1]==0){next()}
  if (dim(runningdat)[1]<20){next()}
  test <- runningdat
  test[test>0] <- 1
  rowSums(test)
  runningstats[1] <- phylum
  runningstats[2] <- dim(runningdat)[1]
  
  #PERMANOVA
  ###
  if(!any(vegdist(t(runningdat), "jaccard")=="NaN")){
    perm <- adonis(vegdist(t(runningdat), "jaccard")~latlongdata$PERMori[match(colnames(runningdat),latlongdata$sitecode)])
    runningstats[3] <- perm$aov.tab$`Pr(>F)`[1]
    
    dispersion <- betadisper(vegdist(t(runningdat), "jaccard"),latlongdata$PERMori[match(colnames(runningdat),latlongdata$sitecode)])
    sig <- anova(dispersion)
    runningstats[4] <- round(sig$`Pr(>F)`[1],4)
    
    #DistDecay
    loop.data.dist <- vegdist(t(runningdat), "jaccard",upper = FALSE)
    loop.data.dist2 <- as.matrix(loop.data.dist)
    loop.data.dist2[upper.tri(loop.data.dist2)] <- NA
    
    pairwise <- melt(loop.data.dist2)
    pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
    pairwise$comp <- rep(NA,length(pairwise[,1]))
    for (row in 1:length(pairwise[,1])){
      pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
    }
    
    pairwise$type <- rep(NA,length(pairwise[,1]))
    for (row in 1:length(pairwise[,1])){
      first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
      second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
      if(length(unique(c(first,second)))>1){
        pairwise$type[row] <- "Mixed"
      } else if (unique(c(first,second))=="m"){
        pairwise$type[row] <- "Artificial"
      } else if (unique(c(first,second))=="n"){
        pairwise$type[row] <- "Natural"}
    }
    
    #We get rid of mixed samples
    pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]
    
    #stats
    site.similarity <- 1-pairwise.nomixed$value
    index.withoutzeros <- which(!site.similarity==0)
    
    if(!length(index.withoutzeros)==0){
      
      site.similarity <- site.similarity[index.withoutzeros]
      distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
      distance <- distance[index.withoutzeros]
      site.type <- pairwise.nomixed$type
      site.type <- site.type[index.withoutzeros]
      
      if(length(unique(site.type))>1 & length(site.similarity)>4){
        
        loopdecay <- summary(lm(log10(site.similarity)~distance*site.type))
        
        runningstats[5] <- round(loopdecay$coefficients,5)[3,4]
        runningstats[6] <- round(loopdecay$coefficients,5)[4,4]
        
      }}
  }
  z18Sphylumstats[match(runningstats[1],rownames(z18Sphylumstats)),] <- runningstats
}

##ProK

ProKphylumstats <- matrix(ncol=length(HeaderOut),nrow=length(as.character(unique(ProKtaxa$Phylum))))
colnames(ProKphylumstats) <- HeaderOut
rownames(ProKphylumstats) <- as.character(unique(ProKtaxa$Phylum))


for (phylum in as.character(unique(ProKtaxa$Phylum))){
  runningstats <- rep(NA,length(HeaderOut))
  temp <- ProKtaxa[!is.na(ProKtaxa$Phylum),]
  runningdat <- rProK[na.omit(match(as.character(temp$OTU[which(temp$Phylum==phylum)]),rownames(r18S))),]
  if (dim(runningdat)[1]==0){next()}
  if (dim(runningdat)[1]<20){next()}
  test <- runningdat
  test[test>0] <- 1
  rowSums(test)
  runningstats[1] <- phylum
  runningstats[2] <- dim(runningdat)[1]
  
  #PERMANOVA
  ###
  if(!any(vegdist(t(runningdat), "jaccard")=="NaN")){
    perm <- adonis(vegdist(t(runningdat), "jaccard")~latlongdata$PERMori[match(colnames(runningdat),latlongdata$sitecode)])
    runningstats[3] <- perm$aov.tab$`Pr(>F)`[1]
    
    dispersion <- betadisper(vegdist(t(runningdat), "jaccard"),latlongdata$PERMori[match(colnames(runningdat),latlongdata$sitecode)])
    sig <- anova(dispersion)
    runningstats[4] <- round(sig$`Pr(>F)`[1],4)
    
    #DistDecay
    loop.data.dist <- vegdist(t(runningdat), "jaccard",upper = FALSE)
    loop.data.dist2 <- as.matrix(loop.data.dist)
    loop.data.dist2[upper.tri(loop.data.dist2)] <- NA
    
    pairwise <- melt(loop.data.dist2)
    pairwise <- pairwise[pairwise$value != 0 & !is.na(pairwise$value),]
    pairwise$comp <- rep(NA,length(pairwise[,1]))
    for (row in 1:length(pairwise[,1])){
      pairwise$comp[row] <- paste(sort(c(as.character(pairwise$Var1[row]),as.character(pairwise$Var2[row]))),collapse=".")
    }
    
    pairwise$type <- rep(NA,length(pairwise[,1]))
    for (row in 1:length(pairwise[,1])){
      first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
      second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
      if(length(unique(c(first,second)))>1){
        pairwise$type[row] <- "Mixed"
      } else if (unique(c(first,second))=="m"){
        pairwise$type[row] <- "Artificial"
      } else if (unique(c(first,second))=="n"){
        pairwise$type[row] <- "Natural"}
    }
    
    #We get rid of mixed samples
    pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]
    
    #stats
    site.similarity <- 1-pairwise.nomixed$value
    index.withoutzeros <- which(!site.similarity==0)
    
    if(!length(index.withoutzeros)==0){
      
      site.similarity <- site.similarity[index.withoutzeros]
      distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
      distance <- distance[index.withoutzeros]
      site.type <- pairwise.nomixed$type
      site.type <- site.type[index.withoutzeros]
      
      if(length(unique(site.type))>1 & length(site.similarity)>4){
        
        loopdecay <- summary(lm(log10(site.similarity)~distance*site.type))
        
        runningstats[5] <- round(loopdecay$coefficients,5)[3,4]
        runningstats[6] <- round(loopdecay$coefficients,5)[4,4]
        
      }}
  }
  ProKphylumstats[match(runningstats[1],rownames(ProKphylumstats)),] <- runningstats
}

##Clean up Tables

COIphylumstats <- COIphylumstats[!is.na(COIphylumstats[,3]),]
z18Sphylumstats <- z18Sphylumstats[!is.na(z18Sphylumstats[,3]),]
ProKphylumstats <- ProKphylumstats[!is.na(ProKphylumstats[,3]),]


write.csv(COIphylumstats,file="model.output/COIphylumstats.csv")
write.csv(z18Sphylumstats,file="model.output/z18Sphylumstats.csv")
write.csv(ProKphylumstats,file="model.output/ProKphylumstats.csv")




