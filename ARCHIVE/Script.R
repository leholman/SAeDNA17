###### Analysis of raw data for Sout Africa
#Started 24.08.2018

#Load in packages
library("dplyr")
library("vegan")
library("reshape")
library("reshape2")
library("ggplot2")
library("stringr")
library("tidyr")
library("RColorBrewer")
library("metacoder")
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
library("RgoogleMaps")
library("reshape")

#Set some variables 
minreads <- 3
items <- NULL

#Set the seed 
set.seed("123456")

#Set wd and load up sediment data 
setwd("~/Desktop/Analysis.SA/rawdata/")
metadat <- read.csv(file="../metadata/locations.csv")
latlongdata <- read.csv("../metadata/mapdata.csv")

#### Data Cleaning ####
files <- system2('ls',stdout=TRUE)

for (file in  files){
  
  rawdat <-read.csv(file=file)
  
  
  #Seperate controls and samples
  samples <- rawdat[substr(colnames(rawdat),3,5) %in% substr(metadat$RealID[metadat$Type=="sample"],3,5)]
  controls <- rawdat[substr(colnames(rawdat),3,5) %in% substr(metadat$RealID[metadat$Type=="control"],3,5)]
  
  
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
  newname <- paste("../cleaned/","Cleaned.",file,sep="")
  write.csv(samples,file=newname)
  
}

#### All samples cleaned ####

## First lest examine overall beta diversity pattersn with some NMDS diagrams

##COI

cleanedCOI <- read.csv("../cleaned/Cleaned.Leraylulu.unoise3.csv",row.names = 1)
rCOI<-
  cleanedCOI%>%
  filter(rowSums(.)>0)%>%
  t(.)%>%
  rrarefy(.,min(rowSums(.)))%>%
  t(.)

MDSCOI <- metaMDS(t(rCOI),distance = "bray")
MDSCOIdat <-as.data.frame(MDSCOI$points)
plot(MDSCOIdat,type="n")
text(MDSCOIdat$MDS1,MDSCOIdat$MDS2,labels=substr(colnames(cleanedCOI),3,5))


##18S

cleaned18S <- read.csv("../cleaned/Cleaned.Zhanlulu.unoise3.csv",row.names = 1)
r18S<-
  cleaned18S%>%
  filter(rowSums(.)>0)%>%
  t(.)%>%
  rrarefy(.,min(rowSums(.)))%>%
  t(.)

MDS18S <- metaMDS(t(r18S),distance = "bray")
MDS18Sdat <-as.data.frame(MDS18S$points)
plot(MDS18Sdat,type="n")
text(MDS18Sdat$MDS1,MDS18Sdat$MDS2,labels=substr(colnames(cleaned18S),3,5))



##OK these look nice but clearly the durban lagoons are not marine! Lets take them out ( and also the natural sites) and replot

cleaned18S<-cbind(cleaned18S[,1:6],cleaned18S[,11:62])

#also get rid of all non marina sites
cleaned18S<-cleaned18S[,substr(colnames(cleaned18S),3,4) %in% unique(substr(metadat$RealID[metadat$sitetype=="m"],3,4))]
#above code keeps the PE.B site in, lets get rid of it
cleaned18S<-cleaned18S[-grep("PE.B",colnames(cleaned18S))]

#now for the COI

cleanedCOI<-cbind(cleanedCOI[,1:6],cleanedCOI[,11:62])

#also get rid of all non marina sites
cleanedCOI<-cleanedCOI[,substr(colnames(cleanedCOI),3,4) %in% unique(substr(metadat$RealID[metadat$sitetype=="m"],3,4))]
#above code keeps the PE.B site in, lets get rid of it
cleanedCOI<-cleanedCOI[-grep("PE.B",colnames(cleanedCOI))]

#Replot
r18S<-
  cleaned18S%>%
  filter(rowSums(.)>0)%>%
  t(.)%>%
  rrarefy(.,min(rowSums(.)))%>%
  t(.)

MDS18S <- metaMDS(t(r18S),distance = "bray")
MDS18Sdat <-as.data.frame(MDS18S$points)
palette(c('#a6cee3','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#8dd3c7','#8dd3c7','#8dd3c7','#b15928'))

pdf(file = "../figures/18SnMDS.pdf",width=7.5,height=6)
plot(-MDS18Sdat$MDS1,MDS18Sdat$MDS2,col=as.factor(substr(colnames(cleaned18S),3,4)),pch=16,cex=2,xlab="MDS1",ylab="MDS2")
textdata <- aggregate(MDS18Sdat, by=list(as.factor(substr(colnames(cleaned18S),3,4))),FUN=mean)
text(-textdata$MDS1-0.07,textdata$MDS2,labels=textdata$Group.1,cex=1.5)
dev.off()

rCOI<-
  cleanedCOI%>%
  filter(rowSums(.)>0)%>%
  t(.)%>%
  rrarefy(.,min(rowSums(.)))%>%
  t(.)

MDSCOI <- metaMDS(t(rCOI),distance = "bray")
MDSCOIdat <-as.data.frame(MDSCOI$points)

pdf(file = "../figures/COInMDS.pdf",width=7.5,height=6)
plot(-MDSCOIdat$MDS1,MDSCOIdat$MDS2,col=as.factor(substr(colnames(cleanedCOI),3,4)),pch=16,cex=2,xlab="MDS1",ylab="MDS2")
textdata <- aggregate(MDSCOIdat, by=list(as.factor(substr(colnames(cleanedCOI),3,4))),FUN=mean)
text(-textdata$MDS1-0.2,textdata$MDS2,labels=textdata$Group.1,cex=1.5)
dev.off()

##These look great, now lets do a permanova to make sure Eco-region is valid 

#lets start by setting up the data
permCOI <- t(cleanedCOI)
sites <- sort(rep(1:13,3))
ecoregions <- latlongdata$PERMori[match(substr(colnames(cleanedCOI),3,4),latlongdata$sitecode)]
ecoPERMCOI <- adonis(permCOI~ecoregions+sites,method="bray",perm=999)

perm18S <- t(cleaned18S)
ecoPERM18S <-  adonis(perm18S~ecoregions+sites,method="bray",perm=999)

##Here are some tests for mulitvariate homogeneity of variance 
permutest(betadisper(vegdist(permCOI,method="bray"),group=ecoregions),pairwise = TRUE)
permutest(betadisper(vegdist(perm18S,method="bray"),group=ecoregions),pairwise = TRUE)

##Taxonomy 
#the blastexpression used is below
#blastn -query Zhan.unoise3.OTUs.fa -db SILVA.fasta -outfmt "6 qseqid stitle pident length evalue bitscore" -num_threads 4 -max_target_seqs 20 -perc_identity 99 -culling_limit 10 -evalue 1e-25 -out resultsone.txt
r18Srawtaxonomy <-read.csv(file="../OTUs/SILVA/resultsone.txt",sep="\t",header=F)
r18Staxonomy <- data.frame("OTU"=unique(r18Srawtaxonomy$V1),"assignment"=as.character(rep(NA,length(unique(r18Srawtaxonomy$V1)))),stringsAsFactors=FALSE)


loopOTU <- "Zotu15"

for (loopOTU in unique(r18Srawtaxonomy$V1)){
  
  loopassignments <- r18Srawtaxonomy[r18Srawtaxonomy$V1==loopOTU & r18Srawtaxonomy$V6==max(r18Srawtaxonomy$V6[r18Srawtaxonomy$V1==loopOTU]),]
  if(length(grep("ncultured|metagenome|environmental|eukaryote|[Mm]arine",loopassignments$species))>0){loopassignments <- loopassignments[-grep("ncultured|metagenome|environmental|eukaryote|[Mm]arine",loopassignments$species),]}
  if(length(unique(loopassignments$species))==1){r18Staxonomy$assignment[r18Staxonomy$OTU==loopOTU] <- unique(loopassignments$species)}
  
}




##Nice now lets create an ascidian dataset
COIrawtaxonomy <- read.csv(file="../OTUs/MIDORI/Leray.UNOISE3.RDP.MIDORI.csv")
COIasc <- cleanedCOI[na.omit(match(as.character(COIrawtaxonomy$out[COIrawtaxonomy$C.name=="Ascidiacea" & COIrawtaxonomy$C.c>0.9]),rownames(cleanedCOI))),]
COIasc$taxonomy <- COIrawtaxonomy$S.name[COIrawtaxonomy$out %in% rownames(COIasc)]

r18Srawtaxonomy <-read.csv(file="../OTUs/SILVA/resultsone.txt",sep="\t",header=F)
#This expression uses regex to take everything after the last semicolon. It is broken by non alpha numeric characters
r18Srawtaxonomy$species  <- gsub(".*;([A-z0-9 -]*)$","\\1",r18Srawtaxonomy$V2)
#Some of the entries ae amnbiguous eg. Flabellual sp. so the extra characters breaks the regex. We can seperate these out now and 
#make a clean dataset
r18Srawtaxonomy <- r18Srawtaxonomy[nchar(r18Srawtaxonomy$species)<50,]
ascidians18Staxonomy <- r18Srawtaxonomy[grep("Ascidiacea",r18Srawtaxonomy$V10),]
asc18S <- cleaned18S[na.omit(match(as.character(ascidians18Staxonomy$V1),rownames(cleaned18S))),]
asc18S$taxonomy <- ascidians18Staxonomy$V12[as.character(ascidians18Staxonomy$V1) %in% rownames(asc18S)]

#now we see how the range interacts with time
histdat <- read.csv("../metadata/historical.csv")
histdat$Year <- as.factor(histdat$Year)
#RABS
pdf(file="../figures/RABSrangeexp.pdf",width = 7,height = 5)
plot(as.numeric(histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Year),histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Range,type="n",xaxt="n",ylab="Distance between most distant sites (kms)",xlab="Survey Date",pch=0,ylim=c(0,2200))
points(histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Range,type="o",pch=0)
points(histdat[histdat$Species=="C. lepadiformis" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="C. lepadiformis" & histdat$Survey !="E",]$Range,type="o",pch=6)
points(histdat[histdat$Species=="M. squamiger" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="M. squamiger" & histdat$Survey !="E",]$Range,type="o",pch=5)
points(histdat[histdat$Species=="S. plicata" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="S. plicata" & histdat$Survey !="E",]$Range,type="o",pch=1)
axis(1,at=1:7,labels=levels(histdat$Year))
legend("topleft",levels(histdat$Species)[-2],pch=c(0,6,5,1),lty=1,text.font=3)
dev.off()

#RABS + eDNA
pdf(file="../figures/RABSrangeexpE.pdf",width = 7,height = 5)
plot(as.numeric(histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Year),histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Range,type="n",xaxt="n",ylab="Distance between most distant sites (kms)",xlab="Survey Date",pch=0,ylim=c(0,2200))
points(histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="A. aspersa" & histdat$Survey !="E",]$Range,type="o",pch=0)
points(histdat[histdat$Species=="C. lepadiformis" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="C. lepadiformis" & histdat$Survey !="E",]$Range,type="o",pch=6)
points(histdat[histdat$Species=="M. squamiger" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="M. squamiger" & histdat$Survey !="E",]$Range,type="o",pch=5)
points(histdat[histdat$Species=="S. plicata" & histdat$Survey !="E",]$Year,histdat[histdat$Species=="S. plicata" & histdat$Survey !="E",]$Range,type="o",pch=1)
axis(1,at=1:7,labels=levels(histdat$Year))
legend("topleft",levels(histdat$Species)[-2],pch=c(0,6,5,1),lty=1,text.font=3)
lines(c(6,7),c(943,1080),col="red",type="b",pch=0)
lines(c(6,7),c(1270,1844),col="red",type="b",pch=5)
dev.off()


#eDNA
pdf(file="../figures/eDNArangeexp.pdf",width = 7,height = 5)
plot(as.numeric(histdat[histdat$Species=="A. aspersa" & histdat$Survey !="R",]$Year),histdat[histdat$Species=="A. aspersa" & histdat$Survey !="R",]$Range,type="b",xaxt="n",ylab="Distance between most distant sites (kms)",xlab="Survey Date",pch=0,ylim=c(0,2200))
points(histdat[histdat$Species=="C. intestinalis" & histdat$Survey !="R",]$Year,histdat[histdat$Species=="C. intestinalis" & histdat$Survey !="R",]$Range,type="o",pch=2)
points(histdat[histdat$Species=="C. lepadiformis" & histdat$Survey !="R",]$Year,histdat[histdat$Species=="C. lepadiformis" & histdat$Survey !="R",]$Range,type="o",pch=6)
points(histdat[histdat$Species=="M. squamiger" & histdat$Survey !="R",]$Year,histdat[histdat$Species=="M. squamiger" & histdat$Survey !="R",]$Range,type="o",pch=5)
points(histdat[histdat$Species=="S. plicata" & histdat$Survey !="R",]$Year,histdat[histdat$Species=="S. plicata" & histdat$Survey !="R",]$Range,type="o",pch=1)
axis(1,at=1:7,labels=levels(histdat$Year))
legend("topleft",levels(histdat$Species),pch=c(0,2,6,5,1),lty=1)
dev.off()

?line

##Now lets make a map
##Read in some sites 
latlongdata <- read.csv("../metadata/mapdata.csv")

#Plot the map
map("world", xlim=c(14,34), ylim=c(-36,-26), col="gray", fill=TRUE)
points(latlongdata$Long[latlongdata$type=="m"],latlongdata$Lat[latlongdata$type=="m"],pch=19,cex=1,col="darkred")
map.scale(ratio=FALSE,cex=0.5)
text(latlongdata$Long[latlongdata$ori=="W" & latlongdata$type=="m"]-0.8,latlongdata$Lat[latlongdata$ori=="W" & latlongdata$type=="m"],labels=latlongdata$sitecode[ latlongdata$type=="m"& latlongdata$ori=="W"],cex=1.5)
text(latlongdata$Long[latlongdata$ori=="S"& latlongdata$type=="m"],latlongdata$Lat[latlongdata$ori=="S"& latlongdata$type=="m"]-0.5,labels=latlongdata$sitecode[latlongdata$ori=="S"& latlongdata$type=="m"],cex=1.5)
text(latlongdata$Long[latlongdata$ori=="E"& latlongdata$type=="m"]+0.8,latlongdata$Lat[latlongdata$ori=="E"& latlongdata$type=="m"],labels=latlongdata$sitecode[latlongdata$ori=="E"& latlongdata$type=="m"],cex=1.5)

#Now lets plot some survey graphics

##Pull in the data
RAdata <- read.csv(file="../RAresults/Rawdata.RABS.csv",header=TRUE)
RAdata.2 <- read.csv(file="../RAresults/Rawdata.eDNA.PA.csv",header=TRUE)
RAdata.3 <- read.csv(file="../RAresults/Rawdata.eDNA.csv",header=TRUE)
RAadjusted <- melt(RAdata, id=c("Species"))
RAadjusted.2 <- melt(RAdata.2, id=c("Species"))
RAadjusted.3 <- melt(RAdata.3, id=c("Species"))
str(RAadjusted)

#RABS
pdf(file="../figures/RABSsurvey.pdf",width=9,height = 5.5)
par(mar=c(2.1,11.1,4.1,2.1))
plot(as.numeric(RAadjusted$variable),as.numeric(RAadjusted$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,9.5))
points(as.numeric(RAadjusted$variable),as.numeric(factor(RAadjusted$Species, levels=rev(levels(RAadjusted$Species)))),cex=(RAadjusted$value)*1.2,pch=19)
axis(3,at=1:11,labels=levels(RAadjusted$variable))
axis(2,at=1:9,labels=rev(levels(RAadjusted$Species)),las=1,font=3)
sizes <- factor(RAadjusted$Species, levels=rev(levels(RAadjusted$Species)))
dev.off()

#RABS (small)
pdf(file="../figures/RABSsurveyS.pdf",width=9,height = 4)
par(mar=c(2.1,11.1,4.1,2.1))
plot(as.numeric(RAadjusted.3$variable),as.numeric(RAadjusted.3$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,9.5))
points(as.numeric(RAadjusted.3$variable),as.numeric(factor(RAadjusted.3$Species, levels=rev(levels(RAadjusted.3$Species)))),cex=(RAadjusted.3$value)*0.7,pch=19)
axis(3,at=1:11,labels=levels(RAadjusted.3$variable))
axis(2,at=1:9,labels=rev(levels(RAadjusted.3$Species)),las=1,font=3)
dev.off()


  #eDNA presence absence
pdf(file="../figures/eDNAsurveyPA.pdf",width=9,height = 5.5)
par(mar=c(2.1,11.1,4.1,2.1))
plot(as.numeric(RAadjusted.2$variable),as.numeric(RAadjusted.2$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,9.5))
points(as.numeric(RAadjusted.2$variable),as.numeric(factor(RAadjusted.2$Species, levels=rev(levels(RAadjusted.2$Species)))),cex=(RAadjusted.2$value)*3,pch=19)
axis(3,at=1:11,labels=levels(RAadjusted.2$variable))
axis(2,at=1:9,labels=rev(levels(RAadjusted.2$Species)),las=1,font=3)
dev.off()

#eDNA semi-quant
pdf(file="../figures/eDNAsurveySQ.pdf",width=9,height = 5.5)
par(mar=c(2.1,11.1,4.1,2.1))
plot(as.numeric(RAadjusted.3$variable),as.numeric(RAadjusted.3$Species),type="n",ylab = "",xlab = "",xaxt="n",yaxt="n",ylim=c(0.5,9.5))
points(as.numeric(RAadjusted.3$variable),as.numeric(factor(RAadjusted.3$Species, levels=rev(levels(RAadjusted.3$Species)))),cex=(RAadjusted.3$value)*0.7,pch=19)
axis(3,at=1:11,labels=levels(RAadjusted.3$variable))
axis(2,at=1:9,labels=rev(levels(RAadjusted.3$Species)),las=1,font=3)
dev.off()

#Is biomass of individuals related to eDNA amount?
summary(lm(RAadjusted$value~RAadjusted.3$value))
summary(lm(test$var1~test$var2))
cor.test(test$var1, test$var2, method="spearman")


##Lets create a subset of data for the natural/nonnatureal questions

##COI
cleanedCOI <- read.csv("../cleaned/Cleaned.Leraylulu.unoise3.csv",row.names = 1)
#required sites
sites <- c("HB","HN","SY","SN","MB","MN","KN","NK","PE","CN","RB","RN")
#new data
cleanedCOI <- cleanedCOI[grep(paste(sites,collapse="|"),names(cleanedCOI))]
#remove PE.B
cleanedCOI <- cleanedCOI[-grep("PE.B",names(cleanedCOI))]

rCOI<-
  cleanedCOI%>%
  filter(rowSums(.)>0)%>%
  t(.)%>%
  rrarefy(.,min(rowSums(.)))%>%
  t(.)

MDSCOI <- metaMDS(t(rCOI),distance = "bray")
MDSCOIdat <-as.data.frame(MDSCOI$points)
palette(c("#af8dc3","#1b7837","#7fbf7b","#b35806","#2166ac","#67a9cf","#fee0b6","#762a83","#8c510a","#d8b365","#5ab4ac","#01665e"))
plot(MDSCOIdat,pch=16,cex=2, col=as.factor( substr(colnames(cleanedCOI),3,4)))
#text(MDSCOIdat$MDS1,MDSCOIdat$MDS2,labels=substr(colnames(cleanedCOI),3,5))

#Now lets do the OTU data
results <- as.data.frame(apply(as.data.frame(rCOI), 2,function(x) length(x[x>0])))
results$location <- substr(row.names(results),3,4)
resultsmean <- aggregate(results$`apply(as.data.frame(rCOI), 2, function(x) length(x[x > 0]))`,by=list(results$location),mean)

barplot(resultsmean$x[c(12,11,2,3,5,6,4,7,8,1,9,10)],col=as.factor(resultsmean$Group.1)[c(12,11,2,3,5,6,4,7,8,1,9,10)],names=resultsmean$Group.1[c(12,11,2,3,5,6,4,7,8,1,9,10)],ylab="mean number of OTUs")

plot(1:length(resultsmean$x[c(12,11,2,3,5,6,4,7,8,1,9,10)]),rep(5,12),pch=16,cex=5,col=as.factor(resultsmean$Group.1)[c(12,11,2,3,5,6,4,7,8,1,9,10)],xlim=c(0,13))
text(1:length(resultsmean$x[c(12,11,2,3,5,6,4,7,8,1,9,10)]),rep(4.6,12),labels=resultsmean$Group.1[c(12,11,2,3,5,6,4,7,8,1,9,10)])

#Plot the map

map("world", xlim=c(14,34), ylim=c(-36,-26), col="gray", fill=TRUE)
points(latlongdata$Long[latlongdata$sitecode %in% substr(colnames(cleanedCOI),3,4)],latlongdata$Lat[latlongdata$sitecode %in% substr(colnames(cleanedCOI),3,4)],pch=19,cex=1,col="darkred")
map.scale(ratio=FALSE,cex=0.5)
text(latlongdata$Long[latlongdata$ori=="W" & latlongdata$type=="m"]-0.8,latlongdata$Lat[latlongdata$ori=="W" & latlongdata$type=="m"],labels=latlongdata$sitecode[ latlongdata$type=="m"& latlongdata$ori=="W"],cex=1.5)
text(latlongdata$Long[latlongdata$ori=="S"& latlongdata$type=="m"],latlongdata$Lat[latlongdata$ori=="S"& latlongdata$type=="m"]-0.5,labels=latlongdata$sitecode[latlongdata$ori=="S"& latlongdata$type=="m"],cex=1.5)
text(latlongdata$Long[latlongdata$ori=="E"& latlongdata$type=="m"]+0.8,latlongdata$Lat[latlongdata$ori=="E"& latlongdata$type=="m"],labels=latlongdata$sitecode[latlongdata$ori=="E"& latlongdata$type=="m"],cex=1.5)




