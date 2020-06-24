### SAeDNA Analysis###
#Luke E. Holman
#16.05.2019
######################
#Packages####

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
library("RgoogleMaps")
library("reshape")

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


#Set A plalette
palette(c('#a6cee3','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#8dd3c7','#8dd3c7','#8dd3c7','#b15928'))



#### Data Cleaning####

files <- list.files("rawdata",full.names = T)
for (file in  files){
  
  rawdat <-read.csv(file=file)
  
  

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


#### Map ####

#First lets build a map with our samples 
pdf("figures/mapv1.pdf",width=7,height=5)
map("world", xlim=c(14,34), ylim=c(-36,-26), col="gray", fill=TRUE)
points(latlongdata$Long[latlongdata$type=="m"],latlongdata$Lat[latlongdata$type=="m"],pch=19,cex=1,col="darkred")
points(latlongdata$Long[latlongdata$type=="n"],latlongdata$Lat[latlongdata$type=="n"],pch=19,cex=0.5,col="darkblue")
map.scale(ratio=FALSE,cex=0.5)
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





### Q1 Patterns of OTU richness across the coast ####
#COI
#RawRichness
COI <- read.csv(file = "cleaned/Cleaned.COI_DADA.csv")

COI2 <- COI[,-1]
rownames(COI2) <- COI$X
rownames(COI2)

##This line rarefys and preserves rownames
rCOI<-t(rrarefy(t(COI2[rowSums(COI2)>0,]),min(colSums(COI2[rowSums(COI2)>0,]))))

COI3 <- as.data.frame(rCOI)
COI3[COI3>1] <- 1
COI.OTUs <- data.frame("ID"=names(COI3),"OTUs"= colSums(COI3))
COI.OTUs$site <- substr(COI.OTUs$ID,3,4)
COI.OTUs <- COI.OTUs[COI.OTUs$site!="DN" & COI.OTUs$ID!="L.PE.B", ]
COI.OTUs$site <- factor(COI.OTUs$site, levels = as.character(latlongdata$sitecode[sort(latlongdata$Lat)]))
COI.OTU2 <- aggregate(COI.OTUs$OTUs,list(site=COI.OTUs$site),mean)

pdf("figures/COI.OTUrichness.pdf",width = 10,height=6)
barplot(COI.OTU2$x,names.arg=COI.OTU2$site,col=as.numeric(latlongdata$type[latlongdata$sitecode %in% COI.OTU2$site]),main="COI",ylab="OTUs")
dev.off()

#ByCoastRichness (here we have chosen only to exmaine OTUs appearing in all three technical reps)
COI3 <- COI3[,colnames(COI3)!="L.PE.B" & colnames(COI3)!="L.DN1"]
substr(colnames(COI3),3,4)

singleCOI <- as.data.frame(t(rowsum(t(COI3),substr(colnames(COI3),3,4))))
singleCOI[singleCOI<2.5] <- 0
singleCOI[singleCOI>2.5] <- 1
singleCOI <- singleCOI[rowSums(singleCOI) > 0,]

#Durban Lagoon labelled DC for some reason - lets label it back!
colnames(singleCOI)[3] <- "DL"


pdf("figures/COIcoastalrichness.pdf",width = 4,height=6)
plot(as.numeric(latlongdata$PERMori[match(colnames(singleCOI),as.character(latlongdata$sitecode))]),colSums(singleCOI),pch=16,cex=1.5,
     col=as.factor(as.character(latlongdata$type[match(colnames(singleCOI),as.character(latlongdata$sitecode))])),
     xaxt="n",
     ylab="OTU Richness",
     xlab="",
     main="COI")
axis(1,1:3,c("West","South","East"))
dev.off()


#18S
#RawRichness
Z18S <- read.csv(file = "cleaned/Cleaned.ZHAN_DADA.csv")
Z18S2 <- Z18S[,-1]
rownames(Z18S2) <- Z18S$X
rZ18S<-t(rrarefy(t(Z18S2[rowSums(Z18S2)>0,]),min(colSums(Z18S2))))

Z18S3 <- as.data.frame(rZ18S)
Z18S3[Z18S3>1] <- 1
Z18S.OTUs <- data.frame("ID"=names(Z18S3),"OTUs"= colSums(Z18S3))
Z18S.OTUs$site <- substr(Z18S.OTUs$ID,3,4)
Z18S.OTUs <- Z18S.OTUs[Z18S.OTUs$site!="Z.DN1" & Z18S.OTUs$ID!="Z.PE.B", ]
Z18S.OTUs$site <- factor(Z18S.OTUs$site, levels = as.character(latlongdata$sitecode[sort(latlongdata$Lat)]))
Z18S.OTU2 <- aggregate(Z18S.OTUs$OTUs,list(site=Z18S.OTUs$site),mean)

pdf("figures/Z18S.OTUrichness.pdf",width = 10,height=6)
barplot(Z18S.OTU2$x,names.arg=Z18S.OTU2$site,col=as.numeric(latlongdata$type[latlongdata$sitecode %in%Z18S.OTU2$site]),main="18S",ylab="OTUs")
dev.off()

#ByCoastRichness

#ByCoastRichness (here we have chosen only to exmaine OTUs appearing in all three technical reps)
Z18S3 <- Z18S3[,colnames(Z18S3)!="Z.PE.B" & colnames(Z18S3 )!="Z.DN1"]
substr(colnames(Z18S3),3,4)

singleZ18S3 <- as.data.frame(t(rowsum(t(Z18S3),substr(colnames(Z18S3),3,4))))
singleZ18S3[singleZ18S3<2.5] <- 0
singleZ18S3[singleZ18S3>2.5] <- 1
singleZ18S3 <- singleZ18S3[rowSums(singleZ18S3) > 0,]

#Durban Lagoon labelled DC for some reason - lets label it back!
colnames(singleZ18S3)[3] <- "DL"

pdf("figures/18Scoastalrichness.pdf",width = 4,height=6)
plot(as.numeric(latlongdata$PERMori[match(colnames(singleZ18S3),as.character(latlongdata$sitecode))]),colSums(singleZ18S3),pch=16,cex=1.2,
     col=as.factor(as.character(latlongdata$type[match(colnames(singleZ18S3),as.character(latlongdata$sitecode))])),
     xaxt="n",
     ylab="OTU Richness",
     xlab="",
     main="18S")
axis(1,1:3,c("West","South","East"))
dev.off()


#ProtK

#RawRichness
ProtK <- read.csv(file = "cleaned/Cleaned.PROK_DADA.csv")
ProtK2 <- ProtK[,-1]
rownames(ProtK2) <- ProtK$X
rProtK<-t(rrarefy(t(ProtK2[rowSums(ProtK2)>0,]),min(colSums(ProtK2))))

ProtK3 <- as.data.frame(rProtK)
ProtK3[ProtK3>1] <- 1
ProtK.OTUs <- data.frame("ID"=names(ProtK3),"OTUs"= colSums(ProtK3))
ProtK.OTUs$site <- substr(ProtK.OTUs$ID,1,2)
ProtK.OTUs <- ProtK.OTUs[ProtK.OTUs$site!="DN" & ProtK.OTUs$ID!="L.PE.B", ]
ProtK.OTUs$site <- factor(ProtK.OTUs$site, levels = as.character(latlongdata$sitecode[sort(latlongdata$Lat)]))
ProtK.OTU2 <- aggregate(ProtK.OTUs$OTUs,list(site=ProtK.OTUs$site),mean)

pdf("figures/ProtK.OTUrichness.pdf",width = 10,height=6)
barplot(ProtK.OTU2$x,names.arg=ProtK.OTU2$site,col=as.numeric(latlongdata$type[latlongdata$sitecode %in% ProtK.OTU2$site]),main="Prokaryotes",ylab="OTUs")
dev.off()

#ByCoast Richness
ProtK3 <- ProtK3[,colnames(ProtK3)!="PE.B" & colnames(ProtK3)!="DN1"]


singleProtK <- as.data.frame(t(rowsum(t(ProtK3),substr(colnames(ProtK3),1,2))))
singleProtK[singleProtK<2.5] <- 0
singleProtK[singleProtK>2.5] <- 1
singleProtK <- singleProtK[rowSums(singleProtK) > 0,]

#Durban Lagoon labelled DC for some reason - lets label it back!
colnames(singleProtK)[3] <- "DL"

pdf("figures/ProtKcoastalrichness.pdf",width = 4,height=6)
plot(as.numeric(latlongdata$PERMori[match(colnames(singleProtK),as.character(latlongdata$sitecode))]),colSums(singleProtK),pch=16,cex=1.2,
     col=as.factor(as.character(latlongdata$type[match(colnames(singleProtK),as.character(latlongdata$sitecode))])),
     xaxt="n",
     ylab="OTU Richness",
     xlab="",
     main="ProtK")
axis(1,1:3,c("West","South","East"))
dev.off()

 #Quick 2 way ANOVA - assumptions totally violated! 
data <- singleCOI
data <- data[colnames(data)!="DL"]

model <- aov(colSums(data)~latlongdata$PERMori[match(colnames(data),as.character(latlongdata$sitecode))]+
             latlongdata$type[match(colnames(data),as.character(latlongdata$sitecode))])
plot(model)
summary(model)
TukeyHSD(model)


#Lets compare to Awad2002
awad <- read.csv("metadata/awad.richness.csv")

pdf("figures/awad.richness.pdf",height = 8,width = 7)
par(mfrow=c(3,2))

#COI
plot(awad$Richness[match(colnames(singleCOI),awad$SiteID)],colSums(singleCOI),
     main="COI stringent",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)
plot(awad$Richness[match(COI.OTU2$site,awad$SiteID)],COI.OTU2$x,
     main="COI relaxed",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)

#18S
plot(awad$Richness[match(colnames(singleZ18S3),awad$SiteID)],colSums(singleZ18S3),
     main="18S stringent",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)
plot(awad$Richness[match(Z18S.OTU2$site,awad$SiteID)],Z18S.OTU2$x,
     main="18S relaxed",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)

#16S
plot(awad$Richness[match(colnames(singleProtK),awad$SiteID)],colSums(singleProtK),
     main="16S stringent",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)
plot(awad$Richness[match(ProtK.OTU2$site,awad$SiteID)],ProtK.OTU2$x,
     main="18S relaxed",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)

dev.off()


#### Q2 is there a difference in beta diversity jaccard and bray-curtis between artificial/natural and west/central/east ####
##For all these analyses 

#Part 1 - nMDS and MDS

##COI
##First we need to collapse technical replicates as these are not ecological replicates
#get rid of these
rCOI2 <- rCOI[,colnames(rCOI)!="L.PE.B" & colnames(rCOI)!="L.DN1" ] 

#use this juicy sapply to average the number of reads across the replicates 
rCOI3 <- sapply(unique(substr(colnames(rCOI2),3,4)),function(x) rowMeans(rCOI2[,substr(colnames(rCOI2),3,4)==x]))

# lets get rid of DL because it is freshwater 
rCOI3 <- rCOI3[,colnames(rCOI3)!="DC"]
rCOI3 <- rCOI3[rowSums(rCOI3) > 0,]


##AllData
pdf("figures/COIbdiv.pdf",height=8,width=10)
par(mfrow=c(2,2))
##PCoA
rCOI.dist <- vegdist(t(rCOI3), "bray")
pcoa <- cmdscale(rCOI.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="COI - PCoA - Bray")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))

rCOI.dist <- vegdist(t(rCOI3), "jaccard")
pcoa <- cmdscale(rCOI.dist)
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="COI - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


#nMDS
test <- metaMDS(t(rCOI3),"bray")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="COI - nMDS - Bray",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

test <- metaMDS(t(rCOI3),"jaccard")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="COI - nMDS - Jaccard",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

dev.off()




##SmallData
rCOI4 <- rCOI3[,match(c("SY","SN","HB","HN","MB","MN","KN","NK","PE","CN","RB","RN"),colnames(rCOI3))]

pdf("figures/COIbdiv.reducedDATA.pdf",height=8,width=10)
par(mfrow=c(2,2))
##PCoA
rCOI.dist <- vegdist(t(rCOI4), "bray")
pcoa <- cmdscale(rCOI.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="COI - PCoA - Bray")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))

rCOI.dist <- vegdist(t(rCOI4), "jaccard")
pcoa <- cmdscale(rCOI.dist)
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="COI - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


#nMDS
test <- metaMDS(t(rCOI4),"bray")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="COI - nMDS - Bray",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

test <- metaMDS(t(rCOI4),"jaccard")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="COI - nMDS - Jaccard",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

dev.off()

#18S

##First we need to collapse technical replicates as these are not ecological replicates
#get rid of these
rZ18S2 <- rZ18S[,colnames(rZ18S)!="Z.PE.B" & colnames(rZ18S)!="Z.DN1" ] 

#use this juicy sapply to average the number of reads across the replicates 
rZ18S3<- sapply(unique(substr(colnames(rZ18S2),3,4)),function(x) rowMeans(rZ18S2[,substr(colnames(rZ18S2),3,4)==x]))

# lets get rid of DL because it is freshwater 
rZ18S3 <- rZ18S3[,colnames(rZ18S3)!="DC"]
rZ18S3 <- rZ18S3[rowSums(rZ18S3) > 0,]

##AllData
pdf("figures/Z18Sbdiv.pdf",height=8,width=10)
par(mfrow=c(2,2))
##PCoA
rZ18S.dist <- vegdist(t(rZ18S3), "bray")
pcoa <- cmdscale(rZ18S.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="18S - PCoA - Bray")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))

rZ18S.dist <- vegdist(t(rZ18S3), "jaccard")
pcoa <- cmdscale(rZ18S.dist)
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="18S - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


#nMDS
test <- metaMDS(t(rZ18S3),"bray")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="18S - nMDS - Bray",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

test <- metaMDS(t(rZ18S3),"jaccard")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="18S - nMDS - Jaccard",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

dev.off()


##SmallData
rZ18S4 <- rZ18S3[,match(c("SY","SN","HB","HN","MB","MN","KN","NK","PE","CN","RB","RN"),colnames(rZ18S3))]

pdf("figures/Z18Sbdiv.reducedDATA.pdf",height=8,width=10)
par(mfrow=c(2,2))
##PCoA
rZ18S.dist <- vegdist(t(rZ18S4), "bray")
pcoa <- cmdscale(rZ18S.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="18S - PCoA - Bray")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))

rZ18S.dist <- vegdist(t(rZ18S4), "jaccard")
pcoa <- cmdscale(rZ18S.dist)
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="18S - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


#nMDS
test <- metaMDS(t(rZ18S4),"bray")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="18S - nMDS - Bray",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

test <- metaMDS(t(rZ18S4),"jaccard")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="18S - nMDS - Jaccard",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

dev.off()


#Finally Prokarotes
##First we need to collapse technical replicates as these are not ecological replicates
#get rid of these
rProtK2 <- rProtK[,colnames(rProtK)!="PE.B" & colnames(rProtK)!="DN1" ] 

#use this juicy sapply to average the number of reads across the replicates 
rProtK3<- sapply(unique(substr(colnames(rProtK2),1,2)),function(x) rowMeans(rProtK2[,substr(colnames(rProtK2),1,2)==x]))

# lets get rid of DL because it is freshwater 
rProtK3 <- rProtK3[,colnames(rProtK3)!="DC"]
rProtK3 <- rProtK3[rowSums(rProtK3) > 0,]

##AllData
pdf("figures/ProtKbdiv.pdf",height=8,width=10)
par(mfrow=c(2,2))
##PCoA
rProtK.dist <- vegdist(t(rProtK3), "bray")
pcoa <- cmdscale(rProtK.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="16S - PCoA - Bray")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))

rProtK.dist <- vegdist(t(rProtK3), "jaccard")
pcoa <- cmdscale(rProtK.dist)
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="16S - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


#nMDS
test <- metaMDS(t(rProtK3),"bray")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="16S - nMDS - Bray",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

test <- metaMDS(t(rProtK3),"jaccard")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="16S - nMDS - Jaccard",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

dev.off()


##SmallData
rProtK4 <- rProtK3[,match(c("SY","SN","HB","HN","MB","MN","KN","NK","PE","CN","RB","RN"),colnames(rProtK3))]

pdf("figures/ProtKbdiv.reducedDATA.pdf",height=8,width=10)
par(mfrow=c(2,2))
##PCoA
rProtK.dist <- vegdist(t(rProtK4), "bray")
pcoa <- cmdscale(rProtK.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="16S - PCoA - Bray")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))

rProtK.dist <- vegdist(t(rProtK4), "jaccard")
pcoa <- cmdscale(rProtK.dist)
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="16S - PCoA - Jaccard")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


#nMDS
test <- metaMDS(t(rProtK4),"bray")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="16S - nMDS - Bray",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

test <- metaMDS(t(rProtK4),"jaccard")
#plot(test$points,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16)
plot(test$points,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,main="16S - nMDS - Jaccard",cex=3)
text(test$points[,1],test$points[,2],labels=rownames(test$points))

dev.off()



#### Now lets test these patterns with PERMANOVA

##Lets try out a bunch of things

c("rCOI3")

data <- rCOI3[,colnames(rProtK3)!= "SB" & colnames(rProtK3)!= "SM"  ]
method1 <- "bray"

#We can write a shitty function that pulls out P values and dispertion tests

tester <- function(data,method1){

#Get natural/nonnatural
sitetype <- latlongdata$type[match(colnames(data),latlongdata$sitecode)]
#Coast
coast <- latlongdata$PERMori[match(colnames(data),latlongdata$sitecode)]

##Get those data out
loop.coast.p <-adonis(t(data)~coast*sitetype,method=method1)$aov.tab["coast",6]
loop.type.p <-adonis(t(data)~coast*sitetype,method=method1)$aov.tab["sitetype",6]
Disper.C <- anova(betadisper(vegdist(t(data),method=method1),coast))$`Pr(>F)`[1]
Disper.S <- anova(betadisper(vegdist(t(data),method=method1),sitetype))$`Pr(>F)`[1]
Disper.CS <- anova(betadisper(vegdist(t(data),method=method1),paste0(coast,sitetype)))$`Pr(>F)`[1]

return(round(c(loop.coast.p,loop.type.p,Disper.C,Disper.S,Disper.CS),digits = 3))
}

PERM.results <- as.data.frame(matrix(nrow=18,ncol = 6))
colnames(PERM.results) <- c("Data.Method","Coast.P","SiteType.P","Disp.Coast.P","Disp.SiteType.P","Disp.Both.P")
PERM.results[1,] <- c("COI.All.Bray",tester(rCOI3,"bray"))
PERM.results[2,] <- c("COI.All.Jacc",tester(rCOI3,"jaccard"))
PERM.results[3,] <- c("COI.Clean.Bray",tester(rCOI3[,colnames(rCOI3)!= "SB" & colnames(rCOI3)!= "SM"  ],"bray"))
PERM.results[4,] <- c("COI.Clean.Jacc",tester(rCOI3[,colnames(rCOI3)!= "SB" & colnames(rCOI3)!= "SM"  ],"jaccard"))
PERM.results[5,] <- c("COI.Neat.Bray",tester(rCOI4,"bray"))
PERM.results[6,] <- c("COI.Neat.Jacc",tester(rCOI4,"jaccard"))
PERM.results[7,] <- c("18S.All.Bray",tester(rZ18S3,"bray"))
PERM.results[8,] <- c("18S.All.Jacc",tester(rZ18S3,"jaccard"))
PERM.results[9,] <- c("18S.Clean.Bray",tester(rZ18S3[,colnames(rZ18S3)!= "SB" & colnames(rZ18S3)!= "SM"  ],"bray"))
PERM.results[10,] <- c("18S.Clean.Jacc",tester(rZ18S3[,colnames(rZ18S3)!= "SB" & colnames(rZ18S3)!= "SM"  ],"jaccard"))
PERM.results[11,] <- c("18S.Neat.Bray",tester(rZ18S4,"bray"))
PERM.results[12,] <- c("18S.Neat.Jacc",tester(rZ18S4,"jaccard"))
PERM.results[13,] <- c("ProK.All.Bray",tester(rProtK3,"bray"))
PERM.results[14,] <- c("ProK.All.Jacc",tester(rProtK3,"jaccard"))
PERM.results[15,] <- c("ProK.Clean.Bray",tester(rProtK3[,colnames(rProtK3)!= "SB" & colnames(rProtK3)!= "SM"  ],"bray"))
PERM.results[16,] <- c("ProK.Clean.Jacc",tester(rProtK3[,colnames(rProtK3)!= "SB" & colnames(rProtK3)!= "SM"  ],"jaccard"))
PERM.results[17,] <- c("ProK.Neat.Bray",tester(rProtK4,"bray"))
PERM.results[18,] <- c("ProK.Neat.Jacc",tester(rProtK4,"jaccard"))

write.csv(PERM.results,file = "PERM.out.csv")


##Lets subset the OTUs based on what remains 

OTUs <- seqinr::read.fasta(file="/Volumes/BackUp2/Ext.HT.Sequence.data/2018.SAeDNA/SAeDNA1.1/5.OTUs/dada2.fasta")
newOTUs <- OTUs[names(OTUs) %in% rownames(rCOI3)]
seqinr::write.fasta(newOTUs,names=names(newOTUs),file="~/Desktop/DADA.OTUs.COI.fasta")

#Pull in the taxonomy 
COI.taxa <- read.csv("Taxonomy/HighConf.COI.csv",header=F)


#replot some stuff with all taxonomy 

COI.taxasubset <- COI3[rownames(COI3) %in% COI.taxa$V1,]
COI.OTUs.taxa <- data.frame("ID"=names(COI.taxasubset),"OTUs"= colSums(COI.taxasubset))
COI.OTUs.taxa$site <- substr(COI.OTUs.taxa$ID,3,4)
COI.OTUs.taxa <- COI.OTUs.taxa[COI.OTUs.taxa$site!="DN" & COI.OTUs.taxa$ID!="L.PE.B", ]
COI.OTUs.taxa$site <- factor(COI.OTUs.taxa$site, levels = as.character(latlongdata$sitecode[sort(latlongdata$Lat)]))
COI.OTUs.taxa <- aggregate(COI.OTUs.taxa$OTUs,list(site=COI.OTUs.taxa$site),mean)
barplot(COI.OTUs.taxa$x,names.arg=COI.OTUs.taxa$site,col=as.numeric(latlongdata$type[latlongdata$sitecode %in% COI.OTUs.taxa$site]),main="COI",ylab="OTUs")

plot(awad$Richness[match(COI.OTUs.taxa$site,awad$SiteID)],COI.OTUs.taxa$x,
     main="COI.taxonomy",
     ylab="eDNA OTU Richness",
     xlab="Awad Richness",
     pch=16)


datasubset <- COI3[rownames(COI3) %in% COI.taxa$V1[COI.taxa$V6=="Ascidiacea"],]

COI.OTUs.taxa <- data.frame("ID"=names(datasubset),"OTUs"= colSums(datasubset))
COI.OTUs.taxa$site <- substr(COI.OTUs.taxa$ID,3,4)
COI.OTUs.taxa <- COI.OTUs.taxa[COI.OTUs.taxa$site!="DN" & COI.OTUs.taxa$ID!="L.PE.B", ]
COI.OTUs.taxa$site <- factor(COI.OTUs.taxa$site, levels = as.character(latlongdata$sitecode[sort(latlongdata$Lat)]))
COI.OTUs.taxa <- aggregate(COI.OTUs.taxa$OTUs,list(site=COI.OTUs.taxa$site),mean)
barplot(COI.OTUs.taxa$x,names.arg=COI.OTUs.taxa$site,col=as.numeric(latlongdata$type[latlongdata$sitecode %in% COI.OTUs.taxa$site]),main="COI",ylab="OTUs")


data <- rCOI3[rownames(rCOI3) %in% COI.taxa$V1,]

sitetype <- latlongdata$type[match(colnames(data),latlongdata$sitecode)]
#Coast
coast <- latlongdata$PERMori[match(colnames(data),latlongdata$sitecode)]

##Get those data out
adonis(t(data)~coast*sitetype,method=method1)
anova(betadisper(vegdist(t(data),method=method1),paste0(coast,sitetype)))

data.dist <- vegdist(t(data), "bray")
pcoa <- cmdscale(data.dist)
#plot(pcoa,col=latlongdata$PERMori[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3)
#text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))
plot(pcoa,col=latlongdata$type[match(rownames(pcoa),latlongdata$sitecode)],pch=16,cex=3,main="COI.subset")
text(pcoa[,1],pcoa[,2],labels=rownames(pcoa))


####Now lets make a distance by beta dissimlairty chart ####
data <- rCOI3

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
pdf("figures/Distance.Similarity.1.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "lightblue")
dev.off()

#Now lets make a character that desrbies the type of comparison
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


plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise$type))
legend("topright",col=1:3,legend = levels(as.factor(pairwise$type)),pch=16)
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],log10(1-pairwise$value),ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= as.factor(pairwise$type))

plot(lm(log10(1-pairwise$value)~Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]))
summary(lm(log10(1-pairwise$value)~Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]))

summary(lm(log10(1-pairwise$value)~Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]*as.factor(pairwise$type)))

####Taxonomy






#####CODE BASEMENT####


#Nuts and bolts
round(MDS$eig*100/sum(MDS$eig),1)


vegdist(t(rCOI3[,colnames(rCOI3) %in% c("SY","TB","HB")]),"bray")
vegdist(t(rCOI3[,colnames(rCOI3) %in% c("MB","KN","PE")]),"bray")
vegdist(t(rCOI3[,colnames(rCOI3) %in% c("BR","EL","DU","RB")]),"bray")



#We rearrange things as so

rCOI3 <- rCOI3[,na.exclude(match(as.character(latlongdata$sitecode),colnames(rCOI3)))]

#Now we produce the beta values as so. 

test <- as.matrix(vegdist(t(test),"jaccard",upper = T))
test[test==0] <- NA

heatmap(test,Colv = NA, Rowv = NA,scale="none")





