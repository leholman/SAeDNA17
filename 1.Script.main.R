#############################################
####==== South African eDNA Analysis ====####
####==== Luke E. Holman====30.03.2020====####
#############################################

###Script 1 - Analysis Script###


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
library("gdata")
library(adespatial)
library(spdep)
library(ade4)
library(prettymapr)

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


#Read in data
rCOI <- read.csv("cleaned/rarefied.COI.csv",row.names = 1)
r18S <- read.csv("cleaned/rarefied.18S.csv",row.names = 1)
rProK <- read.csv("cleaned/rarefied.ProK.csv",row.names = 1)

#Read in taxa
COIassignments <- read.csv("Taxonomy/CleanedTaxonomy/COI.blast.csv")
z18Sassignments <- read.csv("Taxonomy/CleanedTaxonomy/18S.blast.csv")
ProKtaxa <- read.csv("Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.csv")
ProKtaxaALL <- read.csv("Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.ALL.csv")


####====1.0 Maps====####

#Subset the data for only the sites in our data
smallLatLong <- latlongdata[latlongdata$sitecode %in% names(rCOI),]

pdf("figures/Figure1/mapv2.pdf",width=9,height=6.5)
m<- map("world", xlim=c(16,33), ylim=c(-36,-28), col="gray", fill=TRUE,mar=c(4, 4, 4, 4))

#overview ecoregions
rect(15,-38,20,-27,border = NA,col="#78C1E1")
rect(20,-38,26.5,-27,border = NA,col="#139F73")
rect(26.5,-38,34,-27,border = NA,col="#E69F03")

m<- map("world", xlim=c(16,33), ylim=c(-36,-28), col="gray", fill=TRUE,mar=c(4, 4, 4, 4),add=TRUE)
addnortharrow("bottomright",scale=0.75)
addscalebar()


xat <- pretty(m$range[1:2],n = 10)
xlab <- parse(text=degreeLabelsEW(xat))

yat <- pretty(m$range[3:4])
ylab <- parse(text=degreeLabelsNS(yat))

box()
axis(2, las=TRUE, at=yat, labels=ylab)
axis(3, at=xat, labels=xlab)
axis(4, las=TRUE, at=yat, labels=ylab)

points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=1,cex=3,lwd=4,col="white")
points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=1,cex=3,lwd=2,col="darkblue")

points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=4,cex=2,lwd=4,col="white")
points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=4,cex=2,lwd=2,col="darkred")

#Greyscale Os
#points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=16,cex=2.3,col="black")
#points(smallLatLong$Long[smallLatLong$type=="n"],smallLatLong$Lat[smallLatLong$type=="n"],pch=16,cex=2,col="white")

#points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=1,cex=2,lwd=4,col="white")
#points(smallLatLong$Long[smallLatLong$type=="m"],smallLatLong$Lat[smallLatLong$type=="m"],pch=1,cex=2,lwd=2,col="black")

dev.off()

  
####====2.0 Alpha Diversity====####

## Patterns of alpha  taxonomic diversity across the coast 
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

pdf("figures/Figure1/OTU.coast.Richness2.pdf",width = 9,height=3)
par(mar=c(2,4,1,1))
palette(c('#54B4E9','#F0E442','#CC79A7'))
plot(jitter(as.numeric(Coastal.alpha.tall$comb),0.8),
     Coastal.alpha.tall$value,
     type="n",
     xaxt="n",
     ylab="ASVs")
axis(1,at=1:9,labels=rep(c("COI","18S","16S"),3),cex=0.8)

#Add ecoregions
rect(0,0,3.5,1000,border = NA,col="#D3ECF6")
rect(3.5,0,6.5,1000,border = NA,col="#B6DDD3")
rect(6.5,0,10,1000,border = NA,col="#F9E2B4")

#Add points
points(jitter(as.numeric(Coastal.alpha.tall$comb),0.8),
       Coastal.alpha.tall$value,pch=16)

##Lets draw in those mean values
for (num in 1:9){
  run.mean <- mean(Coastal.alpha.tall$value[as.numeric(Coastal.alpha.tall$comb)==num])
  lines(c(num-0.2,num+0.2),c(run.mean,run.mean),lwd=2)
}
box()
dev.off()


####====3.0 Taxonomy====####
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

COI.taxa.count <- minAbundance(CountTable(COIassignments[match(rownames(rCOI),COIassignments$OTU),"phylum"],rCOI,output="Count"),minAbun= 0.02)
z18S.taxa.count <- minAbundance(CountTable(z18Sassignments[match(rownames(r18S),z18Sassignments$OTU),"phylum"],r18S,output="Count"),minAbun= 0.02)
ProK.taxa.count <- minAbundance(CountTable(as.character(ProKtaxa$Phylum),rProK,output="Count"),minAbun= 0.02)


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

#Now lets make an individual one for the big figure

COIplotdat <- prop.table(as.matrix(COI.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(COI.taxa.count)))]),2)[-1,]
z18Splotdat <- prop.table(as.matrix(z18S.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(z18S.taxa.count)))]),2)[-1,]
ProKplotdat <- prop.table(as.matrix(ProK.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.count)))]),2)[-1,]


pdf("figures/Figure1/COI.abun.pdf",width = 11, height = 2.3)
par(mfrow=c(1,1),mai = c(0.1, 0.6, 0.1,0),xpd=TRUE)
barplot(COIplotdat[dim(COIplotdat)[1]:1,],col=rev(getPalette(dim(COI.taxa.count)[1]-1)),axisnames=FALSE,cex.axis=1.35,las=1,
       #you can set the distance between the y axis and the bars by using the two below options, the xlim changes the distances
         xaxs = "i",xlim = c(-0.3, 22))
dev.off()

pdf("figures/Figure1/18S.abun.pdf",width = 11, height = 2.3)
par(mfrow=c(1,1),mai = c(0.1, 0.6, 0.1,0),xpd=TRUE)
barplot(z18Splotdat[dim(z18Splotdat)[1]:1,],col=rev(getPalette(dim(z18S.taxa.count)[1]-1)),axisnames=FALSE,cex.axis=1.35,las=1,
        xaxs = "i",xlim = c(-0.3, 22))
dev.off()

pdf("figures/Figure1/16S.abun.pdf",width = 11, height = 2.3)
par(mfrow=c(1,1),mai = c(0.1, 0.6, 0.1,0),xpd=TRUE)
barplot(ProKplotdat[dim(ProKplotdat)[1]:1,],col=rev(getPalette(dim(ProK.taxa.count)[1]-1)),axisnames=FALSE,cex.axis=1.35,las=1,
        xaxs = "i",xlim = c(-0.3, 22))
dev.off()


#If the current legend in fig.1 breaks we can use the below template to generate a new ones
pdf("figures/Figure1/COI.legend.pdf",width = 2.8, height = 3.4)
par(mai = c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
legend(2,10,rownames(COIplotdat),fill=getPalette(dim(COI.taxa.count)[1]-1),cex=1.2,bty = "n",y.intersp=0.75)
dev.off()

pdf("figures/Figure1/18S.legend.pdf",width = 2.8, height = 3.4)
par(mai = c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
legend(2,10,rownames(z18Splotdat),fill=getPalette(dim(z18S.taxa.count)[1]-1),cex=1.2,bty = "n",y.intersp=0.75)
dev.off()

pdf("figures/Figure1/16S.legend.pdf",width = 2.8, height = 3.4)
par(mai = c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
legend(2,10,rownames(ProKplotdat),fill=getPalette(dim(ProK.taxa.count)[1]-1),cex=1.2,bty = "n",y.intersp=0.75)
dev.off()

# Now lets make a nice version of the above for the supplement

COIplotdat <- prop.table(as.matrix(COI.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(COI.taxa.abun)))]),2)[-1,]
z18Splotdat <- prop.table(as.matrix(z18S.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(z18S.taxa.abun)))]),2)[-1,]
ProKplotdat <- prop.table(as.matrix(ProK.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.abun)))]),2)[-1,]


pdf("figures/SupplAbundphyla.pdf",width = 9, height = 6)
par(mfrow=c(3,1),mai = c(0.4, 0.4, 0.1,2),xpd=TRUE)
barplot(COIplotdat[dim(COIplotdat)[1]:1,],col=rev(getPalette(dim(COI.taxa.count)[1]-1)),axisnames=FALSE,cex.axis=1.35,las=1,xaxs = "i",xlim = c(-0.3, 22))
legend(22,0.8,rownames(COIplotdat),fill=getPalette(dim(COI.taxa.count)[1]-1),cex=1,bty = "n",y.intersp=0.75)
barplot(z18Splotdat[dim(z18Splotdat)[1]:1,],col=rev(getPalette(dim(z18S.taxa.count)[1]-1)),axisnames=FALSE,cex.axis=1.35,las=1,xaxs = "i",xlim = c(-0.3, 22))
legend(22,1,rownames(z18Splotdat),fill=getPalette(dim(z18S.taxa.count)[1]-1),cex=1,bty = "n",y.intersp=0.75)
barplot(ProKplotdat[dim(ProKplotdat)[1]:1,],col=rev(getPalette(dim(ProK.taxa.count)[1]-1)),cex.axis=1.35,las=1,xaxs = "i",xlim = c(-0.3, 22))
legend(22,1,rownames(ProKplotdat),fill=getPalette(dim(ProK.taxa.count)[1]-1),cex=1,bty = "n",y.intersp=0.75)
dev.off()



####====4.0 Beta Diversity====####

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
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/COI.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordiellipse(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            kind = "sd",
            draw = "polygon",
            lty=0,
            col=1:3)
points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m"],
       pch=16,cex=3)

points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n"],
       pch=17,cex=3)
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

#Here we are using smoothed GAMs in the function ordisurf across the nMDS space
##Some thoughts on nMDS ordisurf - https://www.fromthebottomoftheheap.net/2011/06/10/what-is-ordisurf-doing/

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(rCOI),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(rCOI),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(rCOI),envdat$sitecode)])


#Ordisurf - temp
pdf("figures/COI.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - SSS
pdf("figures/COI.SSS.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.SSS<- ordisurf(nMDS,envmatrix$SSS,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()
sink("model.output/smoothedGAM/COI.SSS.txt")
summary(surfout.SSS)
sink()

#Ordisurf - ChlA
pdf("figures/COI.Chl.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.Chl<- ordisurf(nMDS,envmatrix$chla_3yrmean,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.Chl.txt")
summary(surfout.Chl)
sink()


#Ordisurf - impact
pdf("figures/COI.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.impact.txt")
summary(surfout.ipt)
sink()


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

##Lets try and relate the environmental variables to the beta diversity 
#we are using a distance based redudnancy analysis (dbRDA)

#dbRDA
model.full <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean + temp, data = envmatrix)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.COI.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
sink()


model.3.1 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean + temp, data = envmatrix)
model.3.2 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean + temp, data = envmatrix)
model.3.3 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + impact + temp, data = envmatrix)
model.3.4 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean, data = envmatrix)

model.2.1 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + impact, data = envmatrix)
model.2.2 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean, data = envmatrix)
model.2.3 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ SSS + temp, data = envmatrix)
model.2.4 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean, data = envmatrix)
model.2.5 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ impact + temp, data = envmatrix)
model.2.6 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ chla_3yrmean + temp, data = envmatrix)

model.1.1 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~temp, data = envmatrix)
model.1.2 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~impact, data = envmatrix)
model.1.3 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~chla_3yrmean, data = envmatrix)
model.1.4 <- dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~SSS, data = envmatrix)

R2.1.1 <- unname(RsquareAdj(model.1.1)[[2]])*100
names(R2.1.1) <- "SST"
R2.1.2 <- unname(RsquareAdj(model.1.2)[[2]])*100
names(R2.1.2) <- "Impact"
R2.1.3 <- unname(RsquareAdj(model.1.3)[[2]])*100
names(R2.1.3) <- "Chl_a"
R2.1.4 <- unname(RsquareAdj(model.1.4)[[2]])*100
names(R2.1.4) <- "SSS"

R2.2.1 <- unname(RsquareAdj(model.2.1)[[2]])*100 
names(R2.2.1) <- "SSS&Impact"
R2.2.2 <- unname(RsquareAdj(model.2.2)[[2]])*100 
names(R2.2.2) <- "SSS&Chl_a"
R2.2.3 <- unname(RsquareAdj(model.2.3)[[2]])*100 
names(R2.2.3) <- "SSS&SST"
R2.2.4 <- unname(RsquareAdj(model.2.4)[[2]])*100 
names(R2.2.4) <- "Impact&Chl_a"
R2.2.5 <- unname(RsquareAdj(model.2.5)[[2]])*100 
names(R2.2.5) <- "Impact&SST"
R2.2.6 <- unname(RsquareAdj(model.2.6)[[2]])*100 
names(R2.2.6) <- "Chl_a&SST"

R2.3.1 <- unname(RsquareAdj(model.3.1)[[2]])*100 
names(R2.3.1) <- "Impact&Chl_a&SST"
R2.3.2 <- unname(RsquareAdj(model.3.2)[[2]])*100 
names(R2.3.2) <- "SSS&Chl_a&SST"
R2.3.3 <- unname(RsquareAdj(model.3.3)[[2]])*100 
names(R2.3.3) <- "SSS&Impact&SST"
R2.3.4 <- unname(RsquareAdj(model.3.4)[[2]])*100 
names(R2.3.4) <- "SSS&Impact&Chl_a"

R2.full <- unname(RsquareAdj(model.full)[[2]])*100 
names(R2.full) <- "SSS&Impact&Chl_a&SST"

expressionInput <- c(R2.full,R2.3.2,R2.3.1,R2.3.3,R2.3.4,R2.2.6,R2.2.3,R2.2.1,R2.2.5,R2.2.2,R2.2.4,R2.1.1,R2.1.4,R2.1.3,R2.1.2)
                     
R2COI <- expressionInput

#checking variance inflation factor
#vif.cca(dbrda(vegdist(t(rCOI),method = "jaccard",binary = TRUE)~ temp + impact + chla_3yrmean , data = envmatrix))


#Lets replace the upset plot with a bar plot since the order is junk
#barplot(expressionInput,col="black",names=F,yaxt="n",xaxt="n")
#axis(2,cex.axis=1.2,lwd=2)
#axis(1,at=-1:18,labels=F,cex=2,lwd=2,lwd.ticks=0)

pdf("figures/upset.dbRDA.COI.pdf",width = 9,height=5.5)
upset(fromExpression(expressionInput), 
      order.by = c("freq","degree"),
      decreasing=c("T","T"),
      point.size = 3.5, 
      line.size = 2,
      set_size.show =FALSE,
      show.numbers=FALSE,
      mainbar.y.label = "Variance Explained (%)",
      set_size.scale_max=1)
dev.off()

#Upset plots suck, lets try and make the bar ourselves
pdf("figures/Figure 2/bar.COI.pdf",width = 6,height = 4)
par(mar=c(3,3,1,1))
barplot(expressionInput,names=F,col="black",space=rep(0.5,15),yaxt='n',ylim=c(0,25))
axis(2,at=seq(0,25,5),labels=seq(0,25,5))
dev.off()


#how much variation explained per site
inertcomp(mods, proportional = TRUE,display = "sites")


#testing mantel and partial mantel tests

#Lets create a geogrphical distance matrix for all sites  
dist <- read.csv("metadata/distance.csv",row.names = 1)
dist <- as.matrix(dist)
upperTriangle(dist) <- lowerTriangle(dist,byrow = TRUE)
dist <- dist[match(colnames(rCOI),colnames(dist)),]
dist <- dist[,match(colnames(rCOI),colnames(dist))]
dist <- as.dist(dist)

#Lets create a false 2d world for calculations
dist.2d <- read.csv("metadata/distance.csv",row.names = 1)
dist.2d <- as.matrix(dist.2d)[,1]
dist.2d <- dist.2d[match(colnames(rCOI),names(dist.2d))]
dist.2d <- data.frame("Y"=rep(1,18),"X"=dist.2d)

#Lets create distance matrices for each env parameter
SST.dist <- dist(tempdat$TempAvr[match(colnames(rCOI),tempdat$Site.Code)])
SSS.dist <- dist(envdat$SSS[match(colnames(rCOI),envdat$sitecode)])
ChlA.dist <- dist(envdat$chla_3yrmean[match(colnames(rCOI),envdat$sitecode)])
Impact.dist <- dist(envdat$impact[match(colnames(rCOI),envdat$sitecode)])


#Mantel test
m.COI.SST <- mantel.randtest(vegdist(t(rCOI),method = "jaccard",binary = TRUE),SST.dist)
m.COI.SSS <- mantel.randtest(vegdist(t(rCOI),method = "jaccard",binary = TRUE),SSS.dist)
m.COI.ChlA <- mantel.randtest(vegdist(t(rCOI),method = "jaccard",binary = TRUE),ChlA.dist)
m.COI.impact <- mantel.randtest(vegdist(t(rCOI),method = "jaccard",binary = TRUE),Impact.dist)

##Here we run a partial mantel test
###but SST and distance are highly autocorrelated! Need to correct for this
p.COI.SST <- mantel.partial(vegdist(t(rCOI),method = "jaccard",binary = TRUE),SST.dist,dist)
p.COI.SSS <- mantel.partial(vegdist(t(rCOI),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.COI.ChlA <- mantel.partial(vegdist(t(rCOI),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.COI.impact <- mantel.partial(vegdist(t(rCOI),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.COI.SST.c <- msr(m.COI.SST,lw1, 10000)
m.COI.SSS.c <- msr(m.COI.SSS,lw1, 10000)
m.COI.ChlA.c <- msr(m.COI.ChlA,lw1, 10000)
m.COI.impact.c <- msr(m.COI.impact,lw1, 10000)


mantel.p.COI.out <- c(p.COI.SST$statistic,p.COI.SST$signif,
                     p.COI.SSS$statistic,p.COI.SSS$signif,
                     p.COI.ChlA$statistic,p.COI.ChlA$signif,
                     p.COI.impact$statistic,p.COI.impact$signif)
               
mantel.c.COI.out <- unname(c(m.COI.SST.c$obs-m.COI.SST.c$expvar["Expectation"],
                             m.COI.SST.c$pvalue,
                             m.COI.SSS.c$obs-m.COI.SSS.c$expvar["Expectation"],
                             m.COI.SSS.c$pvalue,
                             m.COI.ChlA.c$obs-m.COI.ChlA.c$expvar["Expectation"],
                             m.COI.ChlA.c$pvalue,
                             m.COI.impact.c$obs-m.COI.impact.c$expvar["Expectation"],
                             m.COI.impact.c$pvalue))

##This is how we might go about doing some variance partitioning
##varpart(vegdist(t(rProK),method = "jaccard",binary = TRUE),dist.2d$X, envmatrix)


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
points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m"],
       pch=16,cex=3)
points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n"],
       pch=17,cex=3)
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

##Here we are using smoothed GAMS to explain abiotic variation using the two nMDS axes


#Ordisurf - temp
pdf("figures/18S.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - SSS
pdf("figures/18S.SSS.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.SSS<- ordisurf(nMDS,envmatrix$SSS,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()
sink("model.output/smoothedGAM/18S.SSS.txt")
summary(surfout.SSS)
sink()

#Ordisurf - ChlA
pdf("figures/18S.Chl.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.Chl<- ordisurf(nMDS,envmatrix$chla_3yrmean,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.Chl.txt")
summary(surfout.Chl)
sink()


#Ordisurf - impact
pdf("figures/18S.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.impact.txt")
summary(surfout.ipt)
sink()



#Lets partition the diversity into turnover and nestedness 

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

##Lets try and relate the environmental variables to the beta diversity 
#we are using a distance based redudnancy analysis (dbRDA)

#dbRDA
model.full <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean + temp, data = envmatrix)
plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.18S.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
sink()


model.3.1 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean + temp, data = envmatrix)
model.3.2 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean + temp, data = envmatrix)
model.3.3 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + impact + temp, data = envmatrix)
model.3.4 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean, data = envmatrix)

model.2.1 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + impact, data = envmatrix)
model.2.2 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean, data = envmatrix)
model.2.3 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ SSS + temp, data = envmatrix)
model.2.4 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean, data = envmatrix)
model.2.5 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ impact + temp, data = envmatrix)
model.2.6 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~ chla_3yrmean + temp, data = envmatrix)

model.1.1 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~temp, data = envmatrix)
model.1.2 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~impact, data = envmatrix)
model.1.3 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~chla_3yrmean, data = envmatrix)
model.1.4 <- dbrda(vegdist(t(r18S),method = "jaccard",binary = TRUE)~SSS, data = envmatrix)

R2.1.1 <- unname(RsquareAdj(model.1.1)[[2]])*100
names(R2.1.1) <- "SST"
R2.1.2 <- unname(RsquareAdj(model.1.2)[[2]])*100
names(R2.1.2) <- "Impact"
R2.1.3 <- unname(RsquareAdj(model.1.3)[[2]])*100
names(R2.1.3) <- "Chl_a"
R2.1.4 <- unname(RsquareAdj(model.1.4)[[2]])*100
names(R2.1.4) <- "SSS"

R2.2.1 <- unname(RsquareAdj(model.2.1)[[2]])*100 
names(R2.2.1) <- "SSS&Impact"
R2.2.2 <- unname(RsquareAdj(model.2.2)[[2]])*100 
names(R2.2.2) <- "SSS&Chl_a"
R2.2.3 <- unname(RsquareAdj(model.2.3)[[2]])*100 
names(R2.2.3) <- "SSS&SST"
R2.2.4 <- unname(RsquareAdj(model.2.4)[[2]])*100 
names(R2.2.4) <- "Impact&Chl_a"
R2.2.5 <- unname(RsquareAdj(model.2.5)[[2]])*100 
names(R2.2.5) <- "Impact&SST"
R2.2.6 <- unname(RsquareAdj(model.2.6)[[2]])*100 
names(R2.2.6) <- "Chl_a&SST"

R2.3.1 <- unname(RsquareAdj(model.3.1)[[2]])*100 
names(R2.3.1) <- "Impact&Chl_a&SST"
R2.3.2 <- unname(RsquareAdj(model.3.2)[[2]])*100 
names(R2.3.2) <- "SSS&Chl_a&SST"
R2.3.3 <- unname(RsquareAdj(model.3.3)[[2]])*100 
names(R2.3.3) <- "SSS&Impact&SST"
R2.3.4 <- unname(RsquareAdj(model.3.4)[[2]])*100 
names(R2.3.4) <- "SSS&Impact&Chl_a"

R2.full <- unname(RsquareAdj(model.full)[[2]])*100 
names(R2.full) <- "SSS&Impact&Chl_a&SST"

expressionInput <- c(R2.full,R2.3.2,R2.3.1,R2.3.3,R2.3.4,R2.2.6,R2.2.3,R2.2.1,R2.2.5,R2.2.2,R2.2.4,R2.1.1,R2.1.4,R2.1.3,R2.1.2)

R218S <- expressionInput


pdf("figures/upset.dbRDA.18S.pdf",width = 9,height=5.5)
upset(fromExpression(expressionInput), 
      order.by = c("freq","degree"),
      decreasing=c("T","T"),
      point.size = 3.5, 
      line.size = 2,
      set_size.show =FALSE,
      show.numbers=FALSE,
      mainbar.y.label = "Variance Explained (%)",
      set_size.scale_max=1)
dev.off()


#Upset plots suck, lets try and make the bar ourselves
pdf("figures/Figure 2/bar.18S.pdf",width = 6,height = 4)
par(mar=c(3,3,1,1))
barplot(expressionInput,names=F,col="black",space=rep(0.5,15),ylim=c(0,20),yaxt="n")
axis(2,at=seq(0,20,5),labels=seq(0,20,5))
dev.off()


#Mantel test
m.18S.SST <- mantel.randtest(vegdist(t(r18S),method = "jaccard",binary = TRUE),SST.dist)
m.18S.SSS <- mantel.randtest(vegdist(t(r18S),method = "jaccard",binary = TRUE),SSS.dist)
m.18S.ChlA <- mantel.randtest(vegdist(t(r18S),method = "jaccard",binary = TRUE),ChlA.dist)
m.18S.impact <- mantel.randtest(vegdist(t(r18S),method = "jaccard",binary = TRUE),Impact.dist)

##Here we run a partial mantel test
###but SST and distance are highly autocorrelated! Need to correct for this
p.18S.SST <- mantel.partial(vegdist(t(r18S),method = "jaccard",binary = TRUE),SST.dist,dist)
p.18S.SSS <- mantel.partial(vegdist(t(r18S),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.18S.ChlA <- mantel.partial(vegdist(t(r18S),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.18S.impact <- mantel.partial(vegdist(t(r18S),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.18S.SST.c <- msr(m.18S.SST,lw1, 10000)
m.18S.SSS.c <- msr(m.18S.SSS,lw1, 10000)
m.18S.ChlA.c <- msr(m.18S.ChlA,lw1, 10000)
m.18S.impact.c <- msr(m.18S.impact,lw1, 10000)


mantel.p.18S.out <- c(p.18S.SST$statistic,p.18S.SST$signif,
                      p.18S.SSS$statistic,p.18S.SSS$signif,
                      p.18S.ChlA$statistic,p.18S.ChlA$signif,
                      p.18S.impact$statistic,p.18S.impact$signif)

mantel.c.18S.out <- unname(c(m.18S.SST.c$obs-m.18S.SST.c$expvar["Expectation"],
                             m.18S.SST.c$pvalue,
                             m.18S.SSS.c$obs-m.18S.SSS.c$expvar["Expectation"],
                             m.18S.SSS.c$pvalue,
                             m.18S.ChlA.c$obs-m.18S.ChlA.c$expvar["Expectation"],
                             m.18S.ChlA.c$pvalue,
                             m.18S.impact.c$obs-m.18S.impact.c$expvar["Expectation"],
                             m.18S.impact.c$pvalue))




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
points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m"],
       pch=16,cex=3)
points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n"],
       pch=17,cex=3)
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

##Here we are using smoothed GAMS to explain abiotic variation using the two nMDS axes

#Ordisurf - temp
pdf("figures/ProK.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/ProK.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - SSS
pdf("figures/ProK.SSS.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.SSS<- ordisurf(nMDS,envmatrix$SSS,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()
sink("model.output/smoothedGAM/ProK.SSS.txt")
summary(surfout.SSS)
sink()

#Ordisurf - ChlA
pdf("figures/ProK.Chl.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.Chl<- ordisurf(nMDS,envmatrix$chla_3yrmean,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/ProK.Chl.txt")
summary(surfout.Chl)
sink()


#Ordisurf - impact
pdf("figures/ProK.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/ProK.impact.txt")
summary(surfout.ipt)
sink()




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
par(mfrow=c(1,1))

#Let's exmaine the differences in a dendrogram (UPGMA)
#the below line makes some weights to reorder the leaves as much as possible to match coasts
weights <- match(test$labels,latlongdata$sitecode)
pdf("figures/dend.16S.pdf",width = 9,height=5.5)
plot(reorder(hclust(vegdist(t(rProK), "jaccard"),method="average"),weights),main="16S")
dev.off()

##Lets use dbRDA to examine the effect of env variables
model.full <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean + temp, data = envmatrix)
plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.ProK.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
sink()


model.3.1 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean + temp, data = envmatrix)
model.3.2 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean + temp, data = envmatrix)
model.3.3 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + impact + temp, data = envmatrix)
model.3.4 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean, data = envmatrix)

model.2.1 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + impact, data = envmatrix)
model.2.2 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean, data = envmatrix)
model.2.3 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ SSS + temp, data = envmatrix)
model.2.4 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean, data = envmatrix)
model.2.5 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ impact + temp, data = envmatrix)
model.2.6 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ chla_3yrmean + temp, data = envmatrix)

model.1.1 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~temp, data = envmatrix)
model.1.2 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~impact, data = envmatrix)
model.1.3 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~chla_3yrmean, data = envmatrix)
model.1.4 <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~SSS, data = envmatrix)

R2.1.1 <- unname(RsquareAdj(model.1.1)[[2]])*100
names(R2.1.1) <- "SST"
R2.1.2 <- unname(RsquareAdj(model.1.2)[[2]])*100
names(R2.1.2) <- "Impact"
R2.1.3 <- unname(RsquareAdj(model.1.3)[[2]])*100
names(R2.1.3) <- "Chl_a"
R2.1.4 <- unname(RsquareAdj(model.1.4)[[2]])*100
names(R2.1.4) <- "SSS"

R2.2.1 <- unname(RsquareAdj(model.2.1)[[2]])*100 
names(R2.2.1) <- "SSS&Impact"
R2.2.2 <- unname(RsquareAdj(model.2.2)[[2]])*100 
names(R2.2.2) <- "SSS&Chl_a"
R2.2.3 <- unname(RsquareAdj(model.2.3)[[2]])*100 
names(R2.2.3) <- "SSS&SST"
R2.2.4 <- unname(RsquareAdj(model.2.4)[[2]])*100 
names(R2.2.4) <- "Impact&Chl_a"
R2.2.5 <- unname(RsquareAdj(model.2.5)[[2]])*100 
names(R2.2.5) <- "Impact&SST"
R2.2.6 <- unname(RsquareAdj(model.2.6)[[2]])*100 
names(R2.2.6) <- "Chl_a&SST"

R2.3.1 <- unname(RsquareAdj(model.3.1)[[2]])*100 
names(R2.3.1) <- "Impact&Chl_a&SST"
R2.3.2 <- unname(RsquareAdj(model.3.2)[[2]])*100 
names(R2.3.2) <- "SSS&Chl_a&SST"
R2.3.3 <- unname(RsquareAdj(model.3.3)[[2]])*100 
names(R2.3.3) <- "SSS&Impact&SST"
R2.3.4 <- unname(RsquareAdj(model.3.4)[[2]])*100 
names(R2.3.4) <- "SSS&Impact&Chl_a"

R2.full <- unname(RsquareAdj(model.full)[[2]])*100 
names(R2.full) <- "SSS&Impact&Chl_a&SST"

expressionInput <- c(R2.full,R2.3.2,R2.3.1,R2.3.3,R2.3.4,R2.2.6,R2.2.3,R2.2.1,R2.2.5,R2.2.2,R2.2.4,R2.1.1,R2.1.4,R2.1.3,R2.1.2)

R216S <- expressionInput


pdf("figures/upset.dbRDA.ProK.pdf",width = 9,height=5.5)
upset(fromExpression(expressionInput), 
      order.by = c("freq","degree"),
      decreasing=c("T","T"),
      point.size = 3.5, 
      line.size = 2,
      set_size.show =FALSE,
      show.numbers=FALSE,
      mainbar.y.label = "Variance Explained (%)",
      set_size.scale_max=1)
dev.off()

##


#Upset plots suck, lets try and make the bar ourselves
pdf("figures/Figure 2/bar.ProK.pdf",width = 6,height = 4)
par(mar=c(3,3,1,1))
barplot(expressionInput,names=F,col="black",space=rep(0.5,15),ylim=c(0,25),yaxt="n")
axis(2,at=seq(0,25,5),labels=seq(0,25,5))
dev.off()


#Mantel test
m.ProK.SST <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),SST.dist)
m.ProK.SSS <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),SSS.dist)
m.ProK.ChlA <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),ChlA.dist)
m.ProK.impact <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),Impact.dist)

##Here we run a partial mantel test
###but SST and distance are highly autocorrelated! Need to correct for this
p.ProK.SST <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),SST.dist,dist)
p.ProK.SSS <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.ProK.ChlA <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.ProK.impact <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.ProK.SST.c <- msr(m.ProK.SST,lw1, 10000)
m.ProK.SSS.c <- msr(m.ProK.SSS,lw1, 10000)
m.ProK.ChlA.c <- msr(m.ProK.ChlA,lw1, 10000)
m.ProK.impact.c <- msr(m.ProK.impact,lw1, 10000)


mantel.p.ProK.out <- c(p.ProK.SST$statistic,p.ProK.SST$signif,
                       p.ProK.SSS$statistic,p.ProK.SSS$signif,
                       p.ProK.ChlA$statistic,p.ProK.ChlA$signif,
                       p.ProK.impact$statistic,p.ProK.impact$signif)

mantel.c.ProK.out <- unname(c(m.ProK.SST.c$obs-m.ProK.SST.c$expvar["Expectation"],
                              m.ProK.SST.c$pvalue,
                              m.ProK.SSS.c$obs-m.ProK.SSS.c$expvar["Expectation"],
                              m.ProK.SSS.c$pvalue,
                              m.ProK.ChlA.c$obs-m.ProK.ChlA.c$expvar["Expectation"],
                              m.ProK.ChlA.c$pvalue,
                              m.ProK.impact.c$obs-m.ProK.impact.c$expvar["Expectation"],
                              m.ProK.impact.c$pvalue))



###Lets pull toghether all the mantel tests

mantel.out <-rbind(mantel.p.COI.out,
                   mantel.c.COI.out,
                   mantel.p.18S.out,
                   mantel.c.18S.out,
                   mantel.p.ProK.out,
                   mantel.c.ProK.out)
colnames(mantel.out) <- c("SST R","SST p Value",
                          "SSS R","SSS p Value",
                          "ChlA R","ChlA p Value",
                          "Impact R","Impact p Value")

write.csv(round(mantel.out,3),"model.output/mantel.csv")

#Gather all the dbRDA R2 values

allR2tab <- rbind(R2COI,R218S,R216S)





####====5.0 Distance-Decay====####
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
pdf("figures/Figure3/log10Type.COI.dist.jacc.pdf",width = 6,height=5)
par(mar=c(4,4,1,1))
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type))
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


dist.decay.all <- cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("COI",length(pairwise.nomixed$value)))


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
pdf("figures/Figure3/log10Type.18S.dist.jacc.pdf",width = 6,height=5)
par(mar=c(4,4,1,1))
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=17, col=as.factor(pairwise.nomixed$type))
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=17)

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



dist.decay.all <- rbind(dist.decay.all,cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("18S",length(pairwise.nomixed$value))))





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
pdf("figures/Figure3/log10Type.ProK.dist.jacc.pdf",width = 6,height=5)
par(mar=c(4,4,1,1))
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],log10(1-pairwise.nomixed$value),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=18, col=as.factor(pairwise.nomixed$type))
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=18)

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

dist.decay.all <- rbind(dist.decay.all,cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("16S",length(pairwise.nomixed$value))))

##Plot all the data toghther 

pdf("figures/Figure3/Dist.Decay.all.pdf",width = 6,height=5)
#palette(c('#BF3465','#50B29E','#731683'))
palette(c('grey80','grey50','grey20'))
par(mar=c(4,4,1,1))
plot(dist.decay.all[,1],dist.decay.all[,2],
     ylab="Compositional Similarity (1-Jaccard)",
     xlab="Distance Between Sites (Km)",
     col=as.numeric(as.factor(dist.decay.all[,3])),
     pch=rev(as.numeric(as.factor(dist.decay.all[,3]))+15),
     cex=0.8)
legend("topright",col=c('grey20','grey50','grey80'),legend = rev(levels(as.factor(dist.decay.all[,3]))),pch=c(16,17,18))
dev.off()




####====6.0 Phyla Effcts====####
## Now lets seperate out the phyla and check for the same patterns 


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





####====7.0 Dark Diversity====####

palette(c('#51B9E0','#4bad84','#E8A016'))
#First we make a subset with only the OTUs that have some taxonomy
assigns <- as.character(COIassignments$OTU[COIassignments$assignmentQual!="None"])
assigns <- assigns[!is.na(assigns)]
dCOI <- rCOI[na.omit(match(assigns,rownames(rCOI))),]

#now lets look at these patterns
pdf("figures/dark.COI.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
nMDS <- metaMDS(t(dCOI),"jaccard")
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordiellipse(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            kind = "sd",
            draw = "polygon",
            lty=0,
            col=1:3)
points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="m"],
       pch=16,cex=3)

points(nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",1],
       nMDS$points[latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n",2],
       col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)][latlongdata$type[match(rownames(nMDS$points),latlongdata$sitecode)]=="n"],
       pch=17,cex=3)
text(nMDS$points[,1],nMDS$points[,2],labels=rownames(nMDS$points))
text(-0.45,0.25,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

##PERMANOVA COI
PermDispCOI.d <- betadisper(vegdist(t(dCOI), "jaccard"),latlongdata$PERMori[match(colnames(dCOI),latlongdata$sitecode)])
anova(PermDispCOI.d)
#No sig difference in multivariance homogenity
adonis(vegdist(t(dCOI), "jaccard")~latlongdata$PERMori[match(colnames(dCOI),latlongdata$sitecode)])
pairwise.adonis(vegdist(t(dCOI), "jaccard"),latlongdata$PERMori[match(colnames(dCOI),latlongdata$sitecode)])


#What about distance decay?
palette(c('#D36526','#2671B4'))

data <- dCOI

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
plot(Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)],1-pairwise$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="COI")

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
summary(lm(log10(site.similarity)~distance*site.type))

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.similarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(site.similarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect

pdf("figures/DARK.log10Type.COI.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
#LogPLot
pdf("figures/DARK.log10Type.COI.dist.jacc.pdf",width = 6,height=5)
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


