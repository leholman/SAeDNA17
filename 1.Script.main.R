#############################################
####==== South African eDNA Analysis ====####
####==== Luke E. Holman====19.10.2020====####
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
library("betapart")
library("Hmisc")
library(UpSetR)
library("gdata")
library(adespatial)
library(spdep)
library(ade4)
library(prettymapr)
library(EcolUtils)
library("ape")


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
#By Marker
rCOI <- read.csv("cleaned/rarefied.COI.csv",row.names = 1)
r18S <- read.csv("cleaned/rarefied.18S.csv",row.names = 1)
rProK <- read.csv("cleaned/rarefied.ProK.csv",row.names = 1)

#By Marker then Domain
rCOI.euk <- read.csv("cleaned/rarefied.COI.euk.csv",row.names = 1)
r18S.euk <- read.csv("cleaned/rarefied.18S.euk.csv",row.names = 1)

#By Marker the Kingdom
rCOI.met <- read.csv("cleaned/rarefied.COI.met.csv",row.names=1)
rCOI.pts <- read.csv("cleaned/rarefied.COI.pts.csv",row.names=1)
r18S.met <- read.csv("cleaned/rarefied.18S.met.csv",row.names=1)
r18S.pts <- read.csv("cleaned/rarefied.18S.pts.csv",row.names=1)

#Read in taxa
#BLAST
COIassignments <- read.csv("Taxonomy/CleanedTaxonomy/COI.blast.csv")
z18Sassignments <- read.csv("Taxonomy/CleanedTaxonomy/18S.blast.csv")
#RDP
ProKtaxa <- read.csv("Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.csv")
ProKtaxaALL <- read.csv("Taxonomy/CleanedTaxonomy/ProK.dada2RDPsilva.ALL.csv")
COIassignments.RDP <- read.delim("Taxonomy/RDP.assigns/RDP.class.COI.v4.txt",sep="\t",header=F)
z18Sassignments.RDP <- read.delim("Taxonomy/RDP.assigns/RDP.class.18S.v3.2.txt",sep="\t",header=F)


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

#Hows does our measure of human impact correlate with how we describe the sites?

impact.score <- envdat$impact[match(colnames(rCOI),envdat$sitecode)]
site.type <- latlongdata$type[match(colnames(rCOI),latlongdata$sitecode)]

plot(as.numeric(as.factor(site.type)),impact.score,xlim=c(0.5,2.5),xaxt='n',xlab="",pch=16)
axis(1,c(1,2),labels=c("Artificial","Natural"))
n.mean <- aggregate(impact.score,list(site.type),mean)[1,2]
a.mean <- aggregate(impact.score,list(site.type),mean)[2,2]
lines(c(0.9,1.1),c(n.mean,n.mean),lwd=3)
lines(c(1.9,2.1),c(a.mean,a.mean),lwd=3)


#Normality tests, both not sig diff from normal (but low n) 
shapiro.test(impact.score[site.type=="n"])
shapiro.test(impact.score[site.type=="m"])

#Since we care about if artificial is higher than natural we run one tailed tests 
t.test(impact.score~site.type,alternative="greater")
wilcox.test(impact.score~site.type,alternative="greater")

#Data suggests some association between site type and impact score but no firm evidence. 

####====2.0 Alpha Diversity====####

## Patterns of alpha  taxonomic diversity across the coast 
#SppRichness by marker
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


#SppRichness by marker.domain
temp <- rCOI.euk
temp[temp > 1] <- 1
Coastal.alpha.domain <- data.frame("ID"=names(temp),"COI.euk"= colSums(temp))
temp <- r18S.euk
temp[temp > 1] <- 1
Coastal.alpha.domain$z18S.euk <- colSums(temp)
temp <- rProK
temp[temp > 1] <- 1
Coastal.alpha.domain$ProK <- colSums(temp)
#reorder for plotting
Coastal.alpha.domain <- Coastal.alpha.domain[na.omit(match(latlongdata$sitecode,Coastal.alpha.domain$ID)),]

#SppRichness by marker.kingdom
temp <- rCOI.met
temp[temp > 1] <- 1
Coastal.alpha.kingdom <- data.frame("ID"=names(temp),"COI.met"= colSums(temp))
temp <- r18S.met
temp[temp > 1] <- 1
Coastal.alpha.kingdom$z18S.met <- colSums(temp)
temp <- rCOI.pts
temp[temp > 1] <- 1
Coastal.alpha.kingdom$COI.pts <- colSums(temp)
temp <- r18S.pts
temp[temp > 1] <- 1
Coastal.alpha.kingdom$z18S.pts <- colSums(temp)
temp <- rProK
temp[temp > 1] <- 1
Coastal.alpha.kingdom$ProK <- colSums(temp)
#reorder for plotting
Coastal.alpha.kingdom <- Coastal.alpha.kingdom[na.omit(match(latlongdata$sitecode,Coastal.alpha.kingdom$ID)),]



#richness plot
pdf("figures/OTU.Richness.pdf",width = 9,height=5)
plot(1:18,Coastal.alpha$COI,xaxt='n',xlab="",ylab="OTUrichness",pch=16, col="Darkred",ylim=c(200,950))
points(1:18,Coastal.alpha$z18S,col="Darkblue",pch=16)
points(1:18,Coastal.alpha$ProK,col="Darkgreen",pch=16)
axis(1,at=1:18,labels=Coastal.alpha$ID,cex=0.7)
legend("topright",legend=c("COI","18S","16S"),col=c("Darkred","Darkblue","Darkgreen"),pch=16)
dev.off()

#richness by coast & marker
Coastal.alpha.tall <- melt(Coastal.alpha)
Coastal.alpha.tall$coast <- as.character(latlongdata$PERMori[match(as.character(Coastal.alpha.tall$ID),as.character(latlongdata$sitecode))])
Coastal.alpha.tall$comb <- factor(paste(Coastal.alpha.tall$variable,Coastal.alpha.tall$coast,sep="."),
                                     levels=unique(paste(Coastal.alpha.tall$variable,Coastal.alpha.tall$coast,sep="."))[c(1,4,7,2,5,8,3,6,9)])

#richness by coast & domain
Coastal.alpha.domain.tall <- melt(Coastal.alpha.domain)
Coastal.alpha.domain.tall$coast <- as.character(latlongdata$PERMori[match(as.character(Coastal.alpha.domain.tall$ID),as.character(latlongdata$sitecode))])
Coastal.alpha.domain.tall$comb <- factor(paste(Coastal.alpha.domain.tall$variable,Coastal.alpha.domain.tall$coast,sep="."),
                                  levels=unique(paste(Coastal.alpha.domain.tall$variable,Coastal.alpha.domain.tall$coast,sep="."))[c(1,4,7,2,5,8,3,6,9)])

#richness by coast & kingdom
Coastal.alpha.kingdom.tall <- melt(Coastal.alpha.kingdom)
Coastal.alpha.kingdom.tall$coast <- as.character(latlongdata$PERMori[match(as.character(Coastal.alpha.kingdom.tall$ID),as.character(latlongdata$sitecode))])
Coastal.alpha.kingdom.tall$comb <- factor(paste(Coastal.alpha.kingdom.tall$variable,Coastal.alpha.kingdom.tall$coast,sep="."),
                                         levels=unique(paste(Coastal.alpha.kingdom.tall$variable,Coastal.alpha.kingdom.tall$coast,sep="."))[c(1,4,7,10,13,2,5,8,11,14,3,6,9,12,15)])




#Test for sig diff by coast across markers
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


#Test for sig diff by coast across domains
##ANOVA assump - Are the residuals normally distributed? 

resCOI.euk <- residuals(lm(value~coast,data=Coastal.alpha.domain.tall[Coastal.alpha.domain.tall$variable=="COI.euk",]))
shapiro.test(resCOI.euk)
#yes 
res18S.euk <- residuals(lm(value~coast,data=Coastal.alpha.domain.tall[Coastal.alpha.domain.tall$variable=="z18S.euk",]))
shapiro.test(res18S.euk)
#yes
resProK <- residuals(lm(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="ProK",]))
shapiro.test(resProK)
#yes

## Do coasts have equal variance?
bartlett.test(value~coast,data=Coastal.alpha.domain.tall[Coastal.alpha.domain.tall$variable=="COI.euk",])
#yes
bartlett.test(value~coast,data=Coastal.alpha.domain.tall[Coastal.alpha.domain.tall$variable=="z18S.euk",])
#yes
bartlett.test(value~coast,data=Coastal.alpha.tall[Coastal.alpha.tall$variable=="ProK",])
#yes

ANOVA.COI.euk <- aov(value~coast,data=Coastal.alpha.domain.tall[Coastal.alpha.domain.tall$variable=="COI.euk",])
ANOVA.18S.euk <- aov(value~coast,data=Coastal.alpha.domain.tall[Coastal.alpha.domain.tall$variable=="z18S.euk",])

summary(ANOVA.COI.euk) #Global P =0.095
summary(ANOVA.18S.euk) #Global P =0.174

TukeyHSD(ANOVA.COI.euk)
TukeyHSD(ANOVA.18S.euk)

#Now by marker kingdom
##ANOVA assump - Are the residuals normally distributed? 

resCOI.met <- residuals(lm(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="COI.met",]))
shapiro.test(resCOI.met)
#yes 
resCOI.pts <- residuals(lm(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="COI.pts",]))
shapiro.test(resCOI.pts)
#yes
res18S.met <- residuals(lm(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="z18S.met",]))
shapiro.test(res18S.met)
#yes
res18S.pts <- residuals(lm(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="z18S.pts",]))
shapiro.test(res18S.pts)


## Do coasts have equal variance?
bartlett.test(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="COI.met",])
bartlett.test(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="COI.pts",])
bartlett.test(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="z18S.met",])
bartlett.test(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="z18S.pts",])

ANOVA.COI.met <- aov(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="COI.met",])
ANOVA.COI.pts <- aov(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="COI.pts",])
ANOVA.18S.met <- aov(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="z18S.met",])
ANOVA.18S.pts <- aov(value~coast,data=Coastal.alpha.kingdom.tall[Coastal.alpha.kingdom.tall$variable=="z18S.pts",])

summary(ANOVA.COI.met) #Global P =0.178
summary(ANOVA.COI.pts) #Global P =0.024 *
summary(ANOVA.18S.met) #Global P =0.122
summary(ANOVA.18S.pts) #Global P =0.273


TukeyHSD(ANOVA.COI.pts) 




#Plots
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


pdf("figures/Figure1/OTU.coast.Richness.Domains.pdf",width = 9,height=3)
par(mar=c(2,4,1,1))
palette(c('#54B4E9','#F0E442','#CC79A7'))
plot(jitter(as.numeric(Coastal.alpha.domain.tall$comb),0.8),
     Coastal.alpha.domain.tall$value,
     type="n",
     xaxt="n",
     ylab="ASVs")
axis(1,at=1:9,labels=rep(c("Euk(COI)","Euk(18S)","ProK(16S)"),3),cex=0.8)

#Add ecoregions
rect(0,0,3.5,1000,border = NA,col="#D3ECF6")
rect(3.5,0,6.5,1000,border = NA,col="#B6DDD3")
rect(6.5,0,10,1000,border = NA,col="#F9E2B4")

#Add points
points(jitter(as.numeric(Coastal.alpha.domain.tall$comb),0.8),
       Coastal.alpha.domain.tall$value,pch=16)

##Lets draw in those mean values
for (num in 1:9){
  run.mean <- mean(Coastal.alpha.domain.tall$value[as.numeric(Coastal.alpha.domain.tall$comb)==num])
  lines(c(num-0.2,num+0.2),c(run.mean,run.mean),lwd=2)
}
box()
dev.off()


#test kingdom plot

Coastal.alpha.kingdom.tall.2 <- Coastal.alpha.kingdom.tall[!(Coastal.alpha.kingdom.tall$variable=="COI.pts"| Coastal.alpha.kingdom.tall$variable=="z18S.met"),]
Coastal.alpha.kingdom.tall.2$variable <- droplevels(Coastal.alpha.kingdom.tall.2$variable)



pdf("figures/test.kingdom.richness.pdf",width = 8,height=5)
par(mfrow=c(1,3))
for (dataset in unique(Coastal.alpha.kingdom.tall.2$variable)){
  boxplot(Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset]~Coastal.alpha.kingdom.tall.2$comb[Coastal.alpha.kingdom.tall.2$variable==dataset],
          main=dataset,drop=TRUE,
          xlab="Factor",
          ylab="OTUs",
          col=1:3
  )
  points(c(rep(1,6),rep(2,6),rep(3,6)),Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset],pch=16,cex=1.5)
} 
dev.off()

pdf("figures/Figure1/test.kingdom.richness.2.pdf",width = 9,height=3)
par(mar=c(2,4,1,1),mfrow=c(1,3))
palette(c('#51B9E0','#4bad84','#E8A016'))
for (dataset in unique(Coastal.alpha.kingdom.tall.2$variable) ){
  if (dataset=="COI.met"){assign("ylabel","")}else{assign("ylabel","")}
  plot(jitter(c(rep(1,6),rep(2,6),rep(3,6)),0.8),
       Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset],
       #ylim=c(min(Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset]),
        #      max(Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset])+
        #        20),
       type="n",
       xaxt="n",
       yaxt="n",
       ylab=ylabel,
       bty="L",
       xlim=c(0.5,3.5))
  axis(2,las=1,cex.axis=1.6)
  points(jitter(c(rep(1,6),rep(2,6),rep(3,6)),0.8),
         Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset],
         col=c(rep(1,6),rep(2,6),rep(3,6)),pch=16,cex=2)
  axis(1,1:3,labels=c("West","South","East"),lwd=0,lwd.ticks=1,cex.axis=1.8)
  run.mean <-unname(unlist(aggregate(Coastal.alpha.kingdom.tall.2$value[Coastal.alpha.kingdom.tall.2$variable==dataset]~as.factor(c(rep(1,6),rep(2,6),rep(3,6))), FUN=mean))[4:6])
for (num in 1:3){
  lines(c(num-0.2,num+0.2),c(run.mean[num],run.mean[num]),lwd=3)
}
}
dev.off()

#Add ecoregions
rect(0,0,3.5,1000,border = NA,col="#D3ECF6")
rect(3.5,0,6.5,1000,border = NA,col="#B6DDD3")
rect(6.5,0,10,1000,border = NA,col="#F9E2B4")

#Add points
points(jitter(as.numeric(Coastal.alpha.domain.tall$comb),0.8),
       Coastal.alpha.domain.tall$value,pch=16)

##Lets draw in those mean values
for (num in 1:9){
  run.mean <- mean(Coastal.alpha.domain.tall$value[as.numeric(Coastal.alpha.domain.tall$comb)==num])
  lines(c(num-0.2,num+0.2),c(run.mean,run.mean),lwd=2)
}
box()







####====3.0 Taxonomy====####
#Taxonomic profile
#set NA as "" in taxonomy file

##We first make a function that calculates the count table for us in two ways
CountTable <- function(in.taxonomy,in.data,output="Count",some.unassigned=T){
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
  if(some.unassigned==T){rownames(out.dat)[1] <- "Unassigned"}
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

#Lets make them for the data subsets by kingdom
getPalette = colorRampPalette(brewer.pal(9, "Pastel1"))

COI.met.taxa.abun <- CountTable(COIassignments.RDP$V12[match(rownames(rCOI.met),COIassignments.RDP$V1)],rCOI.met,output="Abundance",some.unassigned=F)
z18S.pts.taxa.abun <- CountTable(z18Sassignments.RDP$V12[match(rownames(r18S.pts),z18Sassignments.RDP$V1)],r18S.pts,output="Abundance",some.unassigned=F)
ProK.taxa.abun.2 <- CountTable(as.character(ProKtaxa$Phylum),rProK,output="Abundance",some.unassigned=F)[-1,]

COI.met.taxa.count <- CountTable(COIassignments.RDP$V12[match(rownames(rCOI.met),COIassignments.RDP$V1)],rCOI.met,output="Count",some.unassigned=F)
z18S.pts.taxa.count <- CountTable(z18Sassignments.RDP$V12[match(rownames(r18S.pts),z18Sassignments.RDP$V1)],r18S.pts,output="Count",some.unassigned=F)
ProK.taxa.count.2 <- CountTable(as.character(ProKtaxa$Phylum),rProK,output="Count",some.unassigned=F)[-1,]

COI.met.taxaprop <- prop.table(as.matrix(minAbundance(COI.met.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(COI.met.taxa.count)))],minAbun= 0.02)),2)
z18S.pts.taxaprop <- prop.table(as.matrix(minAbundance(z18S.pts.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(z18S.pts.taxa.count)))],minAbun= 0.02)),2)
ProK.taxaprop <- prop.table(as.matrix(minAbundance(ProK.taxa.count.2[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.count.2)))],minAbun= 0.02)),2)

COI.met.taxaprop.abun <- prop.table(as.matrix(minAbundance(COI.met.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(COI.met.taxa.count)))],minAbun= 0.01)),2)
z18S.pts.taxaprop.abun <- prop.table(as.matrix(minAbundance(z18S.pts.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(z18S.pts.taxa.count)))],minAbun= 0.01)),2)
ProK.taxaprop.abun <- prop.table(as.matrix(minAbundance(ProK.taxa.abun.2[,na.omit(match(latlongdata$sitecode,colnames(ProK.taxa.count.2)))],minAbun= 0.01)),2)




#The protists are a mess - let's change them to supergroups

phyla2supergroup <- read.csv("Taxonomy/ProtistDesignations/Phyla2SuperGroups.csv")

temp <- as.data.frame(matrix(0,ncol=18,nrow =length(unique(phyla2supergroup$SuperGroup))))
colnames(temp) <- colnames(z18S.pts.taxa.count)
rownames(temp) <- sort(unique(phyla2supergroup$SuperGroup))
for (group in sort(unique(phyla2supergroup$SuperGroup))){
  temp[group,] <- colSums(z18S.pts.taxa.count[phyla2supergroup$SuperGroup==group,])
}
z18S.pts.taxa.count <- temp[rowSums(temp)>1,]
z18S.pts.taxaprop <- prop.table(as.matrix(z18S.pts.taxa.count[,na.omit(match(latlongdata$sitecode,colnames(z18S.pts.taxa.count)))]),2)

temp <- as.data.frame(matrix(0,ncol=18,nrow =length(unique(phyla2supergroup$SuperGroup))))
colnames(temp) <- colnames(z18S.pts.taxa.count)
rownames(temp) <- sort(unique(phyla2supergroup$SuperGroup))
for (group in sort(unique(phyla2supergroup$SuperGroup))){
  temp[group,] <- colSums(z18S.pts.taxa.abun[phyla2supergroup$SuperGroup==group,])
}
z18S.pts.taxa.abun <- temp[rowSums(temp)>1,]
z18S.pts.taxaprop.abun <- prop.table(as.matrix(z18S.pts.taxa.abun[,na.omit(match(latlongdata$sitecode,colnames(z18S.pts.taxa.abun)))]),2)


pdf("figures/Figure1/COI.met.abun.pdf",width = 11, height = 2.3)
par(mfrow=c(1,1),mai = c(0.1, 0.6, 0.1,0),xpd=TRUE)
barplot(COI.met.taxaprop[dim(COI.met.taxaprop)[1]:1,],col=rev(getPalette(dim(COI.met.taxaprop)[1])),axisnames=FALSE,cex.axis=1.35,las=1,
        #you can set the distance between the y axis and the bars by using the two below options, the xlim changes the distances
        xaxs = "i",xlim = c(-0.3, 22))
dev.off()

pdf("figures/Figure1/18S.pts.abun.pdf",width = 11, height = 2.3)
par(mfrow=c(1,1),mai = c(0.1, 0.6, 0.1,0),xpd=TRUE)
barplot(z18S.pts.taxaprop[dim(z18S.pts.taxaprop)[1]:1,],col=rev(getPalette(dim(z18S.pts.taxaprop)[1])),axisnames=FALSE,cex.axis=1.35,las=1,
        #you can set the distance between the y axis and the bars by using the two below options, the xlim changes the distances
        xaxs = "i",xlim = c(-0.3, 22))
dev.off()

pdf("figures/Figure1/ProK.abun.pdf",width = 11, height = 2.3)
par(mfrow=c(1,1),mai = c(0.1, 0.6, 0.1,0),xpd=TRUE)
barplot(ProK.taxaprop[dim(ProK.taxaprop)[1]:1,],col=rev(getPalette(dim(ProK.taxaprop)[1])),axisnames=FALSE,cex.axis=1.35,las=1,
        #you can set the distance between the y axis and the bars by using the two below options, the xlim changes the distances
        xaxs = "i",xlim = c(-0.3, 22))
dev.off()



pdf("figures/Figure1/COI.met.legend.pdf",width = 2.8, height = 3.4)
par(mai = c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
legend(2,10,rownames(COI.met.taxaprop),fill=getPalette(dim(COI.met.taxaprop)[1]),cex=1.2,bty = "n",y.intersp=0.75)
dev.off()

pdf("figures/Figure1/18S.pts.legend.pdf",width = 2.8, height = 3.4)
par(mai = c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
legend(2,10,rownames(z18S.pts.taxaprop),fill=getPalette(dim(z18S.pts.taxaprop)[1]),cex=1.2,bty = "n",y.intersp=0.75)
dev.off()

pdf("figures/Figure1/ProK.legend.pdf",width = 2.8, height = 3.4)
par(mai = c(0,0,0,0))
plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
legend(2,10,rownames(ProK.taxaprop),fill=getPalette(dim(ProK.taxaprop)[1]),cex=1.2,bty = "n",y.intersp=0.75)
dev.off()

#Lets make a figure for the supplement 
pdf("figures/AbundanceTab.RDP.Suppl.pdf",width = 11, height = 7)
par(mfrow=c(3,1),mai = c(0.3, 0.6, 0.1,2.2),xpd=TRUE)
barplot(COI.met.taxaprop.abun[dim(COI.met.taxaprop.abun)[1]:1,],col=rev(getPalette(dim(COI.met.taxaprop.abun)[1])),axisnames=FALSE,cex.axis=1.35,las=1,xaxs = "i",xlim = c(-0.3, 22))
legend(22,1,rownames(COI.met.taxaprop.abun),fill=getPalette(dim(COI.met.taxaprop.abun)[1]),cex=1.2,bty = "n",y.intersp=0.75)
barplot(z18S.pts.taxaprop.abun[dim(z18S.pts.taxaprop.abun)[1]:1,],col=rev(getPalette(dim(z18S.pts.taxaprop.abun)[1])),axisnames=FALSE,cex.axis=1.35,las=1,xaxs = "i",xlim = c(-0.3, 22))
legend(22,0.8,rownames(z18S.pts.taxaprop.abun),fill=getPalette(dim(z18S.pts.taxaprop.abun)[1]),cex=1.2,bty = "n",y.intersp=0.75)
barplot(ProK.taxaprop.abun[dim(ProK.taxaprop.abun)[1]:1,],col=rev(getPalette(dim(ProK.taxaprop.abun)[1])),cex.axis=1.35,las=1, xaxs = "i",xlim = c(-0.3, 22),cex.names=1.4)
legend(21.5,1,rownames(ProK.taxaprop.abun),fill=getPalette(dim(ProK.taxaprop.abun)[1]),cex=1.2,bty = "n",y.intersp=0.75,xpd=NA)
dev.off()



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
####======4.1 Kingdom Analyses====####


#Lets look at the effects by marker kingdom
#COI metazoans first
nMDS <- metaMDS(vegdist(t(rCOI.met),"jaccard",binary=TRUE))

#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/COI.met.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            lty=1,
            lwd=2,
            col=1:3)
ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
text(-0.5,0.3,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

sink("model.output/beta.COI.met.txt")
##PERMANOVA COI.met
PermDispCOI.met <- betadisper(vegdist(t(rCOI.met), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rCOI.met),latlongdata$sitecode)])
anova(PermDispCOI.met)
#No sig difference in multivariate homogenity
adonis(vegdist(t(rCOI.met), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(rCOI.met),latlongdata$sitecode)])
adonis.pair(vegdist(t(rCOI.met), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rCOI.met),latlongdata$sitecode)])
#p<0.001
sink()
#Testing mantel and partial mantel tests

#Lets create a geographic distance matrix for all sites  
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


##Here we run a partial mantel test
#
p.COI.met.SST <- mantel.partial(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),SST.dist,dist)
p.COI.met.SSS <- mantel.partial(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.COI.met.ChlA <- mantel.partial(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.COI.met.impact <- mantel.partial(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141
#Mantel test
m.COI.met.SST <- mantel.randtest(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),SST.dist)
m.COI.met.SSS <- mantel.randtest(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),SSS.dist)
m.COI.met.ChlA <- mantel.randtest(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),ChlA.dist)
m.COI.met.impact <- mantel.randtest(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE),Impact.dist)

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.COI.met.SST.c <- msr(m.COI.met.SST,lw1, 10000)
m.COI.met.SSS.c <- msr(m.COI.met.SSS,lw1, 10000)
m.COI.met.ChlA.c <- msr(m.COI.met.ChlA,lw1, 10000)
m.COI.met.impact.c <- msr(m.COI.met.impact,lw1, 10000)


mantel.p.COI.met.out <- c(p.COI.met.SST$statistic,p.COI.met.SST$signif,
                      p.COI.met.SSS$statistic,p.COI.met.SSS$signif,
                      p.COI.met.ChlA$statistic,p.COI.met.ChlA$signif,
                      p.COI.met.impact$statistic,p.COI.met.impact$signif)

mantel.c.COI.met.out <- unname(c(m.COI.met.SST.c$obs-m.COI.met.SST.c$expvar["Expectation"],
                             m.COI.met.SST.c$pvalue,
                             m.COI.met.SSS.c$obs-m.COI.met.SSS.c$expvar["Expectation"],
                             m.COI.met.SSS.c$pvalue,
                             m.COI.met.ChlA.c$obs-m.COI.met.ChlA.c$expvar["Expectation"],
                             m.COI.met.ChlA.c$pvalue,
                             m.COI.met.impact.c$obs-m.COI.met.impact.c$expvar["Expectation"],
                             m.COI.met.impact.c$pvalue))



#We see no effect of salinity or ChlA after geography is taken into account, 

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(rCOI),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(rCOI),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(rCOI),envdat$sitecode)])

#now let's plot


#Ordisurf - temp
pdf("figures/COI.met.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.met.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - impact
pdf("figures/COI.met.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.met.impact.txt")
summary(surfout.ipt)
sink()


##Now let's use a redundancy analysis to look at the relative role of the env factors that remain important

model.full <- dbrda(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE)~ envmatrix$temp + envmatrix$impact)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.COI.met.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
varpart(vegdist(t(rCOI.met),method = "jaccard",binary = TRUE), envmatrix$temp,envmatrix$impact)
sink()


#COI protists
nMDS <- metaMDS(vegdist(t(rCOI.pts),"jaccard",binary=TRUE))

#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/COI.pts.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
           lty=1,
           lwd=2,
           col=1:3)
ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
text(-0.45,0.2,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

##PERMANOVA COI.pts
sink("model.output/beta.COI.pts.txt")
PermDispCOI.pts <- betadisper(vegdist(t(rCOI.pts), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rCOI.pts),latlongdata$sitecode)])
anova(PermDispCOI.pts) #p = 0.06 difference driven by east to west 
TukeyHSD(PermDispCOI.pts)
#No sig difference in multivariance homogenity
adonis(vegdist(t(rCOI.pts), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(rCOI.pts),latlongdata$sitecode)])
adonis.pair(vegdist(t(rCOI.pts), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rCOI.pts),latlongdata$sitecode)])
#p<0.01
sink()


#Testing mantel and partial mantel tests

#Lets create a geographic distance matrix for all sites  
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


##Here we run a partial mantel test
###
p.COI.pts.SST <- mantel.partial(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),SST.dist,dist)
p.COI.pts.SSS <- mantel.partial(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.COI.pts.ChlA <- mantel.partial(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.COI.pts.impact <- mantel.partial(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141
#Mantel test
m.COI.pts.SST <- mantel.randtest(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),SST.dist)
m.COI.pts.SSS <- mantel.randtest(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),SSS.dist)
m.COI.pts.ChlA <- mantel.randtest(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),ChlA.dist)
m.COI.pts.impact <- mantel.randtest(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE),Impact.dist)

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.COI.pts.SST.c <- msr(m.COI.pts.SST,lw1, 10000)
m.COI.pts.SSS.c <- msr(m.COI.pts.SSS,lw1, 10000)
m.COI.pts.ChlA.c <- msr(m.COI.pts.ChlA,lw1, 10000)
m.COI.pts.impact.c <- msr(m.COI.pts.impact,lw1, 10000)


mantel.p.COI.pts.out <- c(p.COI.pts.SST$statistic,p.COI.pts.SST$signif,
                          p.COI.pts.SSS$statistic,p.COI.pts.SSS$signif,
                          p.COI.pts.ChlA$statistic,p.COI.pts.ChlA$signif,
                          p.COI.pts.impact$statistic,p.COI.pts.impact$signif)

mantel.c.COI.pts.out <- unname(c(m.COI.pts.SST.c$obs-m.COI.pts.SST.c$expvar["Expectation"],
                                 m.COI.pts.SST.c$pvalue,
                                 m.COI.pts.SSS.c$obs-m.COI.pts.SSS.c$expvar["Expectation"],
                                 m.COI.pts.SSS.c$pvalue,
                                 m.COI.pts.ChlA.c$obs-m.COI.pts.ChlA.c$expvar["Expectation"],
                                 m.COI.pts.ChlA.c$pvalue,
                                 m.COI.pts.impact.c$obs-m.COI.pts.impact.c$expvar["Expectation"],
                                 m.COI.pts.impact.c$pvalue))


#We see an effect of temp and impact

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(rCOI),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(rCOI),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(rCOI),envdat$sitecode)])

#now let's plot


#Ordisurf - impact
pdf("figures/COI.pts.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.pts.impact.txt")
summary(surfout.ipt)
sink()

#Ordisurf - temp
pdf("figures/COI.pts.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.pts.temp.txt")
summary(surfout.ipt)
sink()


##Now let's use a redundancy analysis to look at the relative role of the env factors that remain important

model.full <- dbrda(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE)~ envmatrix$temp + envmatrix$impact)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.COI.pts.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
#Var partitioning
varpart(vegdist(t(rCOI.pts),method = "jaccard",binary = TRUE), envmatrix$temp,envmatrix$impact)
sink()


#18S metazoans
nMDS <- metaMDS(vegdist(t(r18S.met),"jaccard",binary=TRUE))

#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/18S.met.beta.txt",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
           lty=1,
           lwd=2,
           col=1:3)
ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
text(-0.5,0.25,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

##PERMANOVA 18S.met
sink("model.output/beta.18S.met.txt")
PermDisp18S.met <- betadisper(vegdist(t(r18S.met), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(r18S.met),latlongdata$sitecode)])
anova(PermDisp18S.met)
#No sig difference in multivariance homogenity
adonis(vegdist(t(r18S.met), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(r18S.met),latlongdata$sitecode)])
adonis.pair(vegdist(t(r18S.met), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(r18S.met),latlongdata$sitecode)])
#p<0.01
sink()

#Testing mantel and partial mantel tests

#Lets create a geographic distance matrix for all sites  
dist <- read.csv("metadata/distance.csv",row.names = 1)
dist <- as.matrix(dist)
upperTriangle(dist) <- lowerTriangle(dist,byrow = TRUE)
dist <- dist[match(colnames(r18S),colnames(dist)),]
dist <- dist[,match(colnames(r18S),colnames(dist))]
dist <- as.dist(dist)

#Lets create a false 2d world for calculations
dist.2d <- read.csv("metadata/distance.csv",row.names = 1)
dist.2d <- as.matrix(dist.2d)[,1]
dist.2d <- dist.2d[match(colnames(r18S),names(dist.2d))]
dist.2d <- data.frame("Y"=rep(1,18),"X"=dist.2d)

#Lets create distance matrices for each env parameter
SST.dist <- dist(tempdat$TempAvr[match(colnames(r18S),tempdat$Site.Code)])
SSS.dist <- dist(envdat$SSS[match(colnames(r18S),envdat$sitecode)])
ChlA.dist <- dist(envdat$chla_3yrmean[match(colnames(r18S),envdat$sitecode)])
Impact.dist <- dist(envdat$impact[match(colnames(r18S),envdat$sitecode)])


##Here we run a partial mantel test
#
p.18S.met.SST <- mantel.partial(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),SST.dist,dist)
p.18S.met.SSS <- mantel.partial(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.18S.met.ChlA <- mantel.partial(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.18S.met.impact <- mantel.partial(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141
#Mantel test
m.18S.met.SST <- mantel.randtest(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),SST.dist)
m.18S.met.SSS <- mantel.randtest(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),SSS.dist)
m.18S.met.ChlA <- mantel.randtest(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),ChlA.dist)
m.18S.met.impact <- mantel.randtest(vegdist(t(r18S.met),method = "jaccard",binary = TRUE),Impact.dist)

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.18S.met.SST.c <- msr(m.18S.met.SST,lw1, 10000)
m.18S.met.SSS.c <- msr(m.18S.met.SSS,lw1, 10000)
m.18S.met.ChlA.c <- msr(m.18S.met.ChlA,lw1, 10000)
m.18S.met.impact.c <- msr(m.18S.met.impact,lw1, 10000)


mantel.p.18S.met.out <- c(p.18S.met.SST$statistic,p.18S.met.SST$signif,
                          p.18S.met.SSS$statistic,p.18S.met.SSS$signif,
                          p.18S.met.ChlA$statistic,p.18S.met.ChlA$signif,
                          p.18S.met.impact$statistic,p.18S.met.impact$signif)

mantel.c.18S.met.out <- unname(c(m.18S.met.SST.c$obs-m.18S.met.SST.c$expvar["Expectation"],
                                 m.18S.met.SST.c$pvalue,
                                 m.18S.met.SSS.c$obs-m.18S.met.SSS.c$expvar["Expectation"],
                                 m.18S.met.SSS.c$pvalue,
                                 m.18S.met.ChlA.c$obs-m.18S.met.ChlA.c$expvar["Expectation"],
                                 m.18S.met.ChlA.c$pvalue,
                                 m.18S.met.impact.c$obs-m.18S.met.impact.c$expvar["Expectation"],
                                 m.18S.met.impact.c$pvalue))


#We see no effect of salinity or ChlA after geography is taken into account, 

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(r18S),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(r18S),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(r18S),envdat$sitecode)])

#now let's plot


#Ordisurf - temp
pdf("figures/18S.met.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.met.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - impact
pdf("figures/18S.met.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.met.impact.txt")
summary(surfout.ipt)
sink()


##Now let's use a redundancy analysis to look at the relative role of the env factors that remain important

model.full <- dbrda(vegdist(t(r18S.met),method = "jaccard",binary = TRUE)~ envmatrix$temp + envmatrix$impact)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.18S.met.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
varpart(vegdist(t(r18S.met),method = "jaccard",binary = TRUE), envmatrix$temp,envmatrix$impact)
sink()


#18S protists
nMDS <- metaMDS(vegdist(t(r18S.pts),"jaccard",binary=TRUE))

#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/18S.pts.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
           lty=1,
           lwd=2,
           col=1:3)
ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
text(-0.32,0.43,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

##PERMANOVA 18S.pts
sink("model.output/beta.18S.pts.txt")
PermDisp18S.pts <- betadisper(vegdist(t(r18S.pts), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(r18S.pts),latlongdata$sitecode)])
anova(PermDisp18S.pts)
#No sig difference in multivariance homogenity
adonis(vegdist(t(r18S.pts), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(r18S.pts),latlongdata$sitecode)])
adonis.pair(vegdist(t(r18S.pts), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(r18S.pts),latlongdata$sitecode)])
#p<0.01
sink()

#Testing mantel and partial mantel tests

#Lets create a geographic distance matrix for all sites  
dist <- read.csv("metadata/distance.csv",row.names = 1)
dist <- as.matrix(dist)
upperTriangle(dist) <- lowerTriangle(dist,byrow = TRUE)
dist <- dist[match(colnames(r18S),colnames(dist)),]
dist <- dist[,match(colnames(r18S),colnames(dist))]
dist <- as.dist(dist)

#Lets create a false 2d world for calculations
dist.2d <- read.csv("metadata/distance.csv",row.names = 1)
dist.2d <- as.matrix(dist.2d)[,1]
dist.2d <- dist.2d[match(colnames(r18S),names(dist.2d))]
dist.2d <- data.frame("Y"=rep(1,18),"X"=dist.2d)

#Lets create distance matrices for each env parameter
SST.dist <- dist(tempdat$TempAvr[match(colnames(r18S),tempdat$Site.Code)])
SSS.dist <- dist(envdat$SSS[match(colnames(r18S),envdat$sitecode)])
ChlA.dist <- dist(envdat$chla_3yrmean[match(colnames(r18S),envdat$sitecode)])
Impact.dist <- dist(envdat$impact[match(colnames(r18S),envdat$sitecode)])


##Here we run a partial mantel test
###
p.18S.pts.SST <- mantel.partial(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),SST.dist,dist)
p.18S.pts.SSS <- mantel.partial(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.18S.pts.ChlA <- mantel.partial(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.18S.pts.impact <- mantel.partial(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141
#Mantel test
m.18S.pts.SST <- mantel.randtest(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),SST.dist)
m.18S.pts.SSS <- mantel.randtest(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),SSS.dist)
m.18S.pts.ChlA <- mantel.randtest(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),ChlA.dist)
m.18S.pts.impact <- mantel.randtest(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE),Impact.dist)

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.18S.pts.SST.c <- msr(m.18S.pts.SST,lw1, 10000)
m.18S.pts.SSS.c <- msr(m.18S.pts.SSS,lw1, 10000)
m.18S.pts.ChlA.c <- msr(m.18S.pts.ChlA,lw1, 10000)
m.18S.pts.impact.c <- msr(m.18S.pts.impact,lw1, 10000)


mantel.p.18S.pts.out <- c(p.18S.pts.SST$statistic,p.18S.pts.SST$signif,
                          p.18S.pts.SSS$statistic,p.18S.pts.SSS$signif,
                          p.18S.pts.ChlA$statistic,p.18S.pts.ChlA$signif,
                          p.18S.pts.impact$statistic,p.18S.pts.impact$signif)

mantel.c.18S.pts.out <- unname(c(m.18S.pts.SST.c$obs-m.18S.pts.SST.c$expvar["Expectation"],
                                 m.18S.pts.SST.c$pvalue,
                                 m.18S.pts.SSS.c$obs-m.18S.pts.SSS.c$expvar["Expectation"],
                                 m.18S.pts.SSS.c$pvalue,
                                 m.18S.pts.ChlA.c$obs-m.18S.pts.ChlA.c$expvar["Expectation"],
                                 m.18S.pts.ChlA.c$pvalue,
                                 m.18S.pts.impact.c$obs-m.18S.pts.impact.c$expvar["Expectation"],
                                 m.18S.pts.impact.c$pvalue))


#We see no effect of temp or salinity  after geography is taken into account, 

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(r18S),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(r18S),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(r18S),envdat$sitecode)])

#now let's plot


#Ordisurf - ChlA
pdf("figures/18S.pts.chla.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$chla_3yrmean,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.pts.chla.txt")
summary(surfout.temp)
sink()

#Ordisurf - impact
pdf("figures/18S.pts.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.pts.impact.txt")
summary(surfout.ipt)
sink()


##Now let's use a redundancy analysis to look at the relative role of the env factors that remain important

model.full <- dbrda(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE)~ envmatrix$chla_3yrmean + envmatrix$impact)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.18S.pts.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
varpart(vegdist(t(r18S.pts),method = "jaccard",binary = TRUE), envmatrix$chla_3yrmean,envmatrix$impact)

sink()


##Here we analyse ProK

nMDS <- metaMDS(vegdist(t(rProK),"jaccard",binary=TRUE))

#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/16S.beta.pdf",width = 5,height=5)
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
           lty=1,
           lwd=2,
           col=1:3)
ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
text(-0.49,0.27,labels=paste0("Stress = ",round(nMDS$stress,3)))
dev.off()

##PERMANOVA 16S.pts
sink("model.output/beta.16s.txt")
PermDisp16S.pts <- betadisper(vegdist(t(rProK), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
anova(PermDisp16S.pts)
TukeyHSD(PermDisp16S.pts)
#Sig difference in multivariance homogenity
adonis(vegdist(t(rProK), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
adonis.pair(vegdist(t(rProK), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
#p<0.01
sink()

#Testing mantel and partial mantel tests

#Lets create a geographic distance matrix for all sites  
dist <- read.csv("metadata/distance.csv",row.names = 1)
dist <- as.matrix(dist)
upperTriangle(dist) <- lowerTriangle(dist,byrow = TRUE)
dist <- dist[match(colnames(rProK),colnames(dist)),]
dist <- dist[,match(colnames(rProK),colnames(dist))]
dist <- as.dist(dist)

#Lets create a false 2d world for calculations
dist.2d <- read.csv("metadata/distance.csv",row.names = 1)
dist.2d <- as.matrix(dist.2d)[,1]
dist.2d <- dist.2d[match(colnames(rProK),names(dist.2d))]
dist.2d <- data.frame("Y"=rep(1,18),"X"=dist.2d)

#Lets create distance matrices for each env parameter
SST.dist <- dist(tempdat$TempAvr[match(colnames(rProK),tempdat$Site.Code)])
SSS.dist <- dist(envdat$SSS[match(colnames(rProK),envdat$sitecode)])
ChlA.dist <- dist(envdat$chla_3yrmean[match(colnames(rProK),envdat$sitecode)])
Impact.dist <- dist(envdat$impact[match(colnames(rProK),envdat$sitecode)])


##Here we run a partial mantel test
###
p.16S.SST <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),SST.dist,dist)
p.16S.SSS <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.16S.ChlA <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.16S.impact <- mantel.partial(vegdist(t(rProK),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141
#Mantel test
m.16S.SST <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),SST.dist)
m.16S.SSS <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),SSS.dist)
m.16S.ChlA <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),ChlA.dist)
m.16S.impact <- mantel.randtest(vegdist(t(rProK),method = "jaccard",binary = TRUE),Impact.dist)

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.16S.SST.c <- msr(m.16S.SST,lw1, 10000)
m.16S.SSS.c <- msr(m.16S.SSS,lw1, 10000)
m.16S.ChlA.c <- msr(m.16S.ChlA,lw1, 10000)
m.16S.impact.c <- msr(m.16S.impact,lw1, 10000)


mantel.p.16S.out <- c(p.16S.SST$statistic,p.16S.SST$signif,
                          p.16S.SSS$statistic,p.16S.SSS$signif,
                          p.16S.ChlA$statistic,p.16S.ChlA$signif,
                          p.16S.impact$statistic,p.16S.impact$signif)

mantel.c.16S.out <- unname(c(m.16S.SST.c$obs-m.16S.SST.c$expvar["Expectation"],
                                 m.16S.SST.c$pvalue,
                                 m.16S.SSS.c$obs-m.16S.SSS.c$expvar["Expectation"],
                                 m.16S.SSS.c$pvalue,
                                 m.16S.ChlA.c$obs-m.16S.ChlA.c$expvar["Expectation"],
                                 m.16S.ChlA.c$pvalue,
                                 m.16S.impact.c$obs-m.16S.impact.c$expvar["Expectation"],
                                 m.16S.impact.c$pvalue))


#We see no effect of ChlA or salinity  after geography is taken into account, 

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(r18S),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(r18S),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(r18S),envdat$sitecode)])

#now let's plot


#Ordisurf - temp
pdf("figures/16S.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/16s.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - impact
pdf("figures/16S.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/16S.pts.impact.txt")
summary(surfout.ipt)
sink()


##Now let's use a redundancy analysis to look at the relative role of the env factors that remain important

model.full <- dbrda(vegdist(t(rProK),method = "jaccard",binary = TRUE)~ envmatrix$temp + envmatrix$impact)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.16S.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
varpart(vegdist(t(rProK),method = "jaccard",binary = TRUE), envmatrix$temp,envmatrix$impact)
sink()


##Now let's output the partial and correlated mantel tests
mantel.out <-rbind(mantel.p.COI.met.out,
                   mantel.c.COI.met.out,
                   mantel.p.18S.pts.out,
                   mantel.c.18S.pts.out,
                   mantel.p.16S.out,
                   mantel.c.16S.out)
colnames(mantel.out) <- c("SST R","SST p Value",
                          "SSS R","SSS p Value",
                          "ChlA R","ChlA p Value",
                          "Impact R","Impact p Value")

write.csv(round(mantel.out,3),"model.output/taxa.mantel.csv")




####======4.2 Marker Analyses====####
### First let's look at beta diversity by marker

##First COI
#PCoA
rCOI.dist <- vegdist(t(rCOI), "jaccard",binary = TRUE)
pcoa <- cmdscale(rCOI.dist)

#nMDS
nMDS <- metaMDS(t(rCOI),"jaccard",binary = TRUE)

##PERMANOVA COI
sink("model.output/beta.COI.txt")
PermDispCOI <- betadisper(vegdist(t(rCOI), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(rCOI),latlongdata$sitecode)])
anova(PermDispCOI)
#No sig difference in multivariance homogenity
adonis(vegdist(t(rCOI), "jaccard",binary = TRUE)~latlongdata$PERMori[match(colnames(rCOI),latlongdata$sitecode)])
adonis.pair(vegdist(t(rCOI), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(rCOI),latlongdata$sitecode)])
#p<0.001
sink()

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


#Testing mantel and partial mantel tests

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
r18S.dist <- vegdist(t(r18S), "jaccard",binary = TRUE)
pcoa <- cmdscale(r18S.dist)

#nMDS
nMDS <- metaMDS(t(r18S),"jaccard",binary = TRUE)

##PERMANOVA 18S
sink("model.output/beta.18S.txt")
PermDisp18S <- betadisper(vegdist(t(r18S), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(r18S),latlongdata$sitecode)])
anova(PermDisp18S)
#No sig difference in multivariance homogenity
adonis(vegdist(t(r18S), "jaccard",binary = TRUE)~latlongdata$PERMori[match(colnames(r18S),latlongdata$sitecode)])
adonis.pair(vegdist(t(r18S), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(r18S),latlongdata$sitecode)])
#p<0.001
sink()

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
ProK.dist <- vegdist(t(rProK), "jaccard",binary = TRUE)
pcoa <- cmdscale(ProK.dist)

#nMDS
nMDS <- metaMDS(t(rProK),"jaccard",binary = TRUE)

##PERMANOVA ProK
sink("model.output/beta.16s(second).txt")
PermDispProK <- betadisper(vegdist(t(rProK), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
anova(PermDispProK)
TukeyHSD(PermDispProK)
plot(TukeyHSD(PermDispProK))
#Difference in variance between East and West 
adonis(vegdist(t(rProK), "jaccard",binary = TRUE)~latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
adonis.pair(vegdist(t(rProK), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(rProK),latlongdata$sitecode)])
#p<0.001
sink()

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


#### Now let's run the analysis by domains 


##First COI.Euk

#nMDS
nMDS <- metaMDS(t(rCOI.euk),"jaccard",binary = TRUE)

##PERMANOVA COI euk
PermDispCOI.euk <- betadisper(vegdist(t(rCOI.euk), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(rCOI.euk),latlongdata$sitecode)])
anova(PermDispCOI.euk)
#No sig difference in multivariance homogenity
adonis(vegdist(t(rCOI.euk), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(rCOI.euk),latlongdata$sitecode)])
adonis.pair(vegdist(t(rCOI.euk), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(rCOI.euk),latlongdata$sitecode)])
#p<0.001


#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/COI.euk.beta.pdf",width = 5,height=5)
palette(c('#51B9E0','#4bad84','#E8A016'))
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


#Here we are using smoothed GAMs in the function ordisurf across the nMDS space
##Some thoughts on nMDS ordisurf - https://www.fromthebottomoftheheap.net/2011/06/10/what-is-ordisurf-doing/

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(rCOI.euk),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(rCOI.euk),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(rCOI.euk),envdat$sitecode)])


#Ordisurf - temp
pdf("figures/COI.euk.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.euk.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - SSS
pdf("figures/COI.euk.SSS.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.SSS<- ordisurf(nMDS,envmatrix$SSS,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()
sink("model.output/smoothedGAM/COI.euk.SSS.txt")
summary(surfout.SSS)
sink()

#Ordisurf - ChlA
pdf("figures/COI.euk.Chl.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.Chl<- ordisurf(nMDS,envmatrix$chla_3yrmean,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.euk.Chl.txt")
summary(surfout.Chl)
sink()


#Ordisurf - impact
pdf("figures/COI.euk.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/COI.euk.impact.txt")
summary(surfout.ipt)
sink()


##Lets try and relate the environmental variables to the beta diversity 
#we are using a distance based redudnancy analysis (dbRDA)

#dbRDA
model.full <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean + temp, data = envmatrix)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.COI.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
sink()


model.3.1 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean + temp, data = envmatrix)
model.3.2 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean + temp, data = envmatrix)
model.3.3 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + impact + temp, data = envmatrix)
model.3.4 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean, data = envmatrix)

model.2.1 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + impact, data = envmatrix)
model.2.2 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean, data = envmatrix)
model.2.3 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ SSS + temp, data = envmatrix)
model.2.4 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean, data = envmatrix)
model.2.5 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ impact + temp, data = envmatrix)
model.2.6 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~ chla_3yrmean + temp, data = envmatrix)

model.1.1 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~temp, data = envmatrix)
model.1.2 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~impact, data = envmatrix)
model.1.3 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~chla_3yrmean, data = envmatrix)
model.1.4 <- dbrda(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE)~SSS, data = envmatrix)

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


#Lets replace the upset plot with a bar plot since the order is junk
#barplot(expressionInput,col="black",names=F,yaxt="n",xaxt="n")
#axis(2,cex.axis=1.2,lwd=2)
#axis(1,at=-1:18,labels=F,cex=2,lwd=2,lwd.ticks=0)

pdf("figures/upset.dbRDA.COI.euk.pdf",width = 9,height=5.5)
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
pdf("figures/Figure 2/bar.COI.euk.pdf",width = 6,height = 4)
par(mar=c(3,3,1,1))
barplot(expressionInput,names=F,col="black",space=rep(0.5,15),yaxt='n',ylim=c(0,25))
axis(2,at=seq(0,25,5),labels=seq(0,25,5))
dev.off()




#testing mantel and partial mantel tests

#Lets create a geogrphical distance matrix for all sites  
dist <- read.csv("metadata/distance.csv",row.names = 1)
dist <- as.matrix(dist)
upperTriangle(dist) <- lowerTriangle(dist,byrow = TRUE)
dist <- dist[match(colnames(rCOI.euk),colnames(dist)),]
dist <- dist[,match(colnames(rCOI.euk),colnames(dist))]
dist <- as.dist(dist)

#Lets create a false 2d world for calculations
dist.2d <- read.csv("metadata/distance.csv",row.names = 1)
dist.2d <- as.matrix(dist.2d)[,1]
dist.2d <- dist.2d[match(colnames(rCOI.euk),names(dist.2d))]
dist.2d <- data.frame("Y"=rep(1,18),"X"=dist.2d)

#Lets create distance matrices for each env parameter
SST.dist <- dist(tempdat$TempAvr[match(colnames(rCOI.euk),tempdat$Site.Code)])
SSS.dist <- dist(envdat$SSS[match(colnames(rCOI.euk),envdat$sitecode)])
ChlA.dist <- dist(envdat$chla_3yrmean[match(colnames(rCOI.euk),envdat$sitecode)])
Impact.dist <- dist(envdat$impact[match(colnames(rCOI.euk),envdat$sitecode)])


#Mantel test
m.COI.euk.SST <- mantel.randtest(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),SST.dist)
m.COI.euk.SSS <- mantel.randtest(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),SSS.dist)
m.COI.euk.ChlA <- mantel.randtest(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),ChlA.dist)
m.COI.euk.impact <- mantel.randtest(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),Impact.dist)

##Here we run a partial mantel test
###but SST and distance are highly autocorrelated! Need to correct for this as below
p.COI.euk.SST <- mantel.partial(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),SST.dist,dist)
p.COI.euk.SSS <- mantel.partial(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.COI.euk.ChlA <- mantel.partial(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.COI.euk.impact <- mantel.partial(vegdist(t(rCOI.euk),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.COI.euk.SST.c <- msr(m.COI.euk.SST,lw1, 10000)
m.COI.euk.SSS.c <- msr(m.COI.euk.SSS,lw1, 10000)
m.COI.euk.ChlA.c <- msr(m.COI.euk.ChlA,lw1, 10000)
m.COI.euk.impact.c <- msr(m.COI.euk.impact,lw1, 10000)


mantel.p.COI.euk.out <- c(p.COI.euk.SST$statistic,p.COI.euk.SST$signif,
                      p.COI.euk.SSS$statistic,p.COI.euk.SSS$signif,
                      p.COI.euk.ChlA$statistic,p.COI.euk.ChlA$signif,
                      p.COI.euk.impact$statistic,p.COI.euk.impact$signif)

mantel.c.COI.euk.out <- unname(c(m.COI.euk.SST.c$obs-m.COI.euk.SST.c$expvar["Expectation"],
                             m.COI.euk.SST.c$pvalue,
                             m.COI.euk.SSS.c$obs-m.COI.euk.SSS.c$expvar["Expectation"],
                             m.COI.euk.SSS.c$pvalue,
                             m.COI.euk.ChlA.c$obs-m.COI.euk.ChlA.c$expvar["Expectation"],
                             m.COI.euk.ChlA.c$pvalue,
                             m.COI.euk.impact.c$obs-m.COI.euk.impact.c$expvar["Expectation"],
                             m.COI.euk.impact.c$pvalue))

# Then 18S Euk

#nMDS
nMDS <- metaMDS(t(r18S.euk),"jaccard",binary = TRUE)

##PERMANOVA 18S euk
PermDisp18S.euk <- betadisper(vegdist(t(r18S.euk), "jaccard",binary = TRUE),latlongdata$PERMori[match(colnames(r18S.euk),latlongdata$sitecode)])
anova(PermDisp18S.euk)
#No sig difference in multivariance homogenity
adonis(vegdist(t(r18S.euk), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(r18S.euk),latlongdata$sitecode)])
adonis.pair(vegdist(t(r18S.euk), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(r18S.euk),latlongdata$sitecode)])
#p<0.001


#Plot
palette(c('#51B9E0','#4bad84','#E8A016'))
pdf("figures/18S.euk.beta.pdf",width = 5,height=5)
palette(c('#51B9E0','#4bad84','#E8A016'))
par(mfrow=c(1,1))
par(mar=c(3,3,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="")
ordiellipse(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
            kind = "sd",
            draw = "polygon",
            lty=0,
            col=1:3)
#ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],draw = "polygon",col=1:3,border=F)
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


#Here we are using smoothed GAMs in the function ordisurf across the nMDS space
##Some thoughts on nMDS ordisurf - https://www.fromthebottomoftheheap.net/2011/06/10/what-is-ordisurf-doing/

#first we make up a matrix containing the data we need ordered the same as our OTU table
envmatrix <- envdat[match(colnames(r18S.euk),envdat$sitecode),match(c("SSS","impact","chla_3yrmean"),colnames(envdat))]
envmatrix <- cbind(envmatrix,tempdat$TempAvr[match(colnames(r18S.euk),tempdat$Site.Code)])
colnames(envmatrix)[4] <- "temp"
rownames(envmatrix) <- as.character(envdat$sitecode[match(colnames(r18S.euk),envdat$sitecode)])


#Ordisurf - temp
pdf("figures/18S.euk.temp.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.temp<- ordisurf(nMDS,envmatrix$temp,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.euk.temp.txt")
summary(surfout.temp)
sink()

#Ordisurf - SSS
pdf("figures/18S.euk.SSS.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.SSS<- ordisurf(nMDS,envmatrix$SSS,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()
sink("model.output/smoothedGAM/18S.euk.SSS.txt")
summary(surfout.SSS)
sink()

#Ordisurf - ChlA
pdf("figures/18S.euk.Chl.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.Chl<- ordisurf(nMDS,envmatrix$chla_3yrmean,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.euk.Chl.txt")
summary(surfout.Chl)
sink()


#Ordisurf - impact
pdf("figures/18S.euk.impact.ordisurf.pdf",width = 2.5,height=2.5)
par(mar=c(1,1,1,1))
plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,main="",cex=0,xlab="",ylab="",xaxt='n',yaxt='n')
points(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=1.5)
surfout.ipt<- ordisurf(nMDS,envmatrix$impact,method="REML",add = TRUE,col="darkred",lwd=3)
dev.off()

sink("model.output/smoothedGAM/18S.euk.impact.txt")
summary(surfout.ipt)
sink()


##Lets try and relate the environmental variables to the beta diversity 
#we are using a distance based redudnancy analysis (dbRDA)

#dbRDA
model.full <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean + temp, data = envmatrix)

plot(model.full)
anova(model.full,permutations = 10000)
sink("model.output/dbRDA.anova.18S.txt")
anova(model.full,by="margin",permutations = 10000)
RsquareAdj(model.full)
sink()


model.3.1 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean + temp, data = envmatrix)
model.3.2 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean + temp, data = envmatrix)
model.3.3 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + impact + temp, data = envmatrix)
model.3.4 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + impact + chla_3yrmean, data = envmatrix)

model.2.1 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + impact, data = envmatrix)
model.2.2 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + chla_3yrmean, data = envmatrix)
model.2.3 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ SSS + temp, data = envmatrix)
model.2.4 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ impact + chla_3yrmean, data = envmatrix)
model.2.5 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ impact + temp, data = envmatrix)
model.2.6 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~ chla_3yrmean + temp, data = envmatrix)

model.1.1 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~temp, data = envmatrix)
model.1.2 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~impact, data = envmatrix)
model.1.3 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~chla_3yrmean, data = envmatrix)
model.1.4 <- dbrda(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE)~SSS, data = envmatrix)

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


#Lets replace the upset plot with a bar plot since the order is junk
#barplot(expressionInput,col="black",names=F,yaxt="n",xaxt="n")
#axis(2,cex.axis=1.2,lwd=2)
#axis(1,at=-1:18,labels=F,cex=2,lwd=2,lwd.ticks=0)

pdf("figures/upset.dbRDA.18S.euk.pdf",width = 9,height=5.5)
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
pdf("figures/Figure 2/bar.18S.euk.pdf",width = 6,height = 4)
par(mar=c(3,3,1,1))
barplot(expressionInput,names=F,col="black",space=rep(0.5,15),yaxt='n',ylim=c(0,25))
axis(2,at=seq(0,25,5),labels=seq(0,25,5))
dev.off()




#testing mantel and partial mantel tests

#Lets create a geogrphical distance matrix for all sites  
dist <- read.csv("metadata/distance.csv",row.names = 1)
dist <- as.matrix(dist)
upperTriangle(dist) <- lowerTriangle(dist,byrow = TRUE)
dist <- dist[match(colnames(r18S.euk),colnames(dist)),]
dist <- dist[,match(colnames(r18S.euk),colnames(dist))]
dist <- as.dist(dist)

#Lets create a false 2d world for calculations
dist.2d <- read.csv("metadata/distance.csv",row.names = 1)
dist.2d <- as.matrix(dist.2d)[,1]
dist.2d <- dist.2d[match(colnames(r18S.euk),names(dist.2d))]
dist.2d <- data.frame("Y"=rep(1,18),"X"=dist.2d)

#Lets create distance matrices for each env parameter
SST.dist <- dist(tempdat$TempAvr[match(colnames(r18S.euk),tempdat$Site.Code)])
SSS.dist <- dist(envdat$SSS[match(colnames(r18S.euk),envdat$sitecode)])
ChlA.dist <- dist(envdat$chla_3yrmean[match(colnames(r18S.euk),envdat$sitecode)])
Impact.dist <- dist(envdat$impact[match(colnames(r18S.euk),envdat$sitecode)])


#Mantel test
m.18S.euk.SST <- mantel.randtest(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),SST.dist)
m.18S.euk.SSS <- mantel.randtest(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),SSS.dist)
m.18S.euk.ChlA <- mantel.randtest(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),ChlA.dist)
m.18S.euk.impact <- mantel.randtest(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),Impact.dist)

##Here we run a partial mantel test
###but SST and distance are highly autocorrelated! Need to correct for this as below
p.18S.euk.SST <- mantel.partial(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),SST.dist,dist)
p.18S.euk.SSS <- mantel.partial(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),SSS.dist,dist)
p.18S.euk.ChlA <- mantel.partial(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),ChlA.dist,dist)
p.18S.euk.impact <- mantel.partial(vegdist(t(r18S.euk),method = "jaccard",binary = TRUE),Impact.dist,dist)


#Lets try the corrected version from Crabot et al. 2019 DOI:10.1111/2041-210X.13141

#first we need to get some weights and neighbors to run the correction
nb1 <- graph2nb(gabrielneigh(as.matrix(dist.2d)), sym = T) 
lw1 <- nb2listw(nb1)

#now we run the msr corrected test on the output of the mantel w/ 10000 perms
m.18S.euk.SST.c <- msr(m.18S.euk.SST,lw1, 10000)
m.18S.euk.SSS.c <- msr(m.18S.euk.SSS,lw1, 10000)
m.18S.euk.ChlA.c <- msr(m.18S.euk.ChlA,lw1, 10000)
m.18S.euk.impact.c <- msr(m.18S.euk.impact,lw1, 10000)


mantel.p.18S.euk.out <- c(p.18S.euk.SST$statistic,p.18S.euk.SST$signif,
                          p.18S.euk.SSS$statistic,p.18S.euk.SSS$signif,
                          p.18S.euk.ChlA$statistic,p.18S.euk.ChlA$signif,
                          p.18S.euk.impact$statistic,p.18S.euk.impact$signif)

mantel.c.18S.euk.out <- unname(c(m.18S.euk.SST.c$obs-m.18S.euk.SST.c$expvar["Expectation"],
                                 m.18S.euk.SST.c$pvalue,
                                 m.18S.euk.SSS.c$obs-m.18S.euk.SSS.c$expvar["Expectation"],
                                 m.18S.euk.SSS.c$pvalue,
                                 m.18S.euk.ChlA.c$obs-m.18S.euk.ChlA.c$expvar["Expectation"],
                                 m.18S.euk.ChlA.c$pvalue,
                                 m.18S.euk.impact.c$obs-m.18S.euk.impact.c$expvar["Expectation"],
                                 m.18S.euk.impact.c$pvalue))






####======4.3 Phyla Analyses====####

#Here we are analysing some select phyla to check the main effects 

MetazoanPhyla <- c("Arthropoda","Cnidaria","Mollusca","Chordata","Porifera")
ProtistPhyla <- c("Dinoflagellata","Ciliophora","Protalveolata","Cercozoa","Prymnesiophyceae")
BacterialPhyla <- c("Proteobacteria","Bacteroidetes","Cyanobacteria","Epsilonbacteraeota","Verrucomicrobia")

summary <- data.frame("phyla"=c(MetazoanPhyla,ProtistPhyla,BacterialPhyla),
                      "group"=c(rep("Metazoa",5),rep("Protist",5),rep("Bacteria",5)),
                      "n"=rep(NA,15),
                      "F.disp"=rep(NA,15),
                      "p.disp"=rep(NA,15),
                      "F.PERM"=rep(NA,15),
                      "p.PERM"=rep(NA,15)
                      )

#For each phyla
#plot nMDS & show results 
pdf("figures/phyla.data/metazoans.ecoregions.pdf",height=6.5,width = 9)
par(mfrow=c(2,3),mar=c(3,3,1,1))

for (phyla in MetazoanPhyla){
  
  loopdata <- rCOI.met[rownames(rCOI.met) %in% COIassignments.RDP$V1[COIassignments.RDP$V12==phyla],]
  summary$n[summary$phyla==phyla] <- length(loopdata[,1])
  
  ##PERMANOVA 
  loopPermDisp <- betadisper(vegdist(t(loopdata), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  temp <- anova(loopPermDisp)
  summary$F.disp[summary$phyla==phyla] <- round(temp$`F value`[1],3)
  loop.disp.p <-temp$`Pr(>F)`[1]
  #No sig difference in multivariate homogenity
  temp <- adonis(vegdist(t(loopdata), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  loop.PERM.p <-temp$aov.tab$`Pr(>F)`[1]
  adonis.pair(vegdist(t(loopdata), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  #p<0.001
  
  summary$p.disp[summary$phyla==phyla] <- round(loop.disp.p,3)
  summary$p.PERM[summary$phyla==phyla] <- round(loop.PERM.p,3)
  summary$F.PERM[summary$phyla==phyla] <- round(temp$aov.tab$F.Model[1],3)
  
  
  nMDS <- metaMDS(vegdist(t(loopdata),"jaccard",binary=TRUE))
  looppcoa <- pcoa(vegdist(t(loopdata),"jaccard",binary=TRUE))
  
  
  
  #Plot nMDS
  palette(c('#51B9E0','#4bad84','#E8A016'))
  
  if (round(nMDS$stress,3)>0.001){
  plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=0,xlab="",ylab="",main=phyla)
  ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
             lty=1,
             lwd=2,
             col=1:3)
  ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
  #text(-0.5,0.3,labels=paste0("Stress = ",round(nMDS$stress,3)))
  legend("topleft",paste0("Stress = ",round(nMDS$stress,3)), bty="n") 
  #legend("topright",paste0("disp=",round(loop.disp.p,3)," PERM=",round(loop.PERM.p,3)), bty="n") 
  #legend("bottomright",phyla, bty="n") 
  
  }else{

  plot(looppcoa$vectors[,1],looppcoa$vectors[,2],col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=0,xlab="",ylab="",main=phyla)
   
    W <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="W",1:2]
    E <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="E",1:2]
    S <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="S",1:2]
    
    polygon(W[chull(W),],col=adjustcolor('#51B9E0', alpha.f = 0.40),border = NA)
    polygon(E[chull(E),],col=adjustcolor('#E8A016', alpha.f = 0.40),border = NA)
    polygon(S[chull(S),],col=adjustcolor('#4bad84', alpha.f = 0.40),border = NA)
    
  points(looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m",1],
         looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m",2],
         col=latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)][latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m"],
         pch=16,cex=3)
  points(looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n",1],
         looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n",2],
         col=latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)][latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n"],
         pch=17,cex=3)
  text(looppcoa$vectors[,1],looppcoa$vectors[,2],labels=rownames(looppcoa$vectors))
  
  #legend("topleft",paste0("PCoA"), bty="n") 
  #legend("topright",paste0("disp=",round(loop.disp.p,3)," PERM=",round(loop.PERM.p,3)), bty="n") 
  #legend("bottomright",phyla, bty="n") 
  
  }
}

dev.off()

pdf("figures/phyla.data/protists.ecoregions.pdf",height=6.5,width = 9)
par(mfrow=c(2,3),mar=c(3,3,1,1))

for (phyla in ProtistPhyla){
  
  loopdata <- r18S.pts[rownames(r18S.pts) %in% z18Sassignments.RDP$V1[z18Sassignments.RDP$V12==phyla],]
  summary$n[summary$phyla==phyla]  <- length(loopdata[,1])
  
  ##PERMANOVA 
  loopPermDisp <- betadisper(vegdist(t(loopdata), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  temp <- anova(loopPermDisp)
  summary$F.disp[summary$phyla==phyla] <- round(temp$`F value`[1],3)
  loop.disp.p <-temp$`Pr(>F)`[1]
  #No sig difference in multivariate homogenity
  temp <- adonis(vegdist(t(loopdata), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  loop.PERM.p <-temp$aov.tab$`Pr(>F)`[1]
  adonis.pair(vegdist(t(loopdata), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  #p<0.001
  
  summary$p.disp[summary$phyla==phyla] <- round(loop.disp.p,3)
  summary$p.PERM[summary$phyla==phyla] <- round(loop.PERM.p,3)
  summary$F.PERM[summary$phyla==phyla] <- round(temp$aov.tab$F.Model[1],3)
  
  
  nMDS <- metaMDS(vegdist(t(loopdata),"jaccard",binary=TRUE))
  looppcoa <- pcoa(vegdist(t(loopdata),"jaccard",binary=TRUE))
  
  
  
  #Plot nMDS
  palette(c('#51B9E0','#4bad84','#E8A016'))
  
  if (round(nMDS$stress,3)>0.001){
    plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=0,xlab="",ylab="",main=phyla)
    ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
               lty=1,
               lwd=2,
               col=1:3)
    ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
    #text(-0.5,0.3,labels=paste0("Stress = ",round(nMDS$stress,3)))
    legend("topleft",paste0("Stress = ",round(nMDS$stress,3)), bty="n") 
    #legend("topright",paste0("disp=",round(loop.disp.p,3)," PERM=",round(loop.PERM.p,3)), bty="n") 
    #legend("bottomright",phyla, bty="n") 
    
  }else{
    
    plot(looppcoa$vectors[,1],looppcoa$vectors[,2],col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=0,xlab="",ylab="",main=phyla)
    
    W <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="W",1:2]
    E <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="E",1:2]
    S <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="S",1:2]
    
    polygon(W[chull(W),],col=adjustcolor('#51B9E0', alpha.f = 0.40),border = NA)
    polygon(E[chull(E),],col=adjustcolor('#E8A016', alpha.f = 0.40),border = NA)
    polygon(S[chull(S),],col=adjustcolor('#4bad84', alpha.f = 0.40),border = NA)
    
    points(looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m",1],
           looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m",2],
           col=latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)][latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m"],
           pch=16,cex=3)
    points(looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n",1],
           looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n",2],
           col=latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)][latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n"],
           pch=17,cex=3)
    text(looppcoa$vectors[,1],looppcoa$vectors[,2],labels=rownames(looppcoa$vectors))
    
    #legend("topleft",paste0("PCoA"), bty="n") 
    #legend("topright",paste0("disp=",round(loop.disp.p,3)," PERM=",round(loop.PERM.p,3)), bty="n") 
    #legend("bottomright",phyla, bty="n") 
    
  }
}

dev.off()


pdf("figures/phyla.data/bacteria.ecoregions.pdf",height=6.5,width = 9)
par(mfrow=c(2,3),mar=c(3,3,1,1))

for (phyla in BacterialPhyla){
  
  loopdata <- rProK[rownames(rProK) %in% ProKtaxa$OTUname[ProKtaxa$Phylum==phyla],]
  summary$n[summary$phyla==phyla]  <- length(loopdata[,1])
  
  ##PERMANOVA 
  loopPermDisp <- betadisper(vegdist(t(loopdata), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  temp <- anova(loopPermDisp)
  summary$F.disp[summary$phyla==phyla] <- round(temp$`F value`[1],3)
  loop.disp.p <-temp$`Pr(>F)`[1]
  #No sig difference in multivariate homogenity
  temp <- adonis(vegdist(t(loopdata), "jaccard",binary=TRUE)~latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  loop.PERM.p <-temp$aov.tab$`Pr(>F)`[1]
  adonis.pair(vegdist(t(loopdata), "jaccard",binary=TRUE),latlongdata$PERMori[match(colnames(loopdata),latlongdata$sitecode)])
  #p<0.001
  
  summary$p.disp[summary$phyla==phyla] <- round(loop.disp.p,3)
  summary$p.PERM[summary$phyla==phyla] <- round(loop.PERM.p,3)
  summary$F.PERM[summary$phyla==phyla] <- round(temp$aov.tab$F.Model[1],3)
  
  
  nMDS <- metaMDS(vegdist(t(loopdata),"jaccard",binary=TRUE))
  looppcoa <- pcoa(vegdist(t(loopdata),"jaccard",binary=TRUE))
  
  
  
  #Plot nMDS
  palette(c('#51B9E0','#4bad84','#E8A016'))
  
  if (round(nMDS$stress,3)>0.001){
    plot(nMDS$points,col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=0,xlab="",ylab="",main=phyla)
    ordispider(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
               lty=1,
               lwd=2,
               col=1:3)
    ordihull(nMDS,latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],
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
    #text(-0.5,0.3,labels=paste0("Stress = ",round(nMDS$stress,3)))
    legend("topleft",paste0("Stress = ",round(nMDS$stress,3)), bty="n") 
    #legend("topright",paste0("disp=",round(loop.disp.p,3)," PERM=",round(loop.PERM.p,3)), bty="n") 
    #legend("bottomright",phyla, bty="n") 
    
  }else{
    
    plot(looppcoa$vectors[,1],looppcoa$vectors[,2],col=latlongdata$PERMori[match(rownames(nMDS$points),latlongdata$sitecode)],pch=16,cex=0,xlab="",ylab="",main=phyla)
    
    W <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="W",1:2]
    E <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="E",1:2]
    S <- looppcoa$vectors[latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="S",1:2]
    
    polygon(W[chull(W),],col=adjustcolor('#51B9E0', alpha.f = 0.40),border = NA)
    polygon(E[chull(E),],col=adjustcolor('#E8A016', alpha.f = 0.40),border = NA)
    polygon(S[chull(S),],col=adjustcolor('#4bad84', alpha.f = 0.40),border = NA)
    
    points(looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m",1],
           looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m",2],
           col=latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)][latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="m"],
           pch=16,cex=3)
    points(looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n",1],
           looppcoa$vectors[latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n",2],
           col=latlongdata$PERMori[match(rownames(looppcoa$vectors),latlongdata$sitecode)][latlongdata$type[match(rownames(looppcoa$vectors),latlongdata$sitecode)]=="n"],
           pch=17,cex=3)
    text(looppcoa$vectors[,1],looppcoa$vectors[,2],labels=rownames(looppcoa$vectors))
    
    #legend("topleft",paste0("PCoA"), bty="n") 
    #legend("topright",paste0("disp=",round(loop.disp.p,3)," PERM=",round(loop.PERM.p,3)), bty="n") 
    #legend("bottomright",phyla, bty="n") 
    
  }
}

dev.off()

write.csv(summary,"figures/phyla.data/table.csv")

pdf("figures/phyla.data/Table.pdf",height=5,width = 8)
grid.table(summary,rows=rep("",15))
dev.off()


####====5.0 Distance-Decay====####
###Now lets explore the dissimilarity x distance relationship 

####======5.1 Kingdom Analyses ====####

palette(c('#D36526','#2671B4'))

##First metazoans
data <- rCOI.met

data.dist <- vegdist(t(data), "jaccard",binary = "TRUE",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA
pairwise <- melt(data.dist2)
pairwise <- pairwise[!is.na(pairwise$value) & pairwise$value!=0,]
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


#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
pairwise$impactMean <- rep(NA,length(pairwise[,1]))
pairwise$impactDiff <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  first.site <- strsplit(pairwise$comp[row],"\\.")[[1]][1]
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  second.site <- strsplit(pairwise$comp[row],"\\.")[[1]][2]
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
  pairwise$impactMean[row] <- mean(c(envdat$impact[envdat$sitecode==first.site],
                                   envdat$impact[envdat$sitecode==second.site]))
  pairwise$impactDiff[row] <- sqrt((envdat$impact[envdat$sitecode==first.site]-
                                   envdat$impact[envdat$sitecode==second.site])^2)
}

#We get rid of mixed samples & comparisons with complete differences (we can't log transform these :( )
pairwise <- pairwise[pairwise$value!=1,]
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]

##Plot time
pairwise.nomixed$comp %in% Kmdatapair$comp

hist(pairwise.nomixed$value,breaks=20) 
pdf("figures/Overall.COI.met.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="Metazoans")
dev.off()

#Lets do some stats! 

#What about an interaction of the various impact scores?
dist.all <- Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]

#plot(lm(log10(1-pairwise$value)~pairwise$impactMean))
#mean is useless 
summary(lm(log10(1-pairwise$value)~dist.all*pairwise$impactMean))
#diff is useful 
impactDiff.1 <- pairwise$impactDiff
summary(lm(log10(1-pairwise$value)~dist.all*impactDiff.1))
summary(lm(log10(1-pairwise$value)~dist.all+impactDiff.1))
m1 <- lm(log10(1-pairwise$value)~dist.all+impactDiff.1)


plot(dist.all,log10(1-pairwise$value),pch=16)
indata <- data.frame("dist.all"=rep(seq(0,2100,10),5),"site.dissimilarity"=rep(NA,1055),"impactDiff.1"=c(rep(0,211),rep(1,211),rep(2,211),rep(3,211),rep(4,211)))
indata <- cbind(indata,predict(m1,indata,interval="confidence"))
par(mfrow=c(1,1))
palette(heat.colors(5,rev=TRUE))
dataset <- 0
plot(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],type="l",col=dataset+1)
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")

for (dataset in 1:4){
  lines(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],col=dataset+1)
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")
  text(300,-0.08,labels=paste("Impact Difference",dataset))
}


palette(heat.colors(51,rev=TRUE))
par(mfrow=c(1,2))
plot(dist.all,1-pairwise$value,col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)
plot(dist.all,log10(1-pairwise$value),col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)

plot(1:50,1:50,col=1:50,pch=16)
palette(c('#D36526','#2671B4'))

##First lets rename everything to make it easier to understand
site.dissimilarity <- pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type

par(mfrow=c(1,1))
#Lets exmaine the diagnostic plots 
#plot(lm(log10(site.dissimilarity)~distance*site.type))
##They look good
model <- summary(lm(log10(1-site.dissimilarity)~distance*site.type))
##No effect 

sink("model.output/Dist.Decay.lm.COI.met.txt")
print(summary(lm(log10(1-site.dissimilarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.dissimilarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(1-site.dissimilarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.COI.met.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/Figure3/log10Type.COI.met.dist.jacc.pdf",width = 6,height=5)
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


dist.decay.all <- cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("COI.met",length(pairwise.nomixed$value)))

##Now protists
data <- r18S.pts

data.dist <- vegdist(t(data), "jaccard",binary = "TRUE",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA
pairwise <- melt(data.dist2)
pairwise <- pairwise[!is.na(pairwise$value) & pairwise$value!=0,]
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


#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
pairwise$impactMean <- rep(NA,length(pairwise[,1]))
pairwise$impactDiff <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  first.site <- strsplit(pairwise$comp[row],"\\.")[[1]][1]
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  second.site <- strsplit(pairwise$comp[row],"\\.")[[1]][2]
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
  pairwise$impactMean[row] <- mean(c(envdat$impact[envdat$sitecode==first.site],
                                     envdat$impact[envdat$sitecode==second.site]))
  pairwise$impactDiff[row] <- sqrt((envdat$impact[envdat$sitecode==first.site]-
                                      envdat$impact[envdat$sitecode==second.site])^2)
}

#We get rid of mixed samples & comparisons with complete differences (we can't log transform these :( )
pairwise <- pairwise[pairwise$value!=1,]
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]

##Plot time
pairwise.nomixed$comp %in% Kmdatapair$comp

hist(pairwise.nomixed$value,breaks=20) 
pdf("figures/Overall.18S.pts.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="Metazoans")
dev.off()

#Lets do some stats! 

#What about an interaction of the various impact scores?
dist.all <- Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]

#plot(lm(log10(1-pairwise$value)~pairwise$impactMean))
#mean is useless 
summary(lm(log10(1-pairwise$value)~dist.all*pairwise$impactMean))
#diff is useful 
impactDiff.1 <- pairwise$impactDiff
summary(lm(log10(1-pairwise$value)~dist.all*impactDiff.1))

m1 <- lm(log10(1-pairwise$value)~dist.all*impactDiff.1)


plot(dist.all,log10(1-pairwise$value),pch=16)
indata <- data.frame("dist.all"=rep(seq(0,2100,10),5),"site.dissimilarity"=rep(NA,1055),"impactDiff.1"=c(rep(0,211),rep(1,211),rep(2,211),rep(3,211),rep(4,211)))
indata <- cbind(indata,predict(m1,indata,interval="confidence"))
par(mfrow=c(1,1))
palette(heat.colors(5,rev=TRUE))
dataset <- 0
plot(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],type="l",col=dataset+1)
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")

for (dataset in 1:4){
  lines(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],col=dataset+1)
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")
  text(300,-0.08,labels=paste("Impact Difference",dataset))
}


palette(heat.colors(51,rev=TRUE))
par(mfrow=c(1,2))
plot(dist.all,1-pairwise$value,col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)
plot(dist.all,log10(1-pairwise$value),col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)

plot(1:50,1:50,col=1:50,pch=16)
palette(c('#D36526','#2671B4'))

##First lets rename everything to make it easier to understand
site.dissimilarity <- pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type

par(mfrow=c(1,1))
#Lets exmaine the diagnostic plots 
#plot(lm(log10(site.dissimilarity)~distance*site.type))
##They look good
model <- summary(lm(log10(1-site.dissimilarity)~distance*site.type))
##No effect 

sink("model.output/Dist.Decay.lm.18S.pts.txt")
print(summary(lm(log10(1-site.dissimilarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.dissimilarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(1-site.dissimilarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.18S.pts.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/Figure3/log10Type.18S.pts.dist.jacc.pdf",width = 6,height=5)
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


dist.decay.all <- rbind(dist.decay.all,cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("18S.pts",length(pairwise.nomixed$value))))

##Finally Bacteria
data <- rProK

data.dist <- vegdist(t(data), "jaccard",binary = "TRUE",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA
pairwise <- melt(data.dist2)
pairwise <- pairwise[!is.na(pairwise$value) & pairwise$value!=0,]
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


#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
pairwise$impactMean <- rep(NA,length(pairwise[,1]))
pairwise$impactDiff <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  first.site <- strsplit(pairwise$comp[row],"\\.")[[1]][1]
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  second.site <- strsplit(pairwise$comp[row],"\\.")[[1]][2]
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
  pairwise$impactMean[row] <- mean(c(envdat$impact[envdat$sitecode==first.site],
                                     envdat$impact[envdat$sitecode==second.site]))
  pairwise$impactDiff[row] <- sqrt((envdat$impact[envdat$sitecode==first.site]-
                                      envdat$impact[envdat$sitecode==second.site])^2)
}

#We get rid of mixed samples & comparisons with complete differences (we can't log transform these :( )
pairwise <- pairwise[pairwise$value!=1,]
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]

##Plot time
pairwise.nomixed$comp %in% Kmdatapair$comp

hist(pairwise.nomixed$value,breaks=20) 
pdf("figures/Overall.16S.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="Metazoans")
dev.off()

#Lets do some stats! 

#What about an interaction of the various impact scores?
dist.all <- Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]

#plot(lm(log10(1-pairwise$value)~pairwise$impactMean))
#mean is useless 
summary(lm(log10(1-pairwise$value)~dist.all*pairwise$impactMean))
#diff is useful 
impactDiff.1 <- pairwise$impactDiff
m1 <- lm(log10(1-pairwise$value)~dist.all*impactDiff.1)


plot(dist.all,log10(1-pairwise$value),pch=16)
indata <- data.frame("dist.all"=rep(seq(0,2100,10),5),"site.dissimilarity"=rep(NA,1055),"impactDiff.1"=c(rep(0,211),rep(1,211),rep(2,211),rep(3,211),rep(4,211)))
indata <- cbind(indata,predict(m1,indata,interval="confidence"))
par(mfrow=c(1,1))
palette(heat.colors(5,rev=TRUE))
dataset <- 0
plot(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],type="l",col=dataset+1)
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")

for (dataset in 1:4){
  lines(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],col=dataset+1)
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")
  text(300,-0.08,labels=paste("Impact Difference",dataset))
}


palette(heat.colors(51,rev=TRUE))
par(mfrow=c(1,2))
plot(dist.all,1-pairwise$value,col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)
plot(dist.all,log10(1-pairwise$value),col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)

plot(1:50,1:50,col=1:50,pch=16)
palette(c('#D36526','#2671B4'))

##First lets rename everything to make it easier to understand
site.dissimilarity <- pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type

par(mfrow=c(1,1))
#Lets exmaine the diagnostic plots 
#plot(lm(log10(site.dissimilarity)~distance*site.type))
##They look good
model <- summary(lm(log10(1-site.dissimilarity)~distance*site.type))
##No effect 

sink("model.output/Dist.Decay.lm.16S.txt")
print(summary(lm(log10(site.dissimilarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.dissimilarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(1-site.dissimilarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.16S.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/Figure3/log10Type.16S.dist.jacc.pdf",width = 6,height=5)
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


dist.decay.all <- rbind(dist.decay.all,cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("16S",length(pairwise.nomixed$value))))

##Now we run the analysis on the unused datasets

##COI protists
data <- rCOI.pts

data.dist <- vegdist(t(data), "jaccard",binary = "TRUE",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA
pairwise <- melt(data.dist2)
pairwise <- pairwise[!is.na(pairwise$value) & pairwise$value!=0,]
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


#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
pairwise$impactMean <- rep(NA,length(pairwise[,1]))
pairwise$impactDiff <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  first.site <- strsplit(pairwise$comp[row],"\\.")[[1]][1]
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  second.site <- strsplit(pairwise$comp[row],"\\.")[[1]][2]
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
  pairwise$impactMean[row] <- mean(c(envdat$impact[envdat$sitecode==first.site],
                                     envdat$impact[envdat$sitecode==second.site]))
  pairwise$impactDiff[row] <- sqrt((envdat$impact[envdat$sitecode==first.site]-
                                      envdat$impact[envdat$sitecode==second.site])^2)
}

#We get rid of mixed samples & comparisons with complete differences (we can't log transform these :( )
pairwise <- pairwise[pairwise$value!=1,]
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]

##Plot time
pairwise.nomixed$comp %in% Kmdatapair$comp

hist(pairwise.nomixed$value,breaks=20) 
pdf("figures/Overall.COI.pts.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="Metazoans")
dev.off()

#Lets do some stats! 

#What about an interaction of the various impact scores?
dist.all <- Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]

#plot(lm(log10(1-pairwise$value)~pairwise$impactMean))
#mean is useless 
summary(lm(log10(1-pairwise$value)~dist.all*pairwise$impactMean))
#diff is useful 
impactDiff.1 <- pairwise$impactDiff
m1 <- lm(log10(1-pairwise$value)~dist.all+impactDiff.1)


plot(dist.all,log10(1-pairwise$value),pch=16)
indata <- data.frame("dist.all"=rep(seq(0,2100,10),5),"site.dissimilarity"=rep(NA,1055),"impactDiff.1"=c(rep(0,211),rep(1,211),rep(2,211),rep(3,211),rep(4,211)))
indata <- cbind(indata,predict(m1,indata,interval="confidence"))
par(mfrow=c(1,1))
palette(heat.colors(5,rev=TRUE))
dataset <- 0
plot(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],type="l",col=dataset+1)
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")

for (dataset in 1:4){
  lines(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],col=dataset+1)
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")
  text(300,-0.08,labels=paste("Impact Difference",dataset))
}


palette(heat.colors(51,rev=TRUE))
par(mfrow=c(1,2))
plot(dist.all,1-pairwise$value,col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)
plot(dist.all,log10(1-pairwise$value),col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)

plot(1:50,1:50,col=1:50,pch=16)
palette(c('#D36526','#2671B4'))

##First lets rename everything to make it easier to understand
site.dissimilarity <- pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type

par(mfrow=c(1,1))
#Lets exmaine the diagnostic plots 
#plot(lm(log10(1-site.dissimilarity)~distance*site.type))
##They look good
model <- summary(lm(log10(1-site.dissimilarity)~distance*site.type))
##No effect 

sink("model.output/Dist.Decay.lm.COI.pts.txt")
print(summary(lm(log10(1-site.dissimilarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.dissimilarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(1-site.dissimilarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.COI.pts.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI.pts")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/Figure3/log10Type.COI.pts.dist.jacc.pdf",width = 6,height=5)
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


dist.decay.all <- rbind(dist.decay.all,cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("COI.pts",length(pairwise.nomixed$value))))


##18S metazoans
data <- r18S.met

data.dist <- vegdist(t(data), "jaccard",binary = "TRUE",upper = FALSE)
data.dist2 <- as.matrix(data.dist)
data.dist2[upper.tri(data.dist2)] <- NA
pairwise <- melt(data.dist2)
pairwise <- pairwise[!is.na(pairwise$value) & pairwise$value!=0,]
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


#Now lets make a character that describes the type of comparison
pairwise$type <- rep(NA,length(pairwise[,1]))
pairwise$impactMean <- rep(NA,length(pairwise[,1]))
pairwise$impactDiff <- rep(NA,length(pairwise[,1]))
for (row in 1:length(pairwise[,1])){
  first <- as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][1]])
  first.site <- strsplit(pairwise$comp[row],"\\.")[[1]][1]
  second <-as.character(latlongdata$type[latlongdata$sitecode==strsplit(pairwise$comp[row],"\\.")[[1]][2]])
  second.site <- strsplit(pairwise$comp[row],"\\.")[[1]][2]
  if(length(unique(c(first,second)))>1){
    pairwise$type[row] <- "Mixed"
  } else if (unique(c(first,second))=="m"){
    pairwise$type[row] <- "Artificial"
  } else if (unique(c(first,second))=="n"){
    pairwise$type[row] <- "Natural"}
  pairwise$impactMean[row] <- mean(c(envdat$impact[envdat$sitecode==first.site],
                                     envdat$impact[envdat$sitecode==second.site]))
  pairwise$impactDiff[row] <- sqrt((envdat$impact[envdat$sitecode==first.site]-
                                      envdat$impact[envdat$sitecode==second.site])^2)
}

#We get rid of mixed samples & comparisons with complete differences (we can't log transform these :( )
pairwise <- pairwise[pairwise$value!=1,]
pairwise.nomixed <- pairwise[pairwise$type!="Mixed",]

##Plot time
pairwise.nomixed$comp %in% Kmdatapair$comp

hist(pairwise.nomixed$value,breaks=20) 
pdf("figures/Overall.18S.met.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col= "darkblue",main="Metazoans")
dev.off()

#Lets do some stats! 

#What about an interaction of the various impact scores?
dist.all <- Kmdatapair$value[match(pairwise$comp,Kmdatapair$comp)]

#plot(lm(log10(1-pairwise$value)~pairwise$impactMean))
#mean is useless 
summary(lm(log10(1-pairwise$value)~dist.all*pairwise$impactMean))
#diff is useful 
impactDiff.1 <- pairwise$impactDiff
m1 <- lm(log10(1-pairwise$value)~dist.all+impactDiff.1)


plot(dist.all,log10(1-pairwise$value),pch=16)
indata <- data.frame("dist.all"=rep(seq(0,2100,10),5),"site.dissimilarity"=rep(NA,1055),"impactDiff.1"=c(rep(0,211),rep(1,211),rep(2,211),rep(3,211),rep(4,211)))
indata <- cbind(indata,predict(m1,indata,interval="confidence"))
par(mfrow=c(1,1))
palette(heat.colors(5,rev=TRUE))
dataset <- 0
plot(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],type="l",col=dataset+1)
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
#lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")

for (dataset in 1:4){
  lines(indata$dist.all[indata$impactDiff.1==dataset],indata$fit[indata$impactDiff.1==dataset],col=dataset+1)
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$lwr[indata$impactDiff.1==dataset],col="red")
  #lines(indata$dist.all[indata$impactDiff.1==dataset],indata$upr[indata$impactDiff.1==dataset],col="red")
  text(300,-0.08,labels=paste("Impact Difference",dataset))
}


palette(heat.colors(51,rev=TRUE))
par(mfrow=c(1,2))
plot(dist.all,1-pairwise$value,col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)
plot(dist.all,log10(1-pairwise$value),col=as.numeric(cut(round(pairwise$impactDiff,digits = 2),50)),pch=16)

plot(1:50,1:50,col=1:50,pch=16)
palette(c('#D36526','#2671B4'))

##First lets rename everything to make it easier to understand
site.dissimilarity <- pairwise.nomixed$value
distance <- Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)]
site.type <- pairwise.nomixed$type

par(mfrow=c(1,1))
#Lets exmaine the diagnostic plots 
#plot(lm(log10(1-site.dissimilarity)~distance*site.type))
##They look good
model <- summary(lm(log10(1-site.dissimilarity)~distance*site.type))
##No effect 

sink("model.output/Dist.Decay.lm.18S.met.txt")
print(summary(lm(log10(1-site.dissimilarity)~distance*site.type)))
sink()

#Build CF intervals for plotting
indata <- data.frame("distance"=rep(seq(0,2100,10),2),"site.dissimilarity"=rep(NA,422),"site.type"=c(rep("Artificial",211),rep("Natural",211)))
indata <- cbind(indata,predict(lm(log10(1-site.dissimilarity)~distance*site.type),indata,interval="confidence"))


#NowPlot to visualise effect
pdf("figures/Type.18S.met.dist.jacc.pdf",width = 6,height=5)
plot(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,ylab="Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)",pch=16, col=as.factor(pairwise.nomixed$type),main="COI.pts")
legend("topright",col=1:2,legend = levels(as.factor(pairwise.nomixed$type)),pch=16)
dev.off()
#LogPLot
pdf("figures/Figure3/log10Type.18S.met.dist.jacc.pdf",width = 6,height=5)
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


dist.decay.all <- rbind(dist.decay.all,cbind(Kmdatapair$value[match(pairwise.nomixed$comp,Kmdatapair$comp)],1-pairwise.nomixed$value,rep("18S.met",length(pairwise.nomixed$value))))






#Now lets plot some comparisons
dist.decay.all <- as.data.frame(dist.decay.all)
dist.decay.all$log10Diss <- log10(as.numeric(dist.decay.all$V2))

plot(dist.decay.all$V1[dist.decay.all$V3=="COI.met"],dist.decay.all$V2[dist.decay.all$V3=="COI.met"],pch=16,col="DarkBlue",ylim=c(0,0.5))
points(dist.decay.all$V1[dist.decay.all$V3=="18S.pts"],dist.decay.all$V2[dist.decay.all$V3=="18S.pts"],pch=16,col="DarkGreen")
points(dist.decay.all$V1[dist.decay.all$V3=="16S"],dist.decay.all$V2[dist.decay.all$V3=="16S"],pch=16,col="Gold")
points(dist.decay.all$V1[dist.decay.all$V3=="COI.pts"],dist.decay.all$V2[dist.decay.all$V3=="COI.pts"],pch=16,col="DarkOrchid")
points(dist.decay.all$V1[dist.decay.all$V3=="18S.met"],dist.decay.all$V2[dist.decay.all$V3=="18S.met"],pch=16,col="DarkCyan")



plot(dist.decay.all$V1[dist.decay.all$V3=="COI.met"],dist.decay.all$log10Diss[dist.decay.all$V3=="COI.met"],pch=16,col="DarkBlue",ylim=c(-2.3,-0.3))
points(dist.decay.all$V1[dist.decay.all$V3=="18S.pts"],dist.decay.all$log10Diss[dist.decay.all$V3=="18S.pts"],pch=16,col="DarkGreen")
points(dist.decay.all$V1[dist.decay.all$V3=="16S"],dist.decay.all$log10Diss[dist.decay.all$V3=="16S"],pch=16,col="Gold")
points(dist.decay.all$V1[dist.decay.all$V3=="COI.pts"],dist.decay.all$log10Diss[dist.decay.all$V3=="COI.pts"],pch=16,col="DarkOrchid")
points(dist.decay.all$V1[dist.decay.all$V3=="18S.met"],dist.decay.all$log10Diss[dist.decay.all$V3=="18S.met"],pch=16,col="DarkCyan")

Diss <- dist.decay.all$log10Diss
Distance <- as.numeric(dist.decay.all$V1)
Taxa <- dist.decay.all$V3

m1 <- lm(Diss~Distance*Taxa)
m2 <- lm(Diss~Distance+Taxa)

AIC(m1,m2)
summary(lm(Diss~Distance+Taxa))


m1 <- lm(Diss~Distance+Taxa)

indata <- data.frame("Distance"=rep(seq(0,2100,10),5),"Taxa"=c(rep("COI.met",211),rep("18S.pts",211),rep("16S",211),rep("COI.pts",211),rep("18S.met",211)))
indata <- cbind(indata,predict(m1, indata,interval="confidence"))

pdf("figures/Figure3/DistanceDecay.All.Data.pdf",width = 8,height = 6.5)

plot(dist.decay.all$V1[dist.decay.all$V3=="COI.met"],dist.decay.all$log10Diss[dist.decay.all$V3=="COI.met"],pch=16,col="DarkBlue",ylim=c(-2.3,-0.3),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)")
points(dist.decay.all$V1[dist.decay.all$V3=="18S.pts"],dist.decay.all$log10Diss[dist.decay.all$V3=="18S.pts"],pch=16,col="DarkGreen")
points(dist.decay.all$V1[dist.decay.all$V3=="16S"],dist.decay.all$log10Diss[dist.decay.all$V3=="16S"],pch=16,col="Gold")
points(dist.decay.all$V1[dist.decay.all$V3=="COI.pts"],dist.decay.all$log10Diss[dist.decay.all$V3=="COI.pts"],pch=16,col="DarkOrchid")
points(dist.decay.all$V1[dist.decay.all$V3=="18S.met"],dist.decay.all$log10Diss[dist.decay.all$V3=="18S.met"],pch=16,col="DarkCyan")

lines(indata$Distance[indata$Taxa=="COI.met"],indata$fit[indata$Taxa=="COI.met"],col="darkblue",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="COI.met"], rev(indata$Distance[indata$Taxa=="COI.met"])),
        y = c(indata$lwr[indata$Taxa=="COI.met"], 
              rev(indata$upr[indata$Taxa=="COI.met"])),
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)

lines(indata$Distance[indata$Taxa=="18S.pts"],indata$fit[indata$Taxa=="18S.pts"],col="darkgreen",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="18S.pts"], rev(indata$Distance[indata$Taxa=="18S.pts"])),
        y = c(indata$lwr[indata$Taxa=="18S.pts"], 
              rev(indata$upr[indata$Taxa=="18S.pts"])),
        col =  adjustcolor("green3", alpha.f = 0.10), border = NA)


lines(indata$Distance[indata$Taxa=="16S"],indata$fit[indata$Taxa=="16S"],col="goldenrod1",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="16S"], rev(indata$Distance[indata$Taxa=="16S"])),
        y = c(indata$lwr[indata$Taxa=="16S"], 
              rev(indata$upr[indata$Taxa=="16S"])),
        col =  adjustcolor("goldenrod1", alpha.f = 0.10), border = NA)

lines(indata$Distance[indata$Taxa=="COI.pts"],indata$fit[indata$Taxa=="COI.pts"],col="darkorchid",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="COI.pts"], rev(indata$Distance[indata$Taxa=="COI.pts"])),
        y = c(indata$lwr[indata$Taxa=="COI.pts"], 
              rev(indata$upr[indata$Taxa=="COI.pts"])),
        col =  adjustcolor("orchid", alpha.f = 0.10), border = NA)

lines(indata$Distance[indata$Taxa=="18S.met"],indata$fit[indata$Taxa=="18S.met"],col="DarkCyan",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="18S.met"], rev(indata$Distance[indata$Taxa=="18S.met"])),
        y = c(indata$lwr[indata$Taxa=="18S.met"], 
              rev(indata$upr[indata$Taxa=="18S.met"])),
        col =  adjustcolor("cyan", alpha.f = 0.10), border = NA)


legend("topright",col=c("DarkBlue","DarkCyan","DarkGreen","DarkOrchid","goldenrod3"),legend = c("Metazoans.COI","Metazoans.18S","Protists.18S","Protists.COI","Bacteria"),pch=16)

dev.off()
pdf("figures/Figure3/Analysis by taxa/DistanceDecay.Alltaxa.Data.pdf",width = 6,height=5)
par(mar=c(4,4,1,1))

plot(dist.decay.all$V1[dist.decay.all$V3=="COI.met"],dist.decay.all$log10Diss[dist.decay.all$V3=="COI.met"],pch=16,col="DarkBlue",ylim=c(-2.3,-0.3),ylab="Log10 Compositional Similarity (1-Jaccard)",xlab="Distance Between Sites (Km)")
points(dist.decay.all$V1[dist.decay.all$V3=="18S.pts"],dist.decay.all$log10Diss[dist.decay.all$V3=="18S.pts"],pch=16,col="DarkGreen")
points(dist.decay.all$V1[dist.decay.all$V3=="16S"],dist.decay.all$log10Diss[dist.decay.all$V3=="16S"],pch=16,col="Gold")

lines(indata$Distance[indata$Taxa=="COI.met"],indata$fit[indata$Taxa=="COI.met"],col="darkblue",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="COI.met"], rev(indata$Distance[indata$Taxa=="COI.met"])),
        y = c(indata$lwr[indata$Taxa=="COI.met"], 
              rev(indata$upr[indata$Taxa=="COI.met"])),
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)

lines(indata$Distance[indata$Taxa=="18S.pts"],indata$fit[indata$Taxa=="18S.pts"],col="darkgreen",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="18S.pts"], rev(indata$Distance[indata$Taxa=="18S.pts"])),
        y = c(indata$lwr[indata$Taxa=="18S.pts"], 
              rev(indata$upr[indata$Taxa=="18S.pts"])),
        col =  adjustcolor("green3", alpha.f = 0.10), border = NA)


lines(indata$Distance[indata$Taxa=="16S"],indata$fit[indata$Taxa=="16S"],col="goldenrod1",lwd=2)
polygon(x = c(indata$Distance[indata$Taxa=="16S"], rev(indata$Distance[indata$Taxa=="16S"])),
        y = c(indata$lwr[indata$Taxa=="16S"], 
              rev(indata$upr[indata$Taxa=="16S"])),
        col =  adjustcolor("goldenrod1", alpha.f = 0.10), border = NA)



legend("topright",col=c("DarkBlue","DarkGreen","goldenrod3"),legend = c("Metazoa","Protists","Bacteria"),pch=16)

dev.off()


####======5.2 Marker Analyses ====####
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
#plot(lm(log10(site.similarity)~distance*site.type))
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
#plot(lm(log10(site.similarity)~distance*site.type))
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
#plot(lm(log10(site.similarity)~distance*site.type))
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




####====6.0 Phyla Effects====####
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
adonis.pair(vegdist(t(dCOI), "jaccard"),latlongdata$PERMori[match(colnames(dCOI),latlongdata$sitecode)])


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
#plot(lm(log10(site.similarity)~distance*site.type))
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


