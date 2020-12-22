#############################################
####==== South African eDNA Analysis ====####
####==== Luke E. Holman====04.04.2020====####
#############################################

###Script 2 - Supplementary Script###


####====0.0 Packages====####
library(seqinr)
library("vegan")
set.seed(12345)

####====1.0 Dark Diversity - are unannotated OTUs/ASVs real?====####

#### 1.1 Do the ASVs from the COI data contain a open reading frame?  
input <- "rawdata/DADA2.COI.OTUs.gz"

#a little function to count the number of stop codons
countey <- function(x){
  length(x[x=="*"])
}

DNA <- read.fasta(input)

TransTabs <- c(1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23)


output <- matrix(ncol=length(TransTabs),nrow=length(DNA))


for (codontab in 1:length(TransTabs)){
  #select codon table 
  code = TransTabs[codontab]
  #translate in all three frames
  F1 <- getTrans(DNA,frame=0,numcode = code)
  F2 <- getTrans(DNA,frame=1,numcode = code)
  F3 <- getTrans(DNA,frame=2,numcode = code)
  #make a table with all three frames
  table.n <- data.frame("F1"=unlist(lapply(X=F1,countey)),
                        "F2"=unlist(lapply(X=F2,countey)),
                        "F3"=unlist(lapply(X=F3,countey)))
  #output smallest number of stop codons from 3 frames
  output[,codontab] <- apply(table.n,1,min)
  #
  print(paste0(((codontab/length(TransTabs))*100),"%"))
}


OTU <- apply(output,1,min)
#Number of stop codons for best fit codon table 
all.best <- table(OTU)
#Number for invert. mito
all.invert <- table(output[,5])


#What about the OTUs retained in the QC dataset?
rCOI <- read.csv("cleaned/rarefied.COI.csv",row.names = 1)
DNA2 <- DNA[rownames(rCOI)]


output <- matrix(ncol=length(TransTabs),nrow=length(DNA2))

for (codontab in 1:length(TransTabs)){
  #select codon table 
  code = TransTabs[codontab]
  #translate in all three frames
  F1 <- getTrans(DNA2,frame=0,numcode = code)
  F2 <- getTrans(DNA2,frame=1,numcode = code)
  F3 <- getTrans(DNA2,frame=2,numcode = code)
  #make a table with all three frames
  table.n <- data.frame("F1"=unlist(lapply(X=F1,countey)),
                        "F2"=unlist(lapply(X=F2,countey)),
                        "F3"=unlist(lapply(X=F3,countey)))
  #output smallest number of stop codons from 3 frames
  output[,codontab] <- apply(table.n,1,min)
  #
  print(paste0(((codontab/length(TransTabs))*100),"%"))
}


OTU <- apply(output,1,min)

#Number of stop codons for best fit codon table 
QC.best <- table(OTU)
#Number for invert. mito
QC.invert <- table(output[,5])


#what is the largest number of stop codons in any OTU?
maximumstop <- max(as.numeric(names(c(all.invert,all.best,QC.invert,QC.best))))

#now we use this silly function to get the lists in order
stupidListTool <- function(x,y){
  landing <- rep(0,y+1)
  names(landing) <- 0:y
  x2 <- x[match(names(landing),names(x))]
  names(x2) <- names(landing)
  result <- Map("+", x2,landing)
  return(result)
}

plotdata <- cbind("QC.best"=unlist(stupidListTool(QC.best,maximumstop)),
                       "QC.invert"=unlist(stupidListTool(QC.invert,maximumstop)),
                       "all.best"=unlist(stupidListTool(all.best,maximumstop)),
                       "all.invert"=unlist(stupidListTool(all.invert,maximumstop)))

plotdata[is.na(plotdata)] <- 0


round(prop.table(plotdata,2),10)*100


#### 1.2 Let's subset 5% of OTUs/ASVs for checking by hand.  

DNA3 <- DNA2[sample(1:4867,round(4867*0.05))]

write.fasta(DNA3,names=names(DNA3),file.out = "rawdata/DADA2.COI.5%.OTUs.fasta",)
test <- data.frame("Names"=names(DNA3),"Seq"=unlist(getSequence(DNA3,as.string = TRUE)))
write.csv(data.frame("Names"=names(DNA3),"Seq"=unlist(getSequence(DNA3,as.string = TRUE))),"supplement/OTUs.csv")


####====2.0 Power Analysis====####


###Setup

#functions
getNichespp <- function(layout){
  world <- rep(0,length(layout))
  world[layout==sample(unique(layout),1)] <- 1
  return(world)
}

changeIncidence <- function(input){
  if(is.numeric(input)==FALSE){stop("Input value is non-numeric")}
  if(input==1){return(0)}else{
    if(input==0){return(1)}else{
      stop("No suitable input")}
  }}

#Commmunity 
#member type proportions
#perfect niche occupiuers
p.niche <- 0.7
#panmixia distribution
p.panmix <- 0.1
#random distribution
p.random <- 0.2
#community size
n.spp <- 2000


#Environment
##Number of regions
regions <- 3
#Number of sites
sites <- 18 
#Layout of 2D world
site.regions <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3)
#randomness
random <- 0.35


#Make a landing pad for outputs
output <- c()
bigoutput <-c()

par(mfrow=c(2,3),mar=c(4,4,1,1))


for (num.spp in c(20,50,100,500,1000,2000)){
  output <- c()
  n.spp <- num.spp
  
  for (simulation.n in 1:100){
    #par(mfrow=c(2,4),mar=c(1,1,1,1))
    for (variable in seq(0.3,0.44,0.02)){
      random <- variable
      
      
      #Construct World
      
      #niche spp
      obs.niche <- t(sapply(1:round(n.spp*p.niche), function(x) getNichespp(site.regions)))
      
      #panmixic spp
      obs.panmix <- t(sapply(1:round(n.spp*p.panmix), function(x) rep(1,sites)))
      
      #random spp
      obs.rand <- t(sapply(1:round(n.spp*p.random), function(x) sample(c(1,0),sites,replace=TRUE)))
      
      #combine
      sim.dat <- as.data.frame(rbind(obs.niche,obs.panmix,obs.rand))
      
      #add random noise
      #all the possible options
      allDat <- expand.grid(1:dim(sim.dat)[1],1:dim(sim.dat)[2])
      nObvs <- nrow(allDat)*ncol(allDat)
      noise <- allDat[sample(nrow(allDat),round(nrow(allDat)*random)),]
      
      for (item in 1:nrow(noise)){
        sim.dat[noise$Var1[item],noise$Var2[item]] <- changeIncidence(sim.dat[noise$Var1[item],noise$Var2[item]])
      }
      
      #Output
      PERM <- adonis(vegdist(t(sim.dat), "jaccard",binary=TRUE)~as.factor(site.regions))
      PERM$aov.tab$`Pr(>F)`[1]
      
      output <- rbind(output,c(PERM$aov.tab$`Pr(>F)`[1],random,simulation.n,n.spp))
      
      #test <- metaMDS(vegdist(t(sim.dat),"jaccard",binary=TRUE))
      #plot(test)
      #ordihull(test,site.regions)
      #legend("topright",bty="n",legend=paste("R =",random))
      #legend("topleft",bty="n",legend=paste("p =",PERM$aov.tab$`Pr(>F)`[1]))
    }
    print(simulation.n)
  }
  output <- as.data.frame(output)
  plot(output$V1~jitter(output$V2),xlab="Randomness",ylab="PERMANOVA.result",pch=16,cex=0.7,ylim=c(0,1))
  abline(h=0.05,col="red",lty=2)
  legend("topleft",legend=paste("N.spp=",n.spp),bty="n")
  bigoutput <- rbind(bigoutput,output)
}

write.csv(bigoutput,"model.output/simulation.csv")

#bigoutput <-read.csv("model.output/simulation.csv")

pdf("figures/simulationPval.pdf",height=5,width=7)
par(mfrow=c(2,3),mar=c(4,4,1,1))
for (size in unique(bigoutput$V4)){
plot(bigoutput$V1[bigoutput$V4==size]~jitter(bigoutput$V2[bigoutput$V4==size]),xlab="Randomness",ylab="PERMANOVA p value",xaxt='n',pch=16,cex=0.7,ylim=c(0,1))
axis(1,at = unique(bigoutput$V2),labels =format(unique(bigoutput$V2),digits=3),las =2)
legend("topleft",bty="n",legend=paste("nSpp =",size))
abline(h=0.05,col="red",lty=2)
abline(h=0.01,col="#301934",lty=2)
}
dev.off()


####====3.0 Control Sample Additional Analysis====####

controlDat <- read.csv("metadata/ControlData.csv")

#COI

COI.control <- read.csv("controls/COIcontrol.csv",row.names = 1)

COI.smol <- colSums(COI.control[,1:10])

COI.smol.b <- COI.control[,1:10]
COI.smol.b[COI.smol.b > 0] <- 1
COI.smol.b <- colSums(COI.smol.b)

names(COI.smol) <- gsub("L.","",names(COI.smol)) 

tapply(COI.smol, controlDat$Type[match(names(COI.smol),controlDat$ID)], FUN=sum)
tapply(COI.smol.b, controlDat$Type[match(names(COI.smol.b),controlDat$ID)], FUN=sum)


#18S

z18S.control <- read.csv("controls/18Scontrol.csv",row.names = 1)

z18S.smol <- colSums(z18S.control[,1:10])

z18S.smol.b <- z18S.control[,1:10]
z18S.smol.b[z18S.smol.b > 0] <- 1
z18S.smol.b <- colSums(z18S.smol.b)

names(z18S.smol) <- gsub("Z.","",names(z18S.smol)) 

tapply(z18S.smol, controlDat$Type[match(names(z18S.smol),controlDat$ID)], FUN=sum)
tapply(z18S.smol.b, controlDat$Type[match(names(z18S.smol),controlDat$ID)], FUN=sum)

#ProK

ProK.control <- read.csv("controls/ProKcontrol.csv",row.names = 1)

ProK.smol <- colSums(ProK.control[,1:11])

ProK.smol.b <- ProK.control[,1:11]
ProK.smol.b[ProK.smol.b > 0] <- 1
ProK.smol.b <- colSums(ProK.smol.b)


tapply(ProK.smol, controlDat$Type[match(names(ProK.smol),controlDat$ID)], FUN=sum)
tapply(ProK.smol.b, controlDat$Type[match(names(ProK.smol.b),controlDat$ID)], FUN=sum)

##Combine data

controldat <- cbind( tapply(COI.smol, controlDat$Type[match(names(COI.smol),controlDat$ID)], FUN=sum),
                     tapply(COI.smol.b, controlDat$Type[match(names(COI.smol),controlDat$ID)], FUN=sum),
                     tapply(z18S.smol, controlDat$Type[match(names(z18S.smol),controlDat$ID)], FUN=sum),
                     tapply(z18S.smol.b, controlDat$Type[match(names(z18S.smol),controlDat$ID)], FUN=sum),
                     tapply(ProK.smol, controlDat$Type[match(names(ProK.smol),controlDat$ID)], FUN=sum),
                     tapply(ProK.smol.b, controlDat$Type[match(names(ProK.smol),controlDat$ID)], FUN=sum))

controlprop <- prop.table(controldat,2)


pdf("figures/contam.pdf",width = 6,height=4)
par(mar=c(3,4,1,6),xpd=TRUE)
palette(c('#a16928','#798234','#2887a1','#d46780'))
barplot(controlprop[c(2,1,3),],ylab="Proportion",names.arg=c("COI.R","COI.S","18S.R","18S.S","16S.R","16S.S"),col=1:3)
legend(7.5,0.8,legend=c("Travel","Lab","PCR1","PCR2"),col=c(3,2,1,4),pch=15,cex=1.3,bty="n",pt.cex = 3)
dev.off()







