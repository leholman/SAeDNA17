#############################################
####==== South African eDNA Analysis ====####
####==== Luke E. Holman====04.04.2020====####
#############################################

###Script 2 - Supplementary Script###


####====0.0 Packages====####
library(seqinr)
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







