#plotCGH <- function(cghTSV,gainedTSV,lossedTSV,pngName){
#
#  plotCGH is an R function for plotting whole genome copy number from tCoNuT.
#  Individual chromosomes are plotted in a figure with 6 x 4 layout.
#
#  Rscript --vanilla plotCGH.R SAMPLE.cna.tsv SAMPLE.amp.tsv SAMPLE.del.tsv SAMPLE.cna
#
#  ARGUMENTS:
#       [1] *.cna.tsv file from tCoNuT
#       [2] *.amp.tsv file from tCoNuT
#	[3] *.del.tsv file from tCoNuT
#	[4] string for name PNG file
#
#  OUTPUTS:
#       Image file named args[4].png
#
# *  [2010] - [2016] Translational Genomics Research Institute (TGen)
# *  All Rights Reserved.
# *
# * Major Contributor(s):
#    Jessica Aldrich
# * Minor Contributor(s):


Sys.setenv("DISPLAY"="hpc-hn01:0.0")

args <- commandArgs(TRUE)

cghTSV=args[1]
gainedTSV=args[2]
lossedTSV=args[3]
pngName=args[4]

#read in TSV files
cgh<-read.table(cghTSV,header=TRUE,sep="\t")
amp<-read.table(gainedTSV,header=TRUE,sep="\t")
del<-read.table(lossedTSV,header=TRUE,sep="\t")

#convert to MB
amp[,2]=amp[,2]/1e6
del[,2]=del[,2]/1e6
cgh[,2]=cgh[,2]/1e6

#create PNG file
fname=paste(pngName,'.png',sep="")
#pdf(file=fname,width=16,height=12)
png(file=fname,width=500*6*3,height=500*4*3,res=300)
par(mfrow=c(6,4))

for (i in 1:unique(cgh[,1])){
  tmpCGH <- subset(cgh, Chr == i)
  tmpAMP <- subset(amp, Chr == i)
  tmpDEL <- subset(del, Chr == i)
  plot(tmpCGH$Position,tmpCGH$Fold.Change,
       type="l",
       ylim=c(-2,2),
       col="black",
       main=paste("Chromosome ",i),
       ylab="Log2(T/N)",
       xlab="Physical Position (Mb)",
       xlim=c(floor(min(tmpCGH$Position)),ceiling(max(tmpCGH$Position))),
       xaxt="n",
       cex.lab=1.35,
       cex.main=1.5,
       cex.axis=1.15
       )  
  points(tmpAMP$Position,tmpAMP$Fold.Change,
         col="dark red"
         )
  points(tmpDEL$Position,tmpDEL$Fold.Change,
         col="dark green"
         )
  #x=seq(floor(min(tmpCGH$Position)),ceiling(max(tmpCGH$Position)),25)
  x=seq(0,250,25)
  axis(1,at=x)
}

invisible(dev.off())

