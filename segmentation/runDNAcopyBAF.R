#runDNAcopyBAF Rscript
#
#  runDNAcopyBAF runs the CBS algorithm implemented in DNAcopy on BAF (baf.txt) 
#  from parseMergeVCF.pl. DNAcopy parameters were optimized for in-house BAF data. 
#  User may need to modify for their data needs.  Requires the DNAcopy from Bioconductor. 
#  https://bioconductor.org/packages/release/bioc/html/DNAcopy.html
#
#  Rscript --vanilla runDNAcopyBAF baf.txt SAMPLE.baf
#
#  INPUTS:
#       [1] baf.txt file from parseMergeVCF.pl.
#       [2] string for name of PNG files
#
#  OUTPUTS:
#       *.seg is a SEG file with segmentation results
#       *.png is a image file with plot from DNAcopy
#
 
Sys.setenv("DISPLAY"="hpc-hn01:0.0")

args <- commandArgs(TRUE)
bafTXT=args[1]
fName=args[2]
 
library('DNAcopy')
baf=read.table(bafTXT,header=TRUE,sep="\t")
baf$BAF=abs(0.5-baf$BAF)
  
##fName=substr(bafTXT,1,nchar(bafTXT)-4)
CNA.object=CNA(baf$BAF,baf$Chr,baf$Position,data.type="logratio",sampleid=fName)
  
segment.CNA.object=segment(CNA.object,min.width=5)
  
write.table(print(segment.CNA.object),file=paste(fName,'.seg',sep=''),sep="\t",row.names=FALSE)
  
#pngName=substr(cghTSV,1,9)
png(file=paste(fName,'.seg.png',sep=''),width=500*4*3,height=500*2*3, res=300)
plot(segment.CNA.object,plot.type="w")
dev.off()

