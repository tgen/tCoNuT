#runDNAcopyExomeV2 Rscript
#
#  runDNAcopyExomeV2 runs the CBS algorithm implemented in 
#  DNAcopy on tCoNuT exome copy number data. DNAcopy parameters
#  were optimized for in-house exome data. User may need to modify for their
#  data needs. Requires the DNAcopy from Bioconductor. 
#  https://bioconductor.org/packages/release/bioc/html/DNAcopy.html
#
#  Rscript --vanilla runDNAcopyExomeV2.R SAMPLE.cna.tsv SAMPLE.seg
#
#  INPUTS:
#	[1] *cna.tsv file from tCoNuT.
#	[2] string for name of PNG files
#
#  OUTPUTS:
#	*.seg is a SEG file with segmentation results
#	*.png is a image file with plot from DNAcopy
#
  
library('DNAcopy')
Sys.setenv("DISPLAY"="hpc-hn01:0.0")

args = commandArgs(TRUE)
cghTSV = args[1]
pngName = args[2]  

cgh=read.table(cghTSV,header=TRUE,sep="\t")
  
fName=substr(cghTSV,1,nchar(cghTSV)-4)
CNA.object=CNA(cgh$Fold.Change,cgh$Chr,cgh$Position,data.type="logratio",sampleid=fName)
  
segment.CNA.object=segment(CNA.object,alpha=0.001,verbose=1,undo.splits="sdundo",undo.SD=10,min.width=3)
  
write.table(print(segment.CNA.object),file=paste(fName,'.seg',sep=''),sep="\t",row.names=FALSE)
  
# pngName=substr(cghTSV,1,9)
png(file=paste(pngName,'.png',sep=''),width=500*4*3,height=500*2*3, res=300)
plot(segment.CNA.object,plot.type="w")
dev.off()

