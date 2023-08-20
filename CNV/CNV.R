library(dplyr)
setwd("E:\\pancancer\\cnv")                 
rt=read.table("Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",header=F,sep="\t",check.names=F,row.names = 1) 
rt=as.data.frame(rt)
row=read.table("genes.txt",header=F,sep="\t",check.names=F)
row2=row[,1]  
newrt=filter(rt, rownames(rt) %in% row2) 
newrt2=rbind(rt[1,],newrt)
write.table(newrt2,file="New_CNV_Matrix.txt",sep="\t",col.names=F,row.names=T,quote=F)
rm(list = ls())

inputFile="New_CNV_Matrix.txt"     
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)   
GAIN=rowSums(rt> 0)     
LOSS=rowSums(rt< 0)      
GAIN=GAIN/ncol(rt)*100     
LOSS=LOSS/ncol(rt)*100     
data=cbind(GAIN, LOSS)
data=as.data.frame(data)
NONE=100-data$GAIN-data$LOSS
data=cbind(data,NONE)
data=data[order(data[,"GAIN"],decreasing = T),]
gene=rownames(data)
data=cbind(gene,data)
write.table(data,"UCS-CNV.txt",sep = "\t",col.names = T,row.names = F,quote = F)
citation (package = "dplyr" )
