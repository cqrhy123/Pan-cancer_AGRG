

library(corrplot)                                           
setwd("E:\\pancancer\\corrplot")        
rt=read.table("Pancancercorexp.txt",sep="\t",header=T,row.names=1,check.names=F)    


pdf("corrplot.pdf",height=10,width=10)            
par(oma=c(0.5,1,1,1.2))
M=cor(t(rt))
col3=colorRampPalette(c("#1E90FF","#FFFACD","#FF0000"))
corrplot(M, order = "original", type = "upper", tl.cex = 0.6, col = col3(100), tl.col = "black", tl.pos = "upper")

dev.off()
