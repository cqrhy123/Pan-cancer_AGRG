setwd("E:\\pancancer\\相关性分析") 
options(stringsAsFactors = F)
hall<- read.table("UCS-Hallmark.txt",head=T,row.names=1)
exp <- read.table("UCSexpNerve.txt",head=T,row.names=1)
exp <- exp[rownames(hall),]
library(psych)
cortest_psy_sdj <- corr.test(hall[,1:ncol(hall)], exp[,1:ncol(exp)], method = "pearson", adjust = "fdr")
cor <- cortest_psy_sdj$r
write.csv(cor," correlation matrix_UCS.csv")
pval <- cortest_psy_sdj$p
write.csv(pval,"adj_pvalue_UCS.csv")

