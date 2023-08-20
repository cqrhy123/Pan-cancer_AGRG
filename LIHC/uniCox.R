

setwd("E:\\pancancer\\LGG-GEO\\LASSO")
library(survival)
clrs <- fpColors(box="green",line="darkblue", summary="royalblue")             #????ɭ??ͼ??ɫ
rt=read.table("LGGexpNervenew.txt",header=T,sep="\t",check.names=F,row.names=1)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
              beta=signif(coxSummary$coef[1]),
              HR=coxSummary$conf.int[,"exp(coef)"],
              HR.95L=coxSummary$conf.int[,"lower .95"],
              HR.95H=coxSummary$conf.int[,"upper .95"],
              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
write.table(outTab,file="uniCoxsLGGexpNervenew.xls",sep="\t",row.names=F,quote=F)


           
dev.off()



