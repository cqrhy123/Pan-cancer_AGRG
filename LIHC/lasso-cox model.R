
library(survival)
library(survminer)
library(glmnet)
setwd("E:\\pancancer\\LGG-GEO\\LASSO")         
rt=read.table("LGGexpNervenewlasso.txt", header=T, sep="\t", check.names=F, row.names=1)     #??ȡ?????ļ?


x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lasso.lambdaLGGexpNervenewlasso.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("lasso.cvfitLGGexpNervenewlasso.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s=cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lasso.SigExpLGGexpNervenewlasso.txt",sep="\t",row.names=F,quote=F)

rt=read.table("lasso.SigExpLGGexpNervenewlasso.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox, direction="both")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multi.CoxLGGexpNervenewlasso.txt",sep="\t",row.names=F,quote=F)

riskScore=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`", "", coxGene)
outCol=c("futime", "fustat", coxGene)
riskOut=cbind(rt[,outCol], riskScore)
riskOut=cbind(id=rownames(riskOut), riskOut)
write.table(riskOut, file="riskScoreLGGexpNervenewlasso.txt", sep="\t", quote=F, row.names=F)


