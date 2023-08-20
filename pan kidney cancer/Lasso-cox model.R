
if(!require("survival")) BiocManager::install("survival",update = F,ask = F)
if(!require("survivalROC")) BiocManager::install("survivalROC",update = F,ask = F)
if(!require("pheatmap")) BiocManager::install("pheatmap",update = F,ask = F)
if(!require("ggplot2")) BiocManager::install("ggplot2",update = F,ask = F)
if(!require("survminer")) BiocManager::install("survminer",update = F,ask = F)
if(!require("DESeq2")) BiocManager::install("DESeq2",update = F,ask = F)
if(!require("preprocessCore")) BiocManager::install("preprocessCore",update = F,ask = F)
if(!require("limma")) BiocManager::install("limma",update = F,ask = F)
list.of.packages <- c("ggsci", 
                      "doParallel",
                      "Hmisc",
                      "parallel",
                      "foreach",
                      "risksetROC", 
                      "cowplot", 
                      "tidyverse",
                      "plyr",
                      "randomForestSRC",
                      "ggRandomForests",
                      "caret",
                      "glmnet",
                      "leaps",
                      "MASS",
                      "timeROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("glmnet")
library("survival")
library("survminer")
library("survivalROC")
library("ggsci")
library("tidyverse")
library("cowplot")
library("pheatmap")
library("ggplot2")
library("foreach")
library("doParallel")
library("randomForestSRC")
library("ggRandomForests")
library("Hmisc")
library("risksetROC")
library("plyr")
library("preprocessCore")
library("caret")
library("MASS")
library("leaps")
library("timeROC")


setwd("E:\\pancancer\\·º°©\\pan kidney cancer")
rm(list=ls())


clin_dd=read.table("pan_kidney cancer.txt",header=T,sep="\t",row.names=1) 
clin_dd <- as.data.frame(clin_dd)
dd <- t(clin_dd[,3:length(clin_dd)])       

dd1 <- DESeq2::varianceStabilizingTransformation(round(dd))
save(dd1,file = "dd1_vst.rda")
load("dd1_vst.rda")


dd1 <- preprocessCore::normalize.quantiles(as.matrix(dd1))
rownames(dd1) <- rownames(dd)
colnames(dd1) <- colnames(dd)
dd <- as.data.frame(t(dd1))

clin_dd <- data.frame(times=clin_dd$times, status=clin_dd$status,dd)


names <- rownames(clin_dd)
set.seed(123)                  
samp <- createDataPartition(clin_dd$status, p = 0.6, list = FALSE)    
train <- clin_dd[samp,]

test <- clin_dd[-samp,]        

source("subsetlasso.R")       

screenresults <- subsetlasso(train, test, coxpcut=0.05, HRcuthigh=1.5, unicox=T, HRcutlow=0.5, nloop=20, R2= 0.5, BootStrap=F, lambda="min", nThreads=4)


best.results<- bestdata(data=collectdata, screenresults=screenresults, testauc.1y=0.6, testauc.3y=0.6, testauc.5y=0.6, 
                        trainauc.1y=0.6, trainauc.3y=0.6, trainauc.5y=0.6,
                        test.cindex=0.6,train.cindex=0.6,n=c(1:5))

save(screenresults, file = "screenresults.rda")
save(best.results, file = "best.results.rda")


library("timeROC")
ttest <- test
ttrain <- train
ttest$riskscore <- as.numeric(strsplit(best.results$subsetvars_Extertest.riskscore,"#")[[1]])
ttrain$riskscore <- as.numeric(strsplit(best.results$subsetvars_Train.riskscore,"#")[[1]])
ROC.ttest <- timeROC(T=ttest$times/365,delta=ttest$status,
                     marker=ttest$riskscore, cause=1,
                     times=quantile(ttest$times/365, probs=seq(0.2,0.9,0.01)),
                     iid=TRUE)
ROC.ttrain<-timeROC(T=ttrain$times/365, delta=ttrain$status,
                    marker=ttrain$riskscore, cause=1,
                    times=quantile(ttrain$times/365, probs=seq(0.2,0.9,0.01)),
                    iid=TRUE)

plotAUCcurve(ROC.ttest, conf.int=TRUE, col="#0073C2FF")
plotAUCcurve(ROC.ttrain, conf.int=TRUE, col="#EFC000FF",add=TRUE)
legend("topleft", c("test","train"), col=c("#0073C2FF","#EFC000FF"), lty=1, lwd=2)


res.cut <- surv_cutpoint(ttrain, time = "times",
                         event = "status",
                         variables = "riskscore")
res.cut <- res.cut$cutpoint[[1]]


Data <- ttrain

risk <- as.vector(ifelse(Data$riskscore >= res.cut,"high","low"))
Data$risk <- risk
Sur <- Surv(Data$times/365*12, Data$status)
sfit <- survfit(Sur ~ risk, data=Data)  
ggsurvplot(sfit, 
           conf.int=F,    
           #fun="pct",
           pval=TRUE,
           palette = "jco",
           pval.method = T,
           risk.table =T, 
           ncensor.plot = T,
           surv.median.line="hv",
           legend.labs=c("high risk","low risk"))+
  labs(x = "Month")


t1 <- 12
t2 <- 36
t3 <- 60

data <- ttrain
data <- ttest 


surt1 <- survivalROC(Stime=data$times/365*12, status=data$status, marker = data$riskscore, 
                     predict.time =t1, method="KM")
surt2 <- survivalROC(Stime=data$times/365*12, status=data$status, marker = data$riskscore, 
                     predict.time =t2, method="KM")
surt3 <- survivalROC(Stime=data$times/365*12, status=data$status, marker = data$riskscore, 
                     predict.time =t3, method="KM")


plot(surt1$FP, surt1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#8F7700FF', 
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
lines(x=surt2$FP, y=surt2$TP, lwd=2, type="l", col="#0073C2FF")
lines(x=surt3$FP, y=surt3$TP, lwd=2, type="l", col="#EFC000FF")
legend.paste <- c(paste0(t1," month AUC=", round(surt1$AUC,3)),
                  paste0(t2," month AUC=", round(surt2$AUC,3)),
                  paste0(t3," month AUC=", round(surt3$AUC,3)))
legend("bottomright", legend.paste, fill=c('#8F7700FF',"#0073C2FF","#EFC000FF"), bty="n", cex=.8,
       border = NA, y.intersp=1, x.intersp=0.2 )
abline(0,1)

Data <- ttrain
Data <- ttest 

bestvars <- strsplit(best.results$subset_vars,"#")[[1]]

fp <- Data$riskscore
names(fp) <- rownames(Data)
fp<-fp[fp<=10]
fp_dat <- data.frame(s=1:length(fp), v=as.numeric(sort(fp)))
fp_dat$Risk <- ifelse(fp_dat$v >= res.cut,"high","low")
sur_dat <- data.frame(s=1:length(fp),
                      t=Data[names(sort(fp)),'times']/365*12,
                      e=Data[names(sort(fp)),'status']) 
sur_dat$Status <- as.factor(ifelse(sur_dat$e==0,'alive','death'))
exp_dat <- Data[names(sort(fp)), which(colnames(Data) %in% bestvars)]

plot.point <- ggplot(fp_dat, aes(x=s,y=v))+
  geom_point(aes(col=Risk), size=0.5)+
  geom_segment(aes(x = sum(fp_dat$Risk=="low"),
                   y = min(fp_dat$v), 
                   xend = sum(fp_dat$Risk=="low"), 
                   yend = res.cut),linetype="dashed")+
  geom_segment(aes(x=0, y=res.cut,
                   xend=sum(fp_dat$Risk=="low"),
                   yend=res.cut),linetype="dashed")+
  geom_text(aes(x=sum(fp_dat$Risk=="low")/2,
                y=res.cut+5,
                label=paste0("Cutoff: ",round(res.cut,3))),
            col ="black",size = 4,alpha=0.8)+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  labs(y="Risk score", fill="Risk")+
  scale_color_jco()

plot.sur <- ggplot(sur_dat,aes(x=s,y=t))+
  geom_point(aes(col=Status), size=0.5)+
  geom_vline(aes(xintercept=sum(fp_dat$Risk=="low")), linetype="dashed")+
  scale_colour_jco()+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(limits = c(0,NA), expand = c(0,0))+
  scale_fill_discrete(labels=c("Alive", "Dead"))+
  labs(y="Survival time(months)")

mycolors <- colorRampPalette(c("royalblue","tomato"), bias = 1.5)(150)
tmp<-t(scale(exp_dat))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
plot.h <- pheatmap::pheatmap(tmp,col=,show_colnames = F, cluster_cols = F)

plot_grid(plot.point, plot.sur, plot.h$gtable,
          labels = c("A", "B","C"),
          label_x=0,
          label_y=1,
          align = 'v',ncol = 1,axis="t")


load("Unitrain_coxpcut-0.05HRcuthigh-1.5_HRcutlow-0.5.Rda")
MultiNames <- as.character(Uni_pvalue$Characteristics) 
dd1 <- train[,colnames(train) %in% MultiNames]
Train <- data.frame(times=train$times,status=train$status,dd1)
seed <- as.numeric(best.results$nloop)
train_num <- c(as.numeric(strsplit(best.results$train_number,"#")[[1]]))
TTrain<-Train[train_num,]
formula <- as.formula(Surv(TTrain$times+1, TTrain$status)~.)
mod <- model.matrix(formula, TTrain)
set.seed(seed)
cv.fitA <- cv.glmnet(mod,Surv(TTrain$times+1, TTrain$status), type.measure = "deviance", family = "cox")
plot(cv.fitA)
plot(cv.fitA$glmnet.fit)
coeA <- coef(cv.fitA$glmnet.fit, s=cv.fitA$lambda.min)  
tA.active.coef <- coeA[which(coeA[,1]!=0)]
(tA.name <- row.names(coeA)[which(coeA[,1]!=0)])              
(lasso_vars <- strsplit(best.results$lasso_vars,"#")[[1]])    


if(as.numeric(best.results$nlasso_vars[[1]])>10){
  Sur <- Surv(train[train_num,]$times+1, train[train_num,]$status) 
  fml <- as.formula(paste0('Sur~',paste(tA.name,collapse = '+')))
  regsub <- regsubsets(fml,data = train[train_num,],nvmax=length(tA.name))
  regsubsum <- summary(regsub)
  print((bestname <- tA.name[regsubsum$which[,-1][which.min(regsubsum$cp),]]))
  print((subset_vars <- strsplit(best.results$subset_vars,"#")[[1]]))
  plot(regsub, scale="Cp")
}else{
  stop("The number of lasso_vars is less than 10!")
}








