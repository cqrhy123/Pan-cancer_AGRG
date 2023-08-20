


subsetlasso <- function(train, test, coxpcut=0.05, HRcuthigh=1.5, unicox=T, HRcutlow=0.5, nloop=20, R2= 0.5, BootStrap=F, lambda="min", nThreads=4){
  commgene <- intersect(colnames(train), colnames(test))      
  train <- as.data.frame(train)[,commgene]
  test <- as.data.frame(test)[,commgene]
  filename1 <- paste0("Unitrain_","coxpcut-",coxpcut,"HRcuthigh-",HRcuthigh,"_HRcutlow-",HRcutlow,".Rda")
  if(unicox==T){
   
    Sur1 <- Surv(train$times+1, train$status)                 
    UniCox <- function(x){
      fml <- as.formula(paste0('Sur1~', x))
      uniMod <- coxph(fml, data = train)
      SuniMod <- summary(uniMod)
      HR <- round(SuniMod$coefficients[,2],2)
      P_Value <- round(SuniMod$coefficients[,5],5)
      UniCox <- data.frame(Characteristics = x,
                           Hazard_Ratio = HR,
                           P_Value = P_Value)
      return(UniCox)
    }
    
    library("parallel")
    cl <- makeCluster(nThreads)
    clusterExport(cl, c("coxph","train","Surv"))
    Uni <- ldply(parLapply(cl, setdiff(colnames(train), c("times","status")), fun = UniCox), .parallel = T)
    stopCluster(cl)
   
    Uni_pvalue <- na.omit(Uni[Uni$P_Value <= coxpcut,])       
    Uni_pvalue <- Uni_pvalue[Uni_pvalue$Hazard_Ratio>=HRcuthigh & Uni_pvalue$Hazard_Ratio<=20 | Uni_pvalue$Hazard_Ratio<=HRcutlow,] # 去掉HR异常值的变量
    save(Uni_pvalue, file = filename1)
  }else{
    load(filename1)
  }
  MultiNames <- as.character(Uni_pvalue$Characteristics) 
  dd1 <- train[,colnames(train) %in% MultiNames]
  Train <- data.frame(times=train$times, status=train$status, dd1)
 
  out <- list() 
  length(out)<-32
  names(out) <- c("nloop","train_number","lasso_vars","lassovars_coef","nlasso_vars","lassovars_Intertrain_AUC","lassovars_Interverify_AUC",
                  "lassovars_Train_AUC","lassovars_Extertest_AUC","lassovars_Intertrain.cindex","lassovars_Interverify.cindex","lassovars_Train.cindex",
                  "lassovars_Extertest.cindex","lassovars_Intertrain.riskscore","lassovars_Interverify.riskscore",
                  "lassovars_Train.riskscore","lassovars_Extertest.riskscore","subset_vars","subsetvars_coef","nsubset_vars","subsetvars_Intertrain_AUC",
                  "subsetvars_Interverify_AUC","subsetvars_Train_AUC","subsetvars_Extertest_AUC","subsetvars_Intertrain.cindex","subsetvars_Interverify.cindex","subsetvars_Train.cindex",
                  "subsetvars_Extertest.cindex","subsetvars_Intertrain.riskscore",
                  "subsetvars_Interverify.riskscore","subsetvars_Train.riskscore","subsetvars_Extertest.riskscore")
  cl <- makeCluster(nThreads)
  registerDoParallel(cl)
  clusterExport(cl, c("cv.glmnet","createResample","createDataPartition","Surv","regsubsets","coxph","survivalROC","rcorrcens"))
  ff <- foreach(i=1:nloop, .combine=rbind,.multicombine=F, .verbose=T )%dopar% {
    out$nloop[i] <- i
    if(BootStrap==T){
      perm <- c(createResample(Train$status, times = 1, list = FALSE))
    }else{
      perm <- c(createDataPartition(Train$status, p = R2, list = FALSE)) 
    }
    out$train_number[i] <- paste(perm, collapse = "#")
    TTrain<-Train[perm,]
    TTest<-Train[-perm,]
    formula <- as.formula(Surv(TTrain$times+1, TTrain$status)~.)
    mod <- model.matrix(formula, TTrain)
    set.seed(i)
    cv.fitA <- cv.glmnet(mod,Surv(TTrain$times+1, TTrain$status), type.measure = "deviance", family = "cox")
    if(lambda=="min"){
      coeA <- coef(cv.fitA$glmnet.fit, s=cv.fitA$lambda.min)   
    }
    if(lambda=="1se"){
      coeA <- coef(cv.fitA$glmnet.fit, s=cv.fitA$lambda.1se)  
    }
    tA.active.coef <- coeA[which(coeA[,1]!=0)]
    tA.name <- row.names(coeA)[which(coeA[,1]!=0)]
    
    if(length(tA.name) >= 3 ){
 
    out$lasso_vars[i] <- paste(tA.name,collapse = "#")
    out$nlasso_vars[i] <- length(tA.name)
    
    Sur <- Surv(TTrain$times+1, TTrain$status) 
    fml <- as.formula(paste0('Sur~', paste(tA.name,collapse = '+')))          # lasso筛选的变量建立cox模型
    train_Cox <- coxph(fml, data = TTrain)
  
    tA.riskscore<-predict(train_Cox, type="risk", newdata=TTrain)# 内部训练集预测
    tA.riskscore<-tA.riskscore[!is.na(tA.riskscore)&is.finite(tA.riskscore)]

    tB.riskscore<-predict(train_Cox, type="risk", newdata=TTest)      # 内部验证集预测
    tB.riskscore<-tB.riskscore[!is.na(tB.riskscore)&is.finite(tB.riskscore)]
    
    test.riskscore<-predict(train_Cox, type="risk", newdata=test)             # 外部测试集预测
    test.riskscore<-test.riskscore[!is.na(test.riskscore)&is.finite(test.riskscore)]
    
    train.riskscore<-predict(train_Cox, type="risk", newdata=train)           # 全部训练集预测
    train.riskscore<-train.riskscore[!is.na(train.riskscore)&is.finite(train.riskscore)]
    
    Coef<-coef(train_Cox)
    out$lassovars_coef[i]<-paste(round(Coef,3), collapse = "#")

    tA_cindex <-1-rcorrcens(Surv(TTrain[names(tA.riskscore),]$times+1, TTrain[names(tA.riskscore),]$status)~tA.riskscore)[1]
    
    
    tB_cindex <- 1-rcorrcens(Surv(TTest[names(tB.riskscore),]$times+1, TTest[names(tB.riskscore),]$status)~tB.riskscore)[1]
    
    
    test_cindex <-1-rcorrcens(Surv(test[names(test.riskscore),]$times+1,test[names(test.riskscore),]$status)~test.riskscore)[1]
    
    
    train_cindex <- 1-rcorrcens(Surv(train[names(train.riskscore),]$times+1, train[names(train.riskscore),]$status)~train.riskscore)[1]
    
    
    out$lassovars_Intertrain.cindex[i] <- round(tA_cindex,3)
    out$lassovars_Interverify.cindex[i] <- round(tB_cindex,3)
    out$lassovars_Extertest.cindex[i] <- round(test_cindex,3)
    out$lassovars_Train.cindex[i] <- round(train_cindex,3)
  
    test_auc_1y <- survivalROC(Stime=test[names(test.riskscore),]$times/365*12, status=test[names(test.riskscore),]$status, marker = test.riskscore, 
                               predict.time =12, method="KM")$AUC
    test_auc_3y <- survivalROC(Stime=test[names(test.riskscore),]$times/365*12, status=test[names(test.riskscore),]$status, marker = test.riskscore, 
                               predict.time =36, method="KM")$AUC
    test_auc_5y <- survivalROC(Stime=test[names(test.riskscore),]$times/365*12, status=test[names(test.riskscore),]$status, marker = test.riskscore, 
                               predict.time =60, method="KM")$AUC
    
    train_auc_1y <- survivalROC(Stime=train[names(train.riskscore),]$times/365*12, status=train[names(train.riskscore),]$status, marker = train.riskscore, 
                                predict.time =12, method="KM")$AUC
    train_auc_3y <- survivalROC(Stime=train[names(train.riskscore),]$times/365*12, status=train[names(train.riskscore),]$status, marker = train.riskscore, 
                                predict.time =36, method="KM")$AUC
    train_auc_5y <- survivalROC(Stime=train[names(train.riskscore),]$times/365*12, status=train[names(train.riskscore),]$status, marker = train.riskscore, 
                                predict.time =60, method="KM")$AUC
    tA_auc_1y <- survivalROC(Stime=TTrain[names(tA.riskscore),]$times/365*12, status=TTrain[names(tA.riskscore),]$status, marker = tA.riskscore, 
                             predict.time =12, method="KM")$AUC
    tA_auc_3y <- survivalROC(Stime=TTrain[names(tA.riskscore),]$times/365*12, status=TTrain[names(tA.riskscore),]$status, marker = tA.riskscore, 
                             predict.time =36, method="KM")$AUC
    tA_auc_5y <- survivalROC(Stime=TTrain[names(tA.riskscore),]$times/365*12, status=TTrain[names(tA.riskscore),]$status, marker = tA.riskscore, 
                             predict.time =60, method="KM")$AUC
    
    tB_auc_1y <- survivalROC(Stime=TTest[names(tB.riskscore),]$times/365*12, status=TTest[names(tB.riskscore),]$status, marker = tB.riskscore, 
                             predict.time =12, method="KM")$AUC
    tB_auc_3y <- survivalROC(Stime=TTest[names(tB.riskscore),]$times/365*12, status=TTest[names(tB.riskscore),]$status, marker = tB.riskscore, 
                             predict.time =36, method="KM")$AUC
    tB_auc_5y <- survivalROC(Stime=TTest[names(tB.riskscore),]$times/365*12, status=TTest[names(tB.riskscore),]$status, marker = tB.riskscore, 
                             predict.time =60, method="KM")$AUC
    
    out$lassovars_Extertest_AUC[i] <- paste(round(c(test_auc_1y, test_auc_3y, test_auc_5y),3), collapse = "#")
    out$lassovars_Train_AUC[i] <- paste(round(c(train_auc_1y, train_auc_3y, train_auc_5y),3), collapse = "#")
    out$lassovars_Interverify_AUC[i] <- paste(round(c(tB_auc_1y, tB_auc_3y, tB_auc_5y),3), collapse = "#")
    out$lassovars_Intertrain_AUC[i] <- paste(round(c(tA_auc_1y, tA_auc_3y, tA_auc_5y),3), collapse = "#")
    out$lassovars_Extertest.riskscore[i] <- paste(round(test.riskscore,3), collapse = "#")
    out$lassovars_Train.riskscore[i] <- paste(round(train.riskscore,3), collapse = "#")
    out$lassovars_Intertrain.riskscore[i] <- paste(round(tA.riskscore,3), collapse = "#")
    out$lassovars_Interverify.riskscore[i] <- paste(round(tB.riskscore,3), collapse = "#")
    
    if(length(tA.name) <= 10 ){ 
      bestname<-tA.name
    }
    if(length(tA.name) > 10 &length(tA.name) <= 30){ 
   
    Sur <- Surv(TTrain$times+1, TTrain$status) 
    fml <- as.formula(paste0('Sur~', paste(tA.name,collapse = '+')))
    regsub <- regsubsets(fml, data = TTrain, method = "exhaustive",really.big=T,nvmax=10)
    regsubsum <- summary(regsub)
    bestname <- tA.name[regsubsum$which[,-1][which.min(regsubsum$cp),]] 
    }
    if(length(tA.name) > 30){  
      
      Sur <- Surv(TTrain$times+1, TTrain$status) 
      fml <- as.formula(paste0('Sur~', paste(tA.name,collapse = '+')))
      regsub <- regsubsets(fml, data = TTrain, nvmax=10,method = "forward",really.big=T)
      regsubsum <- summary(regsub)
      bestname <- tA.name[regsubsum$which[,-1][which.min(regsubsum$cp),]] 
    }
    
    out$subset_vars[i] <- paste(bestname, collapse = "#")
    out$nsubset_vars[i] <- length(bestname)
      
    fml <- as.formula(paste0('Sur~', paste(bestname,collapse = '+'))) 
    train_Cox <- coxph(fml, data =TTrain)
    
    Coef <- coef(train_Cox)
    out$subsetvars_coef[i] <- paste(round(Coef,3), collapse = "#")
      
    tA.riskscore<-predict(train_Cox, type="risk", newdata=TTrain)
    tA.riskscore<-tA.riskscore[!is.na(tA.riskscore)&is.finite(tA.riskscore)]
    
    tB.riskscore<-predict(train_Cox, type="risk", newdata=TTest)      
    tB.riskscore<-tB.riskscore[!is.na(tB.riskscore)&is.finite(tB.riskscore)]
    
    test.riskscore<-predict(train_Cox, type="risk", newdata=test)             
    test.riskscore<-test.riskscore[!is.na(test.riskscore)&is.finite(test.riskscore)]
    
    train.riskscore<-predict(train_Cox, type="risk", newdata=train)          
    train.riskscore<-train.riskscore[!is.na(train.riskscore)&is.finite(train.riskscore)]
  
    tA_cindex <-1-rcorrcens(Surv(TTrain[names(tA.riskscore),]$times+1, TTrain[names(tA.riskscore),]$status)~tA.riskscore)[1]
    
    
    tB_cindex <- 1-rcorrcens(Surv(TTest[names(tB.riskscore),]$times+1, TTest[names(tB.riskscore),]$status)~tB.riskscore)[1]
    
    
    test_cindex <-1-rcorrcens(Surv(test[names(test.riskscore),]$times+1,test[names(test.riskscore),]$status)~test.riskscore)[1]
    
    
    train_cindex <- 1-rcorrcens(Surv(train[names(train.riskscore),]$times+1, train[names(train.riskscore),]$status)~train.riskscore)[1]
    
    
    out$subsetvars_Intertrain.cindex[i] <- round(tA_cindex,3)
    out$subsetvars_Interverify.cindex[i] <- round(tB_cindex,3)
    out$subsetvars_Extertest.cindex[i] <- round(test_cindex,3)
    out$subsetvars_Train.cindex[i] <- round(train_cindex,3)
    
    test_auc_1y <- survivalROC(Stime=test[names(test.riskscore),]$times/365*12, status=test[names(test.riskscore),]$status, marker = test.riskscore, 
                               predict.time =12, method="KM")$AUC
    test_auc_3y <- survivalROC(Stime=test[names(test.riskscore),]$times/365*12, status=test[names(test.riskscore),]$status, marker = test.riskscore, 
                               predict.time =36, method="KM")$AUC
    test_auc_5y <- survivalROC(Stime=test[names(test.riskscore),]$times/365*12, status=test[names(test.riskscore),]$status, marker = test.riskscore, 
                               predict.time =60, method="KM")$AUC
    
    train_auc_1y <- survivalROC(Stime=train[names(train.riskscore),]$times/365*12, status=train[names(train.riskscore),]$status, marker = train.riskscore, 
                                predict.time =12, method="KM")$AUC
    train_auc_3y <- survivalROC(Stime=train[names(train.riskscore),]$times/365*12, status=train[names(train.riskscore),]$status, marker = train.riskscore, 
                                predict.time =36, method="KM")$AUC
    train_auc_5y <- survivalROC(Stime=train[names(train.riskscore),]$times/365*12, status=train[names(train.riskscore),]$status, marker = train.riskscore, 
                                predict.time =60, method="KM")$AUC
    tA_auc_1y <- survivalROC(Stime=TTrain[names(tA.riskscore),]$times/365*12, status=TTrain[names(tA.riskscore),]$status, marker = tA.riskscore, 
                             predict.time =12, method="KM")$AUC
    tA_auc_3y <- survivalROC(Stime=TTrain[names(tA.riskscore),]$times/365*12, status=TTrain[names(tA.riskscore),]$status, marker = tA.riskscore, 
                             predict.time =36, method="KM")$AUC
    tA_auc_5y <- survivalROC(Stime=TTrain[names(tA.riskscore),]$times/365*12, status=TTrain[names(tA.riskscore),]$status, marker = tA.riskscore, 
                             predict.time =60, method="KM")$AUC
    
    tB_auc_1y <- survivalROC(Stime=TTest[names(tB.riskscore),]$times/365*12, status=TTest[names(tB.riskscore),]$status, marker = tB.riskscore, 
                             predict.time =12, method="KM")$AUC
    tB_auc_3y <- survivalROC(Stime=TTest[names(tB.riskscore),]$times/365*12, status=TTest[names(tB.riskscore),]$status, marker = tB.riskscore, 
                             predict.time =36, method="KM")$AUC
    tB_auc_5y <- survivalROC(Stime=TTest[names(tB.riskscore),]$times/365*12, status=TTest[names(tB.riskscore),]$status, marker = tB.riskscore, 
                             predict.time =60, method="KM")$AUC
   
      out$subsetvars_Extertest_AUC[i] <- paste(round(c(test_auc_1y,test_auc_3y,test_auc_5y),3),collapse = "#")
      out$subsetvars_Train_AUC[i] <- paste(round(c(train_auc_1y,train_auc_3y,train_auc_5y),3),collapse = "#")
      out$subsetvars_Interverify_AUC[i] <- paste(round(c(tB_auc_1y,tB_auc_3y,tB_auc_5y),3),collapse = "#")
      out$subsetvars_Intertrain_AUC[i] <- paste(round(c(tA_auc_1y,tA_auc_3y,tA_auc_5y),3),collapse = "#")
      out$subsetvars_Extertest.riskscore[i] <- paste(round(test.riskscore,3),collapse = "#")
      out$subsetvars_Train.riskscore[i] <- paste(round(train.riskscore,3),collapse = "#")
      out$subsetvars_Intertrain.riskscore[i] <- paste(round(tA.riskscore,3),collapse = "#")
      out$subsetvars_Interverify.riskscore[i] <- paste(round(tB.riskscore,3),collapse = "#")
      out
    }
    }

  stopCluster(cl)
  n <- dim(ff)[1]
  ss <- lapply(1:n,function(n){do.call("rbind", ff[n,])[,ncol(do.call("rbind", ff[n,]))]})
  ss <- ldply(ss)
  return(ss)
}

coletdata <- function(data){
  data1 <-data.frame(as.numeric(data$nloop),
                   do.call(rbind, lapply(strsplit(data$subsetvars_Interverify_AUC,"#"),as.numeric)),
                   do.call(rbind, lapply(strsplit(data$subsetvars_Extertest_AUC,"#"),as.numeric)),
                   do.call(rbind, lapply(strsplit(data$subsetvars_Train_AUC,"#"),as.numeric)),
                   as.numeric(data$subsetvars_Interverify.cindex),
                   as.numeric(data$subsetvars_Extertest.cindex),
                   as.numeric(data$subsetvars_Train.cindex))
  colnames(data1) <- c("nloop","Interverify_AUC.1y","Interverify_AUC.3y","Interverify_AUC.5y",
                       "Extertest_AUC.1y","Extertest_AUC.3y","Extertest_AUC.5y",
                       "Train_AUC.1y","Train_AUC.3y","Train_AUC.5y",
                       "Interverify.cindex","Extertest.cindex","Train.cindex")
  
  data1$cindex.mean <- apply(data1[,c("Interverify.cindex", "Extertest.cindex","Train.cindex")],1,mean)
  data1$auc.1y.mean <- apply(data1[,c("Interverify_AUC.1y", "Extertest_AUC.1y","Train_AUC.1y")],1,mean)
  data1$auc.3y.mean <- apply(data1[,c("Interverify_AUC.3y", "Extertest_AUC.3y","Train_AUC.3y")],1,mean)
  data1$auc.5y.mean <- apply(data1[,c("Interverify_AUC.5y", "Extertest_AUC.5y","Train_AUC.5y")],1,mean)
  data1$auc.interverify.mean <- apply(data1[,c("Interverify_AUC.1y","Interverify_AUC.3y","Interverify_AUC.5y")],1,mean)
  data1$auc.extertest.mean <- apply(data1[,c("Extertest_AUC.1y","Extertest_AUC.3y","Extertest_AUC.5y")],1,mean)
  data1$auc.train.mean <- apply(data1[,c("Train_AUC.1y","Train_AUC.3y","Train_AUC.5y")],1,mean)
  data1
}

bestdata <- function(data, screenresults, testauc.1y, testauc.3y, testauc.5y, trainauc.1y, trainauc.3y, trainauc.5y,
                     test.cindex=0.65,train.cindex=0.65,n=1){
              data <- data[data$Extertest.cindex>=test.cindex & data$Train.cindex>=train.cindex, ]
              data <- data[data$Extertest_AUC.1y>=testauc.1y & data$Extertest_AUC.3y>=testauc.3y & data$Extertest_AUC.5y>=testauc.5y,]
              data <- data[data$Train_AUC.1y>=trainauc.1y & data$Train_AUC.3y>=trainauc.3y & data$Train_AUC.5y>=trainauc.5y,]
              x <- rownames(data[order(data$auc.extertest.mean, decreasing = T),])[n]
              best.results <- screenresults[x,]
              best.results
}
