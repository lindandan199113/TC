###===Figure 4===###
rm (list=ls())
gc()

library(ROCit)
library(caret)
library(dplyr)
library(survminer)
library(survival)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ROCit)
library(ggsci)
library(cancerclass)
##===loading data===###


load('results/Machine Learning/validation_set.Rdata')
load('results/Machine Learning/training_set.Rdata')
load('results/Machine Learning/test_set.Rdata')
load('results/Machine Learning/res.Rdata')
load('data/sig/Treg.Sig1.Rdata')
load('results/Machine Learning/AUC.Rdata')


###===Figure 4A===###

# Flow chart


###===Figure 4B===###


models <-  c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')

Treg.sig_1<-Treg.Sig[c(16,20,21,25,28,30)]#0.5
Treg.sig_2<-c("HLA_DPA1","HLA_DPB1","MARCKS","HLA_DQB1","HLA_DMA", "HLA_DMB")
Treg.Sig[c(16,20,21,25,28,30)]<-Treg.sig_2

for (i in 1:length(Treg.sig_1)) {
  names(training)[names(training)==Treg.sig_1[i]]<-Treg.sig_2[i]
  names(validation)[names(validation)==Treg.sig_1[i]]<-Treg.sig_2[i]
  names(test)[names(test)==Treg.sig_1[i]]<-Treg.sig_2[i]
}

roc <- lapply(1:7,function(i){
  if (!i == 7) {
    prob <- predict(res[['model']][[i]],validation[,-1],type = "prob") # use 'nb' model
    pre <- predict(res[['model']][[i]],validation[,-1]) # use 'nb' model
    test_set <- data.frame(obs = validation$response, NR = prob[,'NR'], R = prob[,'R'], pred=pre)
    roc <- rocit(score = test_set$NR,
                 class = test_set$obs,
                 negref = 'R')
  } else{
    vali <- validation[,colnames(validation) %in% c('response',Treg.Sig)]
    pData <- data.frame(class = vali$response, sample = rownames(vali),row.names = rownames(vali))
    phenoData <- new("AnnotatedDataFrame",data=pData)
    Sig.Exp <- t(vali[,-1])
    Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
    prediction <- predict(res[['model']][[i]], Sig.Exp.test,"NR", ngenes=nrow(Sig.Exp), dist = "cor")
    roc <- rocit(as.numeric(prediction@prediction[,'z']),
                 prediction@prediction[,'class_membership'],negref = 'R')
  }
})

colset <- pal_npg("nrc")(8)
colset <- colset[2:8]

p<-plot(roc[[7]], col = colset[7], 
        legend = FALSE, YIndex = F)

for(i in 1:6){
  lines(roc[[i]]$TPR~roc[[i]]$FPR, 
        col = colset[i], lwd = 2)
}

legend("bottomright", inset=c(0.5,1), xpd=TRUE,
       col = colset[1:7],
       paste(models, 'AUC', round(res[['auc']]$ROC,2)), 
       lwd = 2)



###===Figure 4C===###

validation_set <- validation[,colnames(validation) %in% c('response',Treg.Sig)]
testing_set <- test[,colnames(test) %in% c('response',Treg.Sig)]

#ROC plot
rocplot <- function(data){
  pData <- data.frame(class = data$response, sample = rownames(data),row.names = rownames(data))
  phenoData <- new("AnnotatedDataFrame",data=pData)
  Sig.Exp <- t(data[,-1])
  Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
  prediction <- predict(res[['model']][[which.max(res$auc$ROC)]], Sig.Exp.test,"NR", ngenes=nrow(Sig.Exp), dist = "cor")
  roc <- rocit(as.numeric(prediction@prediction[,'z']),
               prediction@prediction[,'class_membership'],negref = 'R')
  plot(roc,legend=F,YIndex = F,col = colset[7])
  title(deparse(substitute(data)))
  text(x=0.5,y=0.6,labels = paste0("AUC: ",round(roc$AUC,2), " (95%CI: ",round(as.numeric(ciAUC(roc)[5]),2),'-',round(as.numeric(ciAUC(roc)[6]),2),")"))
  
}

rocplot(validation_set) # validation
rocplot(testing_set) # testing



###===Figure 4D===###

#Survival Plot

survplot <- function(data){
  data_set <- data[,colnames(data) %in% c('response',Treg.Sig)]
  
  if(all(is.na(data$OS))){return(NA)} 
  
  model.name = deparse(substitute(data))
  
  pData <- data.frame(class = data_set$response, sample = rownames(data_set),row.names = rownames(data_set))
  phenoData <- new("AnnotatedDataFrame",data=pData)
  Sig.Exp <- t(data_set[,-1])
  Sig.Exp.test <- ExpressionSet(assayData=as.matrix(Sig.Exp),phenoData=phenoData)
  pre <- predict(res[['model']][[which.max(res$auc$ROC)]], Sig.Exp.test,"NR", ngenes=nrow(Sig.Exp), dist = "cor")
  data$pre <- pre@prediction[,'class_membership']
  data <- data[!is.na(data$OS),] #remove patient without OS data
  data$OS <- as.numeric(data$OS)
  data$status <- as.numeric(data$status)
  data$pre <- factor(data$pre,levels = c('R','NR'))
  
  surv <- Surv(data$OS, data$status)
  fit <- surv_fit(surv~pre,data =data)
  medianOS <- surv_median(fit)
  survdiff <- survdiff(surv~pre,data = data)
  cox <- coxph(surv~pre, data=data) 
  
  
  pic <- ggsurvplot(fit, data = data,
                    title= model.name,
                    legend.labs = c("Low risk","High risk"),
                    break.x.by = 4,ylab="OS(probability)", xlab = " Time (Months)",
                    palette = c("#2E9FDF","#E7B800"),
                    censor.shape = 3,censor.size = 1,
                    legend.title='',
                    surv.median.line = "hv"
  )
  
  pic$plot <- pic$plot + ggplot2::annotate("text",x = 20, y = 0.68, size =3,
                                           label = paste("HR :",format(round(summary(cox)$conf.int[1],2),nsmall=2))) +
    ggplot2::annotate("text",x = 20, y = 0.62, size = 3,
                      label = paste("(","95%CI: ", format(round(summary(cox)$conf.int[1,3],2),nsmall = 2),"-",format(round(summary(cox)$conf.int[1,4],2),nsmall = 2),")",sep = ""))+
    ggplot2::annotate("text",x = 20, y = 0.95, size =3,
                      label = paste("Median OS"))+
    ggplot2::annotate("text",x = 20, y = 0.89, size =3,
                      label = paste(round(medianOS$median[1],2),"months")) +
    ggplot2::annotate("text",x = 20, y = 0.82, size =3,
                      label = paste(round(medianOS$median[2],2),"months"))+
    ggplot2::annotate("text",x = 20, y = 0.56, size =3,
                      label = paste("Log-rank p:",round(1-pchisq(survdiff$chisq,1),4)))+
    ggplot2::theme(legend.text=element_text(size=8,face = "bold"))+
    theme(axis.line.x = element_line(size=0.8),
          axis.line.y = element_line(size=0.8),
    )+
    font("xy.title",face = "bold")+
    scale_x_continuous(breaks = seq(0, max(data$OS,na.rm=T), by = 12))
  
  return(pic$plot)
}

survplot(validation)
survplot(test)
