###=====Fig S1====####

rm (list=ls())
gc()


library(dplyr)
library(ggpubr)
library(ggrepel)

load('data/TCGA/TMB_Treg.Rdata') 
load('data/TCGA/medianTMB_Treg.Rdata') 
load('data/TIDE/TIDE score/TIDEscore.Rdata')   
source('Code/corplot.R')
meta <- data.table::fread('data/TCGA/Survival_SupplementalTable_S1_20171025_xena_sp')
medianTMB <- medianTMB[,c(2,3,4,1)]

TIDEscore <- t(TIDEscore)
TIDEscore <- as.data.frame(TIDEscore)
colnames(TIDEscore) <- TIDEscore['Patient',]
TIDEscore <- TIDEscore[c(14,11,12,13,15,4),]
rownames <- rownames(TIDEscore)
colnames <- colnames(TIDEscore)
TIDEscore1=as.data.frame(lapply(TIDEscore,as.numeric))
rownames(TIDEscore1) <- rownames
colnames(TIDEscore1) <- colnames
TIDEscore <- as.matrix(TIDEscore1,nrow=length(rownames(TIDEscore1)),byrow=TRUE,
                       dimnames=list(rownames(TIDEscore1),colnames(TIDEscore1)))

TMB <- t(TMB)
TMB <- as.data.frame(TMB)
colnames(TMB) <- TMB['sample',]
TMB <- TMB[c(1,4:5),]
rownames <- rownames(TMB)


TMB <- apply(TMB,2,as.numeric)
TMB <- as.data.frame(TMB)
rownames(TMB) <- rownames

cancer <- data.frame(sample=colnames(TMB))
cancer <- left_join(cancer,meta[,c(1,3)],by='sample')
CA <- unique(cancer$`cancer type abbreviation`)
CA <- CA[!is.na(CA)]

TMB['cancer',] <- cancer$`cancer type abbreviation`


for(i in CA){
  TMB['group_Treg.Sig',TMB['cancer',]==i]  = ifelse(as.numeric(TMB[1,TMB['cancer',]==i]) > medianTMB[medianTMB[,'CancerType'] ==i, 1], 'High Treg.Sig', 'Low Treg.Sig')
  TMB['group_TMB',TMB['cancer',]==i]  = ifelse(as.numeric(TMB[2,TMB['cancer',]==i]) > medianTMB[medianTMB[,'CancerType'] ==i, 2], 'High TMB', 'Low TMB')
}

for(i in colnames(TMB)){
  if (TMB['group_Treg.Sig',i] == 'High Treg.Sig' & TMB['group_TMB',i] == 'High TMB' ){
    TMB['group',i] = 'HTregHT'
  } else if (TMB['group_Treg.Sig',i] == 'High Treg.Sig' & TMB['group_TMB',i] == 'Low TMB'){
    TMB['group',i] = 'HTregLT'
  } else if (TMB['group_Treg.Sig',i] == 'Low Treg.Sig' & TMB['group_TMB',i] == 'High TMB'){
    TMB['group',i] = 'LTregHT'
  } else if (TMB['group_Treg.Sig',i] == 'Low Treg.Sig' & TMB['group_TMB',i] == 'Low TMB'){
    TMB['group',i] = 'LTregLT'
  }
}

ls_coef <- list()
ls_p <- list()
gp <- unique(as.character(TMB['group',]))


for(i in gp){
  df <- TMB[1,TMB['group',] == i]
  
  
  ls_coef[[i]] <- corplot(df,TIDEscore,1,'coef',Im,meta)
  
  ls_p[[i]] <- corplot(df,TIDEscore,1,"pValue",Im,meta)
  
}

id <- intersect(colnames(TIDEscore),colnames(TMB))
mcp <- as.data.frame(TIDEscore[,id])
mcp$cells <- rownames(TIDEscore[,id])
mcp <- reshape2::melt(mcp,value.name='TIDE scores',id='cells',variable.name='id')
# mcp$`TIDE scores` <- log10(mcp$`TIDE scores`)

group <- TMB[c('group','cancer'),]
group <- as.data.frame(t(group))
group$id <- rownames(group)

mcp <- left_join(mcp,group,by='id')
mcp$group2 <- ifelse(grepl('LTreg',mcp$group),'LTreg','HTreg')
mcp$group3 <- ifelse(grepl('LT',mcp$group),'LT','HT')

ls_p <- list()

ls_compare1 <- list(c('HTregHT','HTregLT'),c('HTregHT','LTregHT'),
                    c('HTregHT','LTregLT'),c('HTregLT','LTregHT'),
                    c('HTregLT','LTregLT'),c('LTregHT','LTregLT'))

for(i in unique(mcp$cells)[1:6]){
  ls_p[[i]] <- ggplot(data = mcp[mcp$cells==i,], aes(x = group, y= `TIDE scores`, fill = group), title = i)+
    geom_boxplot()+
    theme_classic()+
    scale_fill_manual(values=c("#2E9FDF","#66CCFF","#E7B800","#FFFF00"))+
    theme(axis.line.x = element_line(size = 0.8))+
    theme(axis.line.y = element_line(size = 0.8))+
    theme(axis.text = element_text(color = "black"))+
    stat_compare_means(comparison=ls_compare1,method = "wilcox.test",label = "p.signif")+
    labs(title = i)
}


ls_p[["Low Treg vs. High Treg"]] <- ggplot(data = mcp, aes(x =cells , y= `TIDE scores`, fill = group2))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values=c("#2E9FDF","#E7B800"))+
  theme(axis.line.x = element_line(size = 0.8))+
  theme(axis.line.y = element_line(size = 0.8))+
  theme(axis.text = element_text(color = "black"))+
  stat_compare_means(method = "wilcox.test",label = "p.signif")+
  labs(title = "Low Treg vs High Treg")+
  coord_flip()


ls_p[["Low TMB vs. High TMB"]] <- ggplot(data = mcp, aes(x =cells , y= `TIDE scores`, fill = group3))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values=c("#2E9FDF","#E7B800"))+
  theme(axis.line.x = element_line(size = 0.8))+
  theme(axis.line.y = element_line(size = 0.8))+
  theme(axis.text = element_text(color = "black"))+
  stat_compare_means(method = "wilcox.test",label = "p.signif")+
  labs(title = "Low TMB vs High TMB")+
  coord_flip()



Fig.S1.A <- ls_p[7];Fig.S1.A
Fig.S1.B <- ls_p[8];Fig.S1.B
Fig.S1.C <- ggarrange(plotlist = ls_p[1:6],common.legend = T);Fig.S1.C

