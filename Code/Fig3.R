####===Figure 3===####
rm (list=ls())
gc()


###=== GSVA ===###
library(GSVA)
library(GSEABase)
library(data.table)
source('Code/corplot.R')
source('Code/bnplot.R')

load('data/TCGA/bulkExpMatrix2.Rdata')
load('results/GSVA/gsvaTregSig2.Rdata')
load('results/GSVA/gsvaH2.Rdata')

####===Fig. 3A===####
library(circlize)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(grImport2)
library(gridBase)
library(dplyr)

meta <- fread('data/TCGA/Survival_SupplementalTable_S1_20171025_xena_sp')
RM <- fread('data/ImmuneGenes/Repressed Markers.csv')  
RM$Genes <- RM$`Repressed Markers`

coef <- corplot(gsvaSig, bulkExpMatrix2, mode = 2, return = 'coef',Im = RM, meta = meta)
p <- corplot(gsvaSig, bulkExpMatrix2, mode = 2, return = 'pValue',Im = RM, meta = meta)

Fig.3A <-  bnplot(t(coef), t(p), mode = 2, size = c(2,4), col = c("#4783B4", "white", "#8D4A2A"), limit = 1)

Fig.3A


####===Fig. 3B===####

load(file=paste0('data/TCGA/median_Is_score_Treg.Rdata'))
Fig.3B <- ggplot(median_Is_score,aes(x=mIss,y=mSig))+
  geom_point(data = median_Is_score,aes(x=mIss,y=mSig),color="#3B7FC7",size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  geom_text_repel(aes(mIss, mSig, label =  CancerType))+
  geom_smooth(method='lm',se=T,show.legend=F,linetype=1,size = 0.6,color='#333333')+
  stat_cor(method = "spearman",show.legend = F)+
  xlab('median Immunosuppressed score')+
  ylab("median GSVA of Treg.Sig")
Fig.3B


####===Fig. 3C===####

coef <- corplot(gsvaSig,gsvaH,1,'coef',RM,meta)
p <- corplot(gsvaSig,gsvaH,1,"pValue",RM,meta)

mean <- apply(coef,1,mean)

mean <- mean[order(mean,decreasing = T)]

abs <- abs(mean)

abs <- abs[order(abs,decreasing = T)]

order <- mean[ names(mean) %in% names(abs)[1:10] ]

order <- names(order)

order <- rev(order)

Fig.3C <- bnplot(coef[order,],p[order,],mode=2,size=c(1,6),col=c('#3B7FC7',"yellow",'#D41E2A'),limit=1)
Fig.3C

####===Fig. 3D and 3E===####
load(file=paste0('data/TCGA/medianITH_Treg.Rdata'))

Fig.3D <- ggplot(medianITH,aes(x=mITH,y=mSig))+
  geom_point(data = medianITH,aes(x=mITH,y=mSig),color="#3B7FC7",size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  geom_text_repel(aes(mITH, mSig, label =  CancerType))+
  #facet_wrap(~CancerType,scales = 'free')+
  # geom_smooth(method='lm',se=F,show.legend=F,linetype=1,color='darkred',size = 0.6)+
  geom_smooth(method='lm',se=T,show.legend=F,linetype=1,size = 0.6,color='#333333')+
  stat_cor(method = "spearman",show.legend = F)+
  xlab('median ITH')+
  ylab("median GSVA of Treg.Sig")

Fig.3D


#Fig.3E#
load(file=paste0('data/TCGA/medianTMB_Treg.Rdata'))
Fig.3E <- ggplot(medianTMB,aes(x=log10mTMB,y=mSig))+
  geom_point(data = medianTMB,aes(x=log10mTMB,y=mSig),color="#3B7FC7",size=3)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  geom_text_repel(aes(log10mTMB, mSig, label =  CancerType))+
  #facet_wrap(~CancerType,scales = 'free')+
  # geom_smooth(method='lm',se=F,show.legend=F,linetype=1,color='darkred',size = 0.6)+
  geom_smooth(method='lm',se=T,show.legend=F,linetype=1,size = 0.6,color='#333333')+
  stat_cor(method = "spearman",show.legend = F)+
  xlab('median log10TMB')+
  ylab("median GSVA of Treg.Sig")
Fig.3E




