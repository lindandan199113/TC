###==Figure 1===###
rm (list=ls())
gc()

#====loading data===####

library(Seurat)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggsci)
library(RColorBrewer)
source('Code/R_rainclouds.R')



####====Fig. 1A and 1B===####
cohort <- "SKCM_GSE120575_aPD1aCTLA4"
load(paste0("data/",cohort,'/scRNA.Rdata'))
load(paste0('data/',cohort,'/phenotype.Rdata'))
load(paste0('data/CIBERSORTx/CIBERSORTx_result/SKCM_GSE120575_aPD1aCTLA4_CIBERSORTx.Rdata'))

cell <- phenotype[phenotype$Celltype_major == 'Treg' ,'Cell']
response <- phenotype[phenotype$Celltype_major == 'Treg' ,'Response']
names(response) <- phenotype[phenotype$Celltype_major == 'Treg' ,'Cell']
scRNA <- scRNA[,colnames(scRNA) %in% cell]
scRNA <- as.data.frame(scRNA)

rownames(CIBERSORTx)<-CIBERSORTx$Mixture
results <- CIBERSORTx[rownames(CIBERSORTx)%in% cell,c("Mixture","Treg")]

####tsne
set.seed(1234)
tsneOut <- Rtsne::Rtsne(t(scRNA),dims=2,
                        PCA=T,
                        perplexity=100,
                        verbose=F,
                        max_iter=500,
                        check_duplicates=F)
tsne <- data.frame(tSNE1=tsneOut$Y[,1],
                   tSNE2=tsneOut$Y[,2],
                   CIBERSORTx = results[,2],
                   Response=factor(response,levels = c("NR","R")))

myPalette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
p1 <- ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=Response))+
  theme_classic()+
  scale_colour_manual(values =  c("#E7B800","#2E9FDF"))+
  scale_fill_manual(values = c("#E7B800","#2E9FDF"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.background = element_rect(fill = "white", size = 1, colour = "white"))

p2 <- ggplot(tsne,aes(tSNE1,tSNE2,color=CIBERSORTx))+
  geom_point()+
  theme_classic()+
  scale_color_lancet()+
  ggplot2::scale_colour_gradientn(name = "CIBERSORTx\nscores",
                                  colours = myPalette(100), 
                                  guide = ggplot2::guide_colourbar(ticks.colour = "black",
                                                                   ticks.linewidth = 1, 
                                                                   frame.colour = "black"),
                                  breaks = seq(0, 1, 0.2), 
                                  labels = c(0.0 ,0.2, 0.4, 0.6, 0.8, 1.0))+
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        legend.position = "right", 
        legend.background = element_rect(fill = "white", size = 1, colour = "white"))


Fig.1A <- plot_grid(p2,p1,nrow = 2, align="v")
Fig.1A


boxdata <- data.frame(Cell = results$Mixture,TregSocre=results$Treg)
boxdata <- left_join(boxdata,phenotype[,c('Cell','Response')],by='Cell')
boxdata$Response <- factor(boxdata$Response,levels = c("NR","R"))

Fig.1B <- ggplot(boxdata, aes(x = Response, y = TregSocre, fill = Response)) +
  geom_flat_violin(aes(fill = Response),position = position_nudge(x = 0.1, y = 0),  trim = TRUE, alpha = .5, colour = NA)+
  geom_point(aes(x = .55, y = TregSocre, colour = Response),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(x = Response, y = TregSocre, fill = Response),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#E7B800","#2E9FDF"))+
  scale_fill_manual(values = c("#E7B800","#2E9FDF"))+
  theme_classic()+
  geom_signif(comparisons = list(c("R","NR")),map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=13))+
  # ylim(0,1)+
  ylab('CIBERSORTx scores')


Fig.1B



####====Fig. 1C and 1D===####

cohort <- "BCC_GSE123813_aPD1"
load(paste0("data/",cohort,'/scRNA.Rdata')) 
load(paste0("data/",cohort,'/phenotype.Rdata'))
load(paste0('data/CIBERSORTx/CIBERSORTx_result/BCC_GSE123813_aPD1_CIBERSORTx.Rdata'))


cell <- phenotype[phenotype$Celltype_major == 'Treg' ,'Cell']
response <- phenotype[phenotype$Celltype_major == 'Treg' ,'Response']
names(response) <- phenotype[phenotype$Celltype_major == 'Treg' ,'Cell']
scRNA <- scRNA[,colnames(scRNA) %in% cell]  
scRNA <- as.data.frame(scRNA)

rownames(CIBERSORTx)<-CIBERSORTx$Mixture
results <- CIBERSORTx[rownames(CIBERSORTx)%in% cell,c("Mixture","Treg")]

set.seed(1234)
tsneOut <- Rtsne::Rtsne(t(scRNA),dims=2,
                        PCA=T,
                        perplexity=100,
                        verbose=F,
                        max_iter=500,
                        check_duplicates=F)
tsne <- data.frame(tSNE1=tsneOut$Y[,1],
                   tSNE2=tsneOut$Y[,2],
                   CIBERSORTx = results[,2],
                   Response=factor(response,levels = c("NR","R")))

myPalette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
p3 <- ggplot(tsne,aes(tSNE1,tSNE2))+
  geom_point(aes(color=Response))+
  theme_classic()+
  scale_colour_manual(values =  c("#E7B800","#2E9FDF"))+
  scale_fill_manual(values = c("#E7B800","#2E9FDF"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.background = element_rect(fill = "white", size = 1, colour = "white"))

p4 <- ggplot(tsne,aes(tSNE1,tSNE2,color=CIBERSORTx))+
  geom_point()+
  theme_classic()+
  scale_color_lancet()+
  ggplot2::scale_colour_gradientn(name = "CIBERSORTx\nscores",
                                  colours = myPalette(100), 
                                  guide = ggplot2::guide_colourbar(ticks.colour = "black",
                                                                   ticks.linewidth = 1, 
                                                                   frame.colour = "black"),
                                  breaks = seq(0, 1, 0.2), 
                                  labels = c(0.0 ,0.2, 0.4, 0.6, 0.8, 1.0))+
  theme(plot.margin = unit(rep(1.5,4),"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        legend.position = "right", 
        legend.background = element_rect(fill = "white", size = 1, colour = "white"))

Fig.1C <- plot_grid(p4,p3,nrow = 2, align="v")
Fig.1C


boxdata <- data.frame(Cell = results$Mixture,TregSocre=results$Treg)
boxdata <- left_join(boxdata,phenotype[,c('Cell','Response')],by='Cell')
boxdata$Response <- factor(boxdata$Response,levels = c("NR","R"))

Fig.1D <- ggplot(boxdata, aes(x = Response, y = TregSocre, fill = Response)) +
  geom_flat_violin(aes(fill = Response),position = position_nudge(x = 0.1, y = 0),  trim = TRUE, alpha = .5, colour = NA)+
  geom_point(aes(x = .55, y = TregSocre, colour = Response),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(x = Response, y = TregSocre, fill = Response),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#E7B800","#2E9FDF"))+
  scale_fill_manual(values = c("#E7B800","#2E9FDF"))+
  theme_classic()+
  geom_signif(comparisons = list(c("R","NR")),map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=13))+
  # ylim(0,1)+
  ylab('CIBERSORTx scores')

Fig.1D
