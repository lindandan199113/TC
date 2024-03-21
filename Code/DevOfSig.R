
####===Development of Stem.Sig===####
library(Seurat)
library(dplyr)
library(stringr)
library(cancerclass)
library(cutpointr)
library(OptimalCutpoints)
library(reshape2)
library(future)

rm(list=ls())
gc()


##################################################################################
####=================================1. load dir==============================####
##################################################################################

cells <- 'Treg'
load(paste0('data/dir/dir_',cells,'.Rdata')) ## loading dir for scRNA-Seq cohorts with both available malignant and stromal/immune cells data. 
dir <- dir_Treg 



##################################################################################
####===============2. Genes of CIBERSORTx score for each scRNA-Seq cohort=================####
##################################################################################

#Plsease Download 32 scRNA-Seq datasets from http://tisch.comp-genomics.org/
##For demonstration, we provided scRNA-data and meta data of BLCA_GSE149652 dataset. 

##Results of CIBERSORTx analysis for each scRNA-Seq dataset were provided in 'data/CIBERSORTx/CIBERSORTx_result/'


for(n in dir$cohort[2]){ #demonstration using dataset BLCA_GSE149652
  scRNA <- Read10X_h5(paste0('data/tisch2_h5/',n,'_expression.h5'))
  phenotype <- read.delim(paste0('data/scRNAmeta/',n,'_CellMetainfo_table.tsv'))
  load(paste0('data/CIBERSORTx/CIBERSORTx_result/',n,'_CIBERSORTx.Rdata'))
  phenotype <- as.data.frame(phenotype)
  cellType <- phenotype[phenotype$Celltype..major.lineage. == 'Treg' ,]
  exprMatrix <- scRNA[ ,colnames(scRNA) %in% cellType$Cell]
  exprMatrix <- as.data.frame(exprMatrix) 
  
  rownames(CIBERSORTx)<-CIBERSORTx$Mixture
  CIBERSORTx_score <- CIBERSORTx[rownames(CIBERSORTx) %in% cellType$Cell,c("Mixture","Treg")]#筛出Treg细胞的CIBERSORTx评分
  #CIBERSORTx_score <- select(CIBERSORTx_score,-1)
  colnames(CIBERSORTx_score) <- c('cell','CIBERSORTx_score')

  #########future multicore##########
  future::plan("multisession", workers = 40)
  ##########################################
   
  ###===Calculate Spearman correlation of Treg CIBERSORTx score genes===####
  Genes <- data.frame(gene=rownames(exprMatrix),coef=NA,p=NA)
  
  for(i in 1:nrow(exprMatrix)){
    cor <- cor.test(as.numeric(exprMatrix[i,]),CIBERSORTx_score$CIBERSORTx_score,method = 'pearson')
    Genes$coef[i] <- cor$estimate
    Genes$p[i] <- cor$p.value
  }
  
  Genes$p.adjust <- p.adjust(Genes$p,method = 'BH')
  
  save(Genes,file=paste0('results/Genes/' , n , '_Genes.Rdata')) 
}




##################################################################################
####==========3. Differential expression genes of Treg cells=============####
##################################################################################

## results were provided in 'results/DE/'

for(n in dir$cohort[2]){ #demonstration using dataset BLCA_GSE149652
  
  scRNA <- Seurat::Read10X_h5(paste0('data/tisch2_h5/',n,'_expression.h5'))
  
  phenotype <- read.delim(paste0('data/scRNAmeta/',n,'_CellMetainfo_table.tsv'))
  phenotype <- as.data.frame(phenotype)
  phenotype <- phenotype[,c(1,4:6)]
  
  table(phenotype$Celltype..major.lineage.)
  
  phenotype[!phenotype$Celltype..major.lineage. == cells, 3] <- "control"
  meta <- phenotype
  
  sce <- CreateSeuratObject(counts = scRNA, 
                            project = "test")
  
  sce@meta.data$group <- meta$Celltype..major.lineage.
  
  #########future multicore##########
  future::plan("multisession", workers = 40)
  ##########################################
  
  DE <- FindMarkers(sce,ident.1 = cells, group.by = 'group',logfc.threshold = 0.25,min.pct = 0.1,base=exp(1)) 
  DE <- DE[DE$p_val_adj<1e-05,]
  
  save(DE,file = paste0('results/DE/' , n , '_DE.Rdata')) 
}


##################################################################################
####===============================4. get Treg.Sig============================####
##################################################################################


getSig <- function(dir){
  
  ls_Gn <- list()
  allGenes <- c()  
  
  for(n in dir$cohort){ 
    load(file=paste0('results/Genes/' , n , '_Genes.Rdata'))
    load(paste0('results/DE/' , n , '_DE.Rdata'))
    
    Gx <- Genes[Genes$coef>0 & Genes$p.adjust < 1e-05,] 
    
    Gy <- rownames(DE[DE$avg_logFC>0,]) ##
    Gy <- rownames(DE) ##
    Gy <- Gy[!grepl("^RP[SL]", Gy,ignore.case = F)] ##(ribosome protein free)
    
    Gn <- Gx[Gx$gene %in% Gy,]
    
    allGenes <- c(allGenes,Gn$gene)
    ls_Gn[[n]] <- Gn
  }
  
  
  allGenes <- unique(allGenes)
  allGenesDf <- data.frame(gene = allGenes)
  
  for(i in dir$cohort){
    allGenesDf <- left_join(allGenesDf,ls_Gn[[i]][,c('gene','coef')],by='gene')
  }
  
  allGenesDf <- allGenesDf[!is.na(allGenesDf$gene),] 
  rownames(allGenesDf) <- allGenesDf$gene
  allGenesDf <- allGenesDf[,-1]
  colnames(allGenesDf) <- dir$cohort
  genelist <- allGenesDf
  genelist$all_gmean <- compositions::geometricmeanRow(genelist[,1:length(dir$cohort)]) ##spearmanR geometric mean
  sig <- genelist[genelist$all_gmean>0.5,] ##filter genes with spearmanR geometric mean > 0.5
  sig <- sig[order(sig$all_gmean,decreasing = T),]
  sig <- rownames(sig)
  return(sig)
  
}


Treg.Sig <- getSig(dir)

#save(Treg.Sig,file='data/sig/Treg.Sig.Rdata')

