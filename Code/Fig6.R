####===Figure 6===####
rm (list=ls())
gc()


####==loading data==####
library(dplyr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggplot2)
library(ggradar)
load('data/CRISPR/crispr_new.Rdata')  
crispr <- crispr_new3
rownames(crispr)<-crispr$Gene
crispr<-crispr[,-18]

load('data/sig/Treg.Sig.Rdata') 
source('data/sig/sig_others.R')


####===Figure 6A===####

df <- crispr
df$Mean <- apply(df,1,function(z){mean(z,na.rm=T)})
df <- df[order(df$Mean),]  
df$Rank <- seq(1,length(df$Mean))
df <- df[c(1:8, (nrow(df)-7):nrow(df)),]
range(df,na.rm=T)

ann_colors =list(cohort =c(  'Freeman 2019--NK--Melanoma--B16' = "#6BA38B",
                             'Freeman 2019--OT1--Melanoma--B16' = "#BA7D7D",
                             'Kearney 2018--NK--Colon--MC38' = "#9C8FB3",
                             'Kearney 2018--OT1+IgG--Colon--MC38' = "#C47D62",
                             'Kearney 2018--OT1+PD1--Colon--MC38' = "#6CA2B0",
                             'Lawson 2020--Mid_CTL--Breast--EMT6HA' =  '#8EBFB2',
                             'Lawson 2020--Mid_CTL--Breast--X4T1HA' = '#7294AD',
                             'Lawson 2020--Mid_CTL--Colon--CT26HA' = '#CC9BA3',
                             'Lawson 2020--Mid_CTL--Colon--MC38OVA' = '#D4C486',
                             'Lawson 2020--Mid_CTL--Melanoma--B16OVA' = '#92BDD1',
                             'Lawson 2020--Mid_CTL--Renal--RencaHA' = '#A1E3D5',
                             'Manguso 2017--GVAX--Melanoma--B16' = '#E6A8A8',
                             'Manguso 2017--GVAX+PD1--Melanoma--B16' = '#C1ABDB',
                             'Pan 2018--OT1--Melanoma--B16' = '#E3AE96',
                             'Pan 2018--Pmel-1_T--Melanoma--B16' = '#A3C2CC',
                             'Patel 2017--NY-ESO-1_MART-1_T--Melanoma--Mel624' = '#74AD9F',
                             'Vredevoogd 2019--MART-1_T--Melanoma--D10 IFNGR1-' = '#839DB5',
                             'Gu 2023--EPI+OT1--Pancreas--KPC3OVA' = '#CCAFB3',
                             'Gu 2023--MES+OT1--Pancreas--KPC3OVA' = '#DED099',
                             'Griffin 2021--GVAX+PD1--Melanoma--B16' = '#B4CCD9',
                             'Griffin 2021--PD1+CTLA4--Lung--LLC' = '#C0DBD0',
                             'Dubrot 2022--GVAX+PD1--Melanoma--B16' = '#F2A9A9',
                             'Dubrot 2022--PD1+CTLA4--Colon--CT26' = '#D9D1E6',
                             'Dubrot 2022--PD1+CTLA4--Pancreas--KPC' = '#F5B9A1',
                             'Dubrot 2022--PD1+CTLA4--Lung--LLC' = '#CEE6F2',
                             'Dubrot 2022--PD1--Colon--MC38' = '#DCF5F0',
                             'Dubrot 2022--PD1+CTLA4--Pancreas--Panc02' = '#B3D5EB',
                             'Dubrot 2022--PD1+CTLA4--Renal--Renca' = '#F5D2D0',
                             'Dubrot 2022--PD1+CTLA4--Melanoma--YUMMER' = '#E0D8F0',
                             Mean = '#ffffff',
                             Rank = '#ffffff'),
                 cancer = c(Melanoma="#4778a7",Colon="#8ac1e3",Breast="#f08c37",
                            Renal="#ecca61",Pancreas="#ac96bc",Lung = '#822773',
                            Mean = '#ffffff',Rank = '#ffffff'))
column_ha = HeatmapAnnotation(cohort = colnames(df),
                              cancer=c('Melanoma','Melanoma','Colon','Colon',
                                       'Colon','Breast','Breast','Colon','Colon',
                                       'Melanoma','Renal', rep('Melanoma',6),'Pancreas',
                                       'Pancreas','Melanoma','Lung','Melanoma',
                                       'Colon','Pancreas','Lung','Colon','Pancreas',
                                       'Renal','Melanoma',"Mean","Rank"),
                              col = ann_colors,
                              gp = gpar(col = "black", lwd = 2)
)

figure.6A<-Heatmap(as.matrix(df), name = "z scores",
                   column_title = NULL,
                   row_title = NULL,
                   cluster_rows = F,
                   cluster_columns = F,
                   col = colorRamp2(c(-2, 0, 2), c("#5FB54A", "grey", "#d52b30")),
                   show_row_names = T, show_column_names = F,
                   row_names_gp = gpar(fontsize = 10),
                   rect_gp = gpar(col = "black", lwd = 2),
                   width = ncol(df)*unit(5, "mm"), 
                   height = nrow(df)*unit(5, "mm"),
                   na_col = 'white',
                   column_names_side = c('top'),
                   row_split = c(rep('a',8),rep('b',8)),
                   column_split = c(rep('a',29),'b','c'),
                   top_annotation = column_ha
)

figure.6A



###===ranking of genes ===###
rank <- apply(crispr,1,function(z){ mean(z,na.rm=T)})
rank <- data.frame(meanZ=rank,row.names = rownames(crispr))
rank$genes <- rownames(crispr)
rank <- rank[order(rank$meanZ),]
rank$order <- order(rank$meanZ)


####==Fig.6B==####

#signatures with immune resistant genes
num <- c('Treg.Sig',
         'TcellExc.Sig',
         'IMS.Sig',
         'LRRC15.CAF.Sig',
         'CRMA.Sig',
         'ImmuneCells.Sig')

compare_crs <- data.frame(row.names = num)
row <- data.frame(row.names =  rank$genes)
ft <- data.frame()

for(i in c(0.03,0.06,0.09,0.12,0.15,0.18)){
  for(j in num){
    if (j == 'TcellExc.Sig'){
      load("data/sig/FINAL_Tcell_Exclusion.sig.RData")
      sig <- exc.sig$exc.up
    } else if (j == 'IMS.Sig'){
      sig <- IMS.Sig
    } else if (j == 'LRRC15.CAF.Sig'){
      sig <- LRRC15.CAF.Sig
    } else if (j == "CRMA.Sig"){
      sig <- CRMA.Sig
    } else if (j == 'Treg.Sig') {
      sig <- Treg.Sig
    } else if ( j == "ImmuneCells.Sig"){
      sig <- readRDS('data/sig/ImSig.rds')
    } 
    
    row[,j] <- ifelse(rank$genes %in% sig,'p','n')
    
    r <- rank[rank$genes %in% sig,]
    
    
    compare_crs[j,paste(i)] <- sum(r$order<=nrow(rank)*i) / nrow(r)
    
    ft[paste(i,'1'),j] <- sum(r$order<=nrow(rank)*i)
    ft[paste(i,'0'),j] <- nrow(r)-sum(r$order<=nrow(rank)*i)
    
    print(paste(j , i*100,'%' , sum(r$order<=nrow(rank)*i),nrow(r)))
    
  }
}

compare_crs <- as.data.frame(t(compare_crs))

compare_crs$sig <- as.character(paste('top',as.numeric(rownames(compare_crs))*100,'%','top-ranked genes'))


compare_crs$sig <- factor(compare_crs$sig,levels = compare_crs$sig)

compare_crs <- compare_crs[,c(ncol(compare_crs),1:(ncol(compare_crs)-1))]

compare_crs


figure.6B <- ggradar2(compare_crs,
                      all.radar=c(0,0.036,0.071,0.11,0.18,0.23),
                      group.colours = c('#b27573','#bd9249','#669982',
                                        '#6e7a5b','#3b545c','#16332B'),
                      legend.position = "left")


figure.6B


##Fig. 6C
df_rank <- crispr

r <- rank[rank$genes %in% Treg.Sig,]
df_rank <- df_rank[rownames(r)[r$order<nrow(rank)*0.18],] #top 18% genes
df_rank$`mean Z Scores` <- rank[rownames(df_rank),'meanZ']
df_rank <- t(df_rank)

figure.6C <- Heatmap(df_rank, name = "CRISPR immune Scores",
                     row_title = NULL,
                     cluster_rows = F,
                     cluster_columns = F,
                     col = colorRamp2(c(-2, 0, 2), c("#5FB54A", "grey", "#d52b30")),
                     show_row_names = T, show_column_names = T,
                     width = ncol(df_rank)*unit(5, "mm"), 
                     height = nrow(df_rank)*unit(5, "mm"),
                     rect_gp = gpar(col = "black", lwd = 2),
                     na_col = 'white',
                     column_names_side = c('top'),
                     row_split = c(rep('a',29),'b')
)

figure.6C
