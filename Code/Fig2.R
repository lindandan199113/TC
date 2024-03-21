####===Figure 2===####
rm(list=ls())
gc()


####===Fig. 2A===####
library(circlize)
library(RColorBrewer)
library(png)
library(graphics)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(grImport2)
library(gridBase)
library(ggplot2)

dir <- data.table::fread('data/dir/dir_Treg1.csv')
circos.clear() 
Set3 <- brewer.pal(12,"Set3") 
Set2 <- brewer.pal(8,"Set2")
col <- c(Set3,Set2)

sectors = dir$dataset
names(sectors) <- sectors

dataset <- data.frame(dataset = sectors)
rownames(dataset) <- sectors

col_cancer <- c(col[1], col[2], rep(col[3],2),rep(col[4],3),col[5],rep(col[6],2),
                col[7],rep(col[8],2),rep(col[9],3),col[10],rep(col[11],2),
                rep(col[12],7),
                col[13],col[14],col[15],rep(col[16],2),col[19])
names(col_cancer) <- sectors

image = 'data/png/venn.png'
image = as.raster(readPNG(image))


circos.clear()
windows()
circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.01, 0.01)) 
circos.heatmap(dataset, split = dataset$dataset, col = col_cancer, 
               track.height = 0.02,rownames.side = 'outside')

circos.track(ylim = c(0, 1),track.height=0.15,bg.border = "#EEEEEE", panel.fun = function(x, y) {
  circos.raster(image, CELL_META$xcenter, CELL_META$ycenter, 
                width = "0.7cm",
                facing = "inside")
  
})

circos.track(ylim = c(0, 1), sectors = sectors,
             bg.col = col_cancer, bg.border = "#EEEEEE" , track.height = 0.25)

circos.trackText(x = rep(0.5, 32), y = rep(0.5, 32),
                 labels = paste0(rep('G',32), 1:32),
                 cex = 0.5, sectors = sectors, col = "white", font = 2, facing = "clockwise",
                 niceFacing=T)

draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360,
            rou1 = 0.25, col = "#5ea3c4", border = "#EEEEEE")

text(0,0,"Treg.sig",col = 'white',cex = 2,font = 2)

cancer <- unique(col_cancer)
names(cancer) <- unique(dir$cancer)
lgd_cancer = Legend(title = "Cancer", at = names(cancer),
                    legend_gp = gpar(fill = cancer))
draw(lgd_cancer, x = unit(1.2, "snpc"), just = "left")




####===Fig. 2B===####
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

load('data/sig/Treg.Sig.Rdata')

sig <- Treg.Sig

sig_ID <- bitr(sig,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")

eReac <- enrichPathway(gene = sig_ID$ENTREZID,
                       organism = 'human',
                       pvalueCutoff = 0.05)

bardata <- eReac@result[1:16,]
bardata <- bardata[order(bardata$qvalue),]
bardata <- bardata[order(bardata$Description),]
bardata$Description <- factor(bardata$Description,levels=rev(bardata$Description))
barplot <- ggplot(bardata,aes(y=Description,x=Count))+
  geom_bar(stat = "identity",aes(fill=qvalue))+
  scale_fill_gradientn(colours =  c("#4596B4", "#BAC2CC"))+theme_bw()

eReacx <- setReadable(eReac, 'org.Hs.eg.db', 'ENTREZID')

dropAsis <- function(x){
  cls <- class(x)
  structure(x, class = setdiff(cls, "AsIs"))
}
rescale.AsIs <- function(x, ...){
  
  dropAsis <- function(x){
    cls <- class(x)
    structure(x, class = setdiff(cls, "AsIs"))
  }
  
  scales:::rescale(dropAsis(x), ...)
}

cnetplot <-  cnetplot(eReacx, circular = TRUE, colorEdge = TRUE,
                      categorySize="pvalue",showCategory = 30,layout = 'kk')


barplot 
cnetplot

