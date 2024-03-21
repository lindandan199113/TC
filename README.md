<<<<<<< HEAD
# TC
=======
---
output:
  html_document: default
  pdf_document: default
  word_document: default
---

### Environment setup

#### 1. open "TC.Rproj" 
#### 2. run the following codes to setup R environment

```R
# install external packages
renv::restore(exclude=c('ggradar','hdf5r','msigdb.data','multiGSEA','reactome.db','ReactomePA'))
# install local R packages 
install.packages('renv/local/yulab.utils-master.zip', repos = NULL, type="source")
library(yulab.utils)#Load yulab.utils package
install_zip("renv/local/ggradar-master.zip")#Use the install_zip function in the yulab.utils package to load the ggradar package
install.packages('renv/local/hdf5r_1.3.5.tar.gz', repos = NULL, type="source")
install.packages('renv/local/lianos-msigdb.data-68db424.tar.gz', repos = NULL, type="source")
install.packages('renv/local/lianos-multiGSEA-prebioc-61-g85719e2.tar.gz', repos = NULL, type="source")
install.packages('renv/local/reactome.db_1.70.0.tar.gz', repos = NULL, type="source")
install.packages('renv/local/ReactomePA_1.32.0.tar.gz', repos = NULL, type="source")
# check R environment
renv::status()
```
### Folder structures 

```bash
│ GM.Rproj
├─Code # R scripts
│      bnplot.R 
│      corplot.R 
│      DevOfSig.R # For development of Treg.Sig and Treg.SKCM.Sig
│      Fig1.R # "Fig1.R", "Fig2.R"..."FigS2.R" are R scripts used to generate figures presented in original study
│      Fig2.R
│      Fig3.R
│      Fig4.R
│      Fig5.R
│      Fig6.R
│      FigS1.R
│      FigS2.R
│      Model_Treg.Sig.R # For model training, validation, and independent testing of Treg.Sig
│      R_rainclouds.R
│    
├─data # data used to generate results
├─renv # version information
└─results # intermediate variables and results 
```

### Key files explained

#### 1) "DevOfSig.R"

1. structure of DevOfSig.R, code is divided by comments into the following parts.

```R
####===1. load dir===####
####===2. Genes of CIBERSORTx score for each scRNA-Seq cohort===####
####===3. Differential expression genes of Treg cells===####
####===4. get Treg.Sig===####
```

2. Run the codes. Treg.Sig can be obtained in Part.4 respectively.

3. Treg.Sig is derived from **pan-cancer scRNA-Seq analysis**.

4. Development of Treg.Sig: 
- **Gx**: we performed Spearman correlation analysis between gene expression level and CIBERSORTx scores for malignant cells among pan-cancer scRNA datasets. Genes that were positively correlated with CIBERSORTx scores (Spearman R > 0 & FDR < 1e-05) were regarded as Gx. 
- **Gy**: Genes that were differentially up-regulated in Treg cells were regarded as Gy. 
- **Gn**: To obtain up-regulated specific genes that were positively associated with Treg cells, Gx and Gy were intersected to give rise to Gn for each dataset. **For example**, G1 consisted of genes derived from the intersection of Gx and Gy in the first scRNA-Seq dataset. Geometric mean of spearman R was calculated for each gene across G1-G32. 
- **Treg.Sig**: genes with geometric mean of spearman R > 0.5 (moderate to strong correlation) were pooled as Treg.Sig.



####  2) "Model_Treg.Sig.R" 

1. structure of Model Prediction.R

```R
####===1. loading data===####
####===2. remove batch effect===####
####===3. Model Training and Validation with Stem.Sig===####
####===4. Model Testing with Stem.Sig===####
####===5. Compare Stem.Sig with other Signatures in testing datasets===####
```

2. run the codes, model training process has been implemented in the function _"CompareModel"_ in Part.3. 

- Details of _"CompareModel"_

```R
#parameters: 1)training set; 2)validation set; 3)algorithms; 4)gene signature

res <- CompareModel(training = training set, validation = validation set, method = algorithms, sig = gene signature)

#return: 1) res[['model']]: models derived from each specific algorithms after training; 2)res[['auc']]: validation auc for each algorithm.

#We pick the algorithm with highest validation auc as the best Stem.Sig model ("cancerclass" algorithm for this study)

```

- Prediction results of other public signatures could be found in part.5. Corresponding reference and algorithms for each public signature were documented as well.


### sessionInfo()

```R
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8
[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
  [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
  [1] multiGSEA_1.1.99       recipes_0.1.17         plotrix_3.8-2          lmerTest_3.1-3
[5] lme4_1.1-27.1          Matrix_1.4-0           mixtools_1.2.0         rms_6.2-0
[9] SparseM_1.81           Hmisc_4.6-0            Formula_1.2-4          ROCR_1.0-11
[13] ppcor_1.1              MASS_7.3-54            plyr_1.8.6             matrixStats_0.61.0
[17] doParallel_1.0.16      iterators_1.0.13       foreach_1.5.1          pROC_1.18.0
[21] e1071_1.7-9            sva_3.42.0             BiocParallel_1.28.3    genefilter_1.76.0
[25] mgcv_1.8-38            nlme_3.1-153           gridExtra_2.3          future_1.23.0
[29] OptimalCutpoints_1.1-5 cutpointr_1.1.1        scales_1.1.1           ggradar_0.2
[33] ggrepel_0.9.1          corrplot_0.92          MCPcounter_1.1.0       curl_4.3.2
[37] pheatmap_1.0.12        data.table_1.14.2      GSEABase_1.56.0        graph_1.72.0
[41] annotate_1.72.0        XML_3.99-0.8           GSVA_1.42.0            ReactomePA_1.32.0
[45] org.Hs.eg.db_3.14.0    AnnotationDbi_1.56.2   IRanges_2.28.0         S4Vectors_0.32.3
[49] clusterProfiler_4.2.1  gridBase_0.4-7         grImport2_0.2-0        stringr_1.4.0
[53] png_0.1-7              CytoTRACE_0.3.3        SeuratObject_4.0.4     Seurat_4.0.6
[57] cancerclass_1.38.0     binom_1.1-1            Biobase_2.54.0         BiocGenerics_0.40.0
[61] ggsci_2.9              ComplexHeatmap_2.10.0  circlize_0.4.13        RColorBrewer_1.1-2
[65] reshape2_1.4.4         survival_3.2-13        survminer_0.4.9        ggpubr_0.4.0
[69] dplyr_1.0.7            caret_6.0-90           lattice_0.20-45        ggplot2_3.3.5
[73] ROCit_2.1.1

loaded via a namespace (and not attached):
  [1] graphlayouts_0.7.2          pbapply_1.5-0               haven_2.4.3
[4] vctrs_0.3.8                 graphite_1.40.0             blob_1.2.2
[7] prodlim_2019.11.13          spatstat.data_2.1-2         later_1.3.0
[10] nloptr_1.2.2.3              DBI_1.1.2                   SingleCellExperiment_1.16.0
[13] rappdirs_0.3.3              uwot_0.1.11                 jpeg_0.1-9
[16] zlibbioc_1.40.0             MatrixModels_0.5-0          htmlwidgets_1.5.4
[19] mvtnorm_1.1-3               GlobalOptions_0.1.2         hdf5r_1.3.5
[22] leiden_0.3.9                irlba_2.3.5                 DEoptimR_1.0-9
[25] tidygraph_1.2.0             Rcpp_1.0.7                  KernSmooth_2.23-20
[28] promises_1.2.0.1            DelayedArray_0.20.0         limma_3.50.0
[31] fastmatch_1.1-3             digest_0.6.29               polspline_1.1.19
[34] ncdf4_1.19                  sctransform_0.3.2           scatterpie_0.1.7
[37] cowplot_1.1.1               DOSE_3.20.1                 here_1.0.1
[40] ggraph_2.0.5                pkgconfig_2.0.3             GO.db_3.14.0
[43] DelayedMatrixStats_1.16.0   gower_0.2.2                 nnls_1.4
[46] minqa_1.2.4                 reticulate_1.22             SummarizedExperiment_1.24.0
[49] GetoptLong_1.0.5            xfun_0.29                   zoo_1.8-9
[52] tidyselect_1.1.1            purrr_0.3.4                 kernlab_0.9-29
[55] ica_1.0-2                   labelled_2.9.0              pcaPP_1.9-74
[58] viridisLite_0.4.0           rlang_0.4.12                glue_1.6.0
[61] MatrixGenerics_1.6.0        lava_1.6.10                 ggsignif_0.6.3
[64] labeling_0.4.2              httpuv_1.6.4                class_7.3-19
[67] TH.data_1.1-0               reactome.db_1.70.0          DO.db_2.9
[70] jsonlite_1.7.2              XVector_0.34.0              bit_4.0.4
[73] mime_0.12                   stringi_1.7.6               spatstat.sparse_2.1-0
[76] scattermore_0.7             yulab.utils_0.0.4           bitops_1.0-7
[79] rhdf5filters_1.6.0          RSQLite_2.2.9               randomForest_4.6-14
[82] tidyr_1.1.4                 rstudioapi_0.13             qvalue_2.26.0
[85] locfit_1.5-9.4              listenv_0.8.0               miniUI_0.1.1.1
[88] gridGraphics_0.5-1          survMisc_0.5.5              segmented_1.3-4
[91] lifecycle_1.0.1             questionr_0.7.5             timeDate_3043.102
[94] kknn_1.3.1                  munsell_0.5.0               caTools_1.18.2
[97] codetools_0.2-18            GenomeInfoDb_1.30.0         lmtest_0.9-39
[100] htmlTable_2.3.0             xtable_1.8-4                BiocManager_1.30.16
[103] abind_1.4-5                 farver_2.1.0                parallelly_1.30.0
[106] km.ci_0.5-2                 RANN_2.6.1                  aplot_0.1.1
[109] klaR_0.6-15                 ggtree_3.2.1                GenomicRanges_1.46.1
[112] RcppAnnoy_0.0.19            goftest_1.2-3               patchwork_1.1.1
[115] tibble_3.1.6                cluster_2.1.2               future.apply_1.8.1
[118] tidytree_0.3.6              ellipsis_0.3.2              lubridate_1.8.0
[121] ggridges_0.5.3              igraph_1.2.10               fgsea_1.20.0
[124] spatstat.utils_2.3-0        htmltools_0.5.2             utf8_1.2.2
[127] plotly_4.10.0               ModelMetrics_1.2.2.2        foreign_0.8-81
[130] withr_2.4.3                 fitdistrplus_1.1-6          bit64_4.0.5
[133] multcomp_1.4-17             robustbase_0.93-9           Biostrings_2.62.0
[136] spatstat.core_2.3-2         combinat_0.0-8              fastAdaboost_1.0.0
[139] GOSemSim_2.20.0             rsvd_1.0.5                  ScaledMatrix_1.2.0
[142] memoise_2.0.1               forcats_0.5.1               HiClimR_2.2.0
[145] fansi_0.5.0                 highr_0.9                   tensor_1.5
[148] conquer_1.2.1               edgeR_3.36.0                checkmate_2.0.0
[151] renv_0.14.0                 cachem_1.0.6                deldir_1.0-6
[154] tensorA_0.36.2              rjson_0.2.20                rstatix_0.7.0
[157] clue_0.3-60                 rprojroot_2.0.2             tools_4.1.1
[160] sandwich_3.0-1              magrittr_2.0.1              RCurl_1.98-1.5
[163] proxy_0.4-26                car_3.0-12                  ape_5.6
[166] bayesm_3.1-4                ggplotify_0.1.0             httr_1.4.2
[169] boot_1.3-28                 globals_0.14.0              R6_2.5.1
[172] Rhdf5lib_1.16.0             nnet_7.3-16                 KEGGREST_1.34.0
[175] treeio_1.18.1               shape_1.4.6                 beachmat_2.10.0
[178] HDF5Array_1.22.1            BiocSingular_1.10.0         rhdf5_2.38.0
[181] splines_4.1.1               carData_3.0-4               ggfun_0.0.4
[184] colorspace_2.0-2            generics_0.1.1              base64enc_0.1-3
[187] compositions_2.0-2          pillar_1.6.4                tweenr_1.0.2
[190] ccaPP_0.3.3                 GenomeInfoDbData_1.2.7      egg_0.4.5
[193] gtable_0.3.0                knitr_1.37                  latticeExtra_0.6-29
[196] shadowtext_0.1.0            fastmap_1.1.0               quantreg_5.86
[199] msigdb.data_0.99.5          broom_0.7.10                backports_1.4.1
[202] ipred_0.9-12                enrichplot_1.14.1           hms_1.1.1
[205] ggforce_0.3.3               Rtsne_0.15                  shiny_1.7.1
[208] KMsurv_0.1-5                polyclip_1.10-0             numDeriv_2016.8-1.1
[211] lazyeval_0.2.2              crayon_1.4.2                downloader_0.4
[214] sparseMatrixStats_1.6.0     viridis_0.6.2               rpart_4.1-15
[217] compiler_4.1.1              spatstat.geom_2.3-1
```







>>>>>>> 75d305a (提交文件)
