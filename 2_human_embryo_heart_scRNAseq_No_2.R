# GSE181346

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(1)

dat <- readRDS("260228_Heart_ATAC_F6_8_12_19_merge_No0.rds")
# dat
# An object of class Seurat 
# 68120 features across 17617 samples within 1 assay 
# Active assay: RNA (68120 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, umap, harmony

dat$time <- factor(
  x = dat$time,
  levels = c("PCW6", "PCW8", "PCW12", "PCW19")
)

Figure1A <- DimPlot(
  object = dat,
  group.by = 'time',
  label = FALSE,
  label.size = 6) + ggtitle('Sample') + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure1A # W700 H670

Figure1B <- DimPlot(
  object = dat,
  group.by = 'Cluster_names',
  label = TRUE,
  repel = TRUE,
  label.size = 6) + ggtitle('Cell type') + 
  theme(legend.position = "none")  + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure1B # W700 H670

# markers_OLD <- c("TTN", "MYH7", "ANKRD1", "TNNT2", "NPPA", "MYL2", "DSP",
#              "PKP2", "HAND1", "MYOCD", "KCNH2", "TNNI3", "ALPK2",
#              "SCN5A", "MYL4", "MYL3", "GJA5", "TBX5", "MYL7", "MYBPC3",
#              "KCNJ5", "TCAP", "MYH6",
#              "PECAM1", "CDH5", "APLNR", "SOX17", "KDR", "DLL4", "FLT1",
#              "GATA2", "SELE", "EMCN", "MEOX2", "TEK", "NDNF", "JCAD",
#              "ROBO4", "MYDGF", "COL1A1", "MFAP4", "TWIST1","GDNF", "TAGLN",
#              "PDGFRB", "MYH11", "ACTA2",
#              "DCN", "LUM", "TBX18", "TCF21",
#              "WT1", "COL3A1", "VCAM1", "DDR2")

Figure1C <- FeaturePlot(
  object = dat,
  features = "CDH5",
  pt.size = 0.1
) + theme(
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 20),
  axis.text = element_text(size = 20)
)
Figure1C

FeaturePlot(dat, features = "CDH11")

markers <- c("MYH7", "MYH6", "NPPA", "HAND1",  "KCNE11", "MYL2", 
             "GJA5",  "TBX5", "MYL4", "MYL3",  "MYOCD",  "MYBPC3",
             "CDH5", "PECAM1", "APLNR", "SOX17", "KDR", "DLL4", "FLT1",  
             "ROBO4", "RGCC", "NRP2", "TGFB1", "FABP4", "NFATC1", "CDH11",　"KLF1", "TAL1", "VIM", 
             "COL3A1","COL1A1", "RSPO3", "MFAP4", "TAGLN","PDGFRB", "MYH11", "ACTA2",
             "DCN", "LUM", "TBX18", "TCF21", "WT1", "DDR2")

Figure1D <- DotPlot(dat, features = markers) + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("aCM", "vCM","Endo1", "lEC", "aEC", "Cap" ,"vEC", "NC", "OFT", "FB2", "CF",
               "SMC", "PC", "EPC")
  )
Figure1D # W1600 H670

notch <- c("NOTCH1", "NOTCH2","NOTCH3","NOTCH4", 
           "HES1", "HES2","HES3","HES4", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HEYL")

Figure1E <- DotPlot(dat, features = notch) + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("aCM", "vCM","Endo1", "lEC", "aEC", "Cap" ,"vEC", "NC", "OFT", "FB2", "CF",
               "SMC", "PC", "EPC")
  )
Figure1E # W900 H670

Figure1D_data <- Figure1D$data
head(Figure1D_data)
Figure1D_data_2 <- Figure1D_data[,2:5]
head(Figure1D_data_2)

Figure1E_data <- Figure1E$data
head(Figure1E_data)
Figure1E_data_2 <- Figure1E_data[,2:5]
head(Figure1E_data_2)

human_scRNAseq_avgexp_pctexp <- rbind(Figure1D_data_2, Figure1E_data_2)
head(human_scRNAseq_avgexp_pctexp)
write.csv(human_scRNAseq_avgexp_pctexp, "human_scRNAseq_avgexp_pctexp.csv")

cluster_data <- table(Idents(dat), dat$Sample)
head(cluster_data)
write.csv(cluster_data, "human_scRNAseq_cluster_data.csv")

Figure1F <- FeaturePlot(
  object = dat,
  features = "CDH5",
  pt.size = 0.1
) + theme(
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 20),
  axis.text = element_text(size = 20)
)
Figure1F # W700 H670

Figure1G <- FeaturePlot(
  object = dat,
  features = "NOTCH1",
  pt.size = 0.1
) + theme(
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 20),
  axis.text = element_text(size = 20)
)
Figure1G # W500 H530

Figure1H <- FeaturePlot(
  object = dat,
  features = "HES1",
  pt.size = 0.1
) + theme(
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 20),
  axis.text = element_text(size = 20)
)
Figure1H # W500 H530

Figure1I <- FeaturePlot(
  object = dat,
  features = "HES4",
  pt.size = 0.1
) + theme(
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 20),
  axis.text = element_text(size = 20)
)
Figure1I # W500 H530

Figure1J <- FeaturePlot(
  object = dat,
  features = "HEY1",
  pt.size = 0.1
) + theme(
  axis.title = element_text(size = 20),
  legend.text = element_text(size = 20),
  axis.text = element_text(size = 20)
)
Figure1J # W500 H530

###############################

cell_counts_by_time <- table(dat$time)
print(cell_counts_by_time)
# PCW12 PCW19  PCW6  PCW8 
# 11260  2364  3567   426 

cell_counts_by_cluster <- table(dat$Cluster_names)
print(cell_counts_by_cluster)
# aCM   aEC   Cap    CF Endo1   EPC   FB2   lEC    NC   OFT    PC   SMC   vCM 
# 499   315  1076  7049  1992   129   862    90   186   818   897   371  3249 
# vEC 
# 84 

levels(dat$time)
dat$time <- factor(x = dat$time,levels = c("PCW6", "PCW8", "PCW12", "PCW19"))

Figure2A1 <- FeaturePlot(
  object = dat,
  features = "NOTCH1",
  split.by = "time"
) 
Figure2A1

Figure2A2 <- FeaturePlot(
  object = dat,
  features = "HES1",
  split.by = "time"
) 
Figure2A2

Figure2A3 <- FeaturePlot(
  object = dat,
  features = "HES4",
  split.by = "time"
) 
Figure2A3

Figure2A4 <- FeaturePlot(
  object = dat,
  features = "HEY1",
  split.by = "time"
) 
Figure2A4

sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin24.4.0
# Running under: macOS Tahoe 26.3
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /opt/homebrew/Cellar/r/4.5.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1           patchwork_1.3.2   
# [5] ggplot2_4.0.2      Signac_1.16.0     
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-9           deldir_2.0-4           pbapply_1.7-4         
# [4] gridExtra_2.3          rlang_1.1.7            magrittr_2.0.4        
# [7] RcppAnnoy_0.0.23       otel_0.2.0             spatstat.geom_3.7-0   
# [10] matrixStats_1.5.0      ggridges_0.5.7         compiler_4.5.1        
# [13] reshape2_1.4.5         png_0.1-8              vctrs_0.7.1           
# [16] stringr_1.6.0          pkgconfig_2.0.3        crayon_1.5.3          
# [19] fastmap_1.2.0          XVector_0.50.0         labeling_0.4.3        
# [22] Rsamtools_2.26.0       promises_1.5.0         UCSC.utils_1.6.1      
# [25] purrr_1.2.1            GenomeInfoDb_1.46.2    jsonlite_2.0.0        
# [28] goftest_1.2-3          later_1.4.7            spatstat.utils_3.2-1  
# [31] BiocParallel_1.44.0    irlba_2.3.7            parallel_4.5.1        
# [34] cluster_2.1.8.2        R6_2.6.1               ica_1.0-3             
# [37] spatstat.data_3.1-9    stringi_1.8.7          RColorBrewer_1.1-3    
# [40] reticulate_1.45.0      spatstat.univar_3.1-6  parallelly_1.46.1     
# [43] lmtest_0.9-40          GenomicRanges_1.62.1   scattermore_1.2       
# [46] Rcpp_1.1.1             Seqinfo_1.0.0          RcppRoll_0.3.1        
# [49] tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15            
# [52] IRanges_2.44.0         sctransform_0.4.3      httpuv_1.6.16         
# [55] Matrix_1.7-4           splines_4.5.1          igraph_2.2.2          
# [58] tidyselect_1.2.1       abind_1.4-8            spatstat.random_3.4-4 
# [61] spatstat.explore_3.7-0 codetools_0.2-20       miniUI_0.1.2          
# [64] listenv_0.10.0         plyr_1.8.9             lattice_0.22-9        
# [67] tibble_3.3.1           shiny_1.13.0           withr_3.0.2           
# [70] S7_0.2.1               ROCR_1.0-12            Rtsne_0.17            
# [73] future_1.69.0          fastDummies_1.7.5      survival_3.8-6        
# [76] polyclip_1.10-7        fitdistrplus_1.2-6     Biostrings_2.78.0     
# [79] pillar_1.11.1          KernSmooth_2.23-26     stats4_4.5.1          
# [82] plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0        
# [85] S4Vectors_0.48.0       scales_1.4.0           globals_0.19.0        
# [88] xtable_1.8-8           glue_1.8.0             lazyeval_0.2.2        
# [91] tools_4.5.1            data.table_1.18.2.1    RSpectra_0.16-2       
# [94] RANN_2.6.2             dotCall64_1.2          fastmatch_1.1-8       
# [97] cowplot_1.2.0          grid_4.5.1             tidyr_1.3.2           
# [100] nlme_3.1-168           cli_3.6.5              spatstat.sparse_3.1-0 
# [103] spam_2.11-3            viridisLite_0.4.3      dplyr_1.2.0           
# [106] uwot_0.2.4             gtable_0.3.6           digest_0.6.39         
# [109] progressr_0.18.0       BiocGenerics_0.56.0    ggrepel_0.9.7         
# [112] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.9       
# [115] lifecycle_1.0.5        httr_1.4.8             mime_0.13             
# [118] MASS_7.3-65
