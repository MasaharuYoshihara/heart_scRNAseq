# embryonic heart at 12.5/14.5/16.5 days, Illumina NovaSeq 6000
# Ren et al., SciData, 2023, doi: 10.1038/s41597-023-02333-6

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(1)

dat_merge <- readRDS("20260425_mouse_embryo_heart_scRNAseq_after_pca.rds")

ElbowPlot(dat_merge)

dat_merge <- FindNeighbors(dat_merge, dims = 1:15)

# dat_merge_0.5 <- FindClusters(dat_merge, resolution = 0.5)
# dat_merge_0.5 <- RunUMAP(dat_merge_0.5, dims = 1:15)
# DimPlot(dat_merge_0.5, reduction = "umap", label = TRUE)

dat_merge_0.6 <- FindClusters(dat_merge, resolution = 0.6)
dat_merge_0.6 <- RunUMAP(dat_merge_0.6, dims = 1:15)
DimPlot(dat_merge_0.6, reduction = "umap", label = TRUE)

# dat_merge_0.65 <- FindClusters(dat_merge, resolution = 0.65)
# dat_merge_0.65 <- RunUMAP(dat_merge_0.65, dims = 1:15)
# DimPlot(dat_merge_0.65, reduction = "umap", label = TRUE)

# dat_merge_0.7 <- FindClusters(dat_merge, resolution = 0.7)
# dat_merge_0.7 <- RunUMAP(dat_merge_0.7, dims = 1:15)
# DimPlot(dat_merge_0.7, reduction = "umap", label = TRUE)

DimPlot(dat_merge_0.6, reduction = "umap", group.by = "orig.ident")

saveRDS(dat_merge_0.6, "20260425_mouse_embryo_heart_scRNAseq_after_UMAP.rds")

# sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin24.4.0
# Running under: macOS Tahoe 26.4.1
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
#   [1] future_1.69.0      ggplot2_4.0.2      dplyr_1.2.0        Seurat_5.4.0      
# [5] SeuratObject_5.3.0 sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
# [4] rlang_1.1.7            magrittr_2.0.4         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
# [10] compiler_4.5.1         spatstat.geom_3.7-0    png_0.1-8             
# [13] vctrs_0.7.1            reshape2_1.4.5         stringr_1.6.0         
# [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.5.0         purrr_1.2.1            jsonlite_2.0.0        
# [22] goftest_1.2-3          later_1.4.7            spatstat.utils_3.2-1  
# [25] irlba_2.3.7            parallel_4.5.1         cluster_2.1.8.2       
# [28] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [31] RColorBrewer_1.1-3     spatstat.data_3.1-9    reticulate_1.45.0     
# [34] parallelly_1.46.1      spatstat.univar_3.1-6  lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1          
# [40] future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.16          Matrix_1.7-4           splines_4.5.1         
# [46] igraph_2.2.2           tidyselect_1.2.1       abind_1.4-8           
# [49] codetools_0.2-20       spatstat.random_3.4-4  miniUI_0.1.2          
# [52] spatstat.explore_3.7-0 listenv_0.10.0         lattice_0.22-9        
# [55] tibble_3.3.1           plyr_1.8.9             withr_3.0.2           
# [58] shiny_1.13.0           S7_0.2.1               ROCR_1.0-12           
# [61] Rtsne_0.17             fastDummies_1.7.5      survival_3.8-6        
# [64] polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1         
# [67] KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4        
# [70] RcppHNSW_0.6.0         scales_1.4.0           globals_0.19.0        
# [73] xtable_1.8-8           glue_1.8.0             lazyeval_0.2.2        
# [76] tools_4.5.1            data.table_1.18.2.1    RSpectra_0.16-2       
# [79] RANN_2.6.2             dotCall64_1.2          cowplot_1.2.0         
# [82] grid_4.5.1             tidyr_1.3.2            nlme_3.1-168          
# [85] patchwork_1.3.2        cli_3.6.5              spatstat.sparse_3.1-0 
# [88] spam_2.11-3            viridisLite_0.4.3      uwot_0.2.4            
# [91] gtable_0.3.6           digest_0.6.39          progressr_0.18.0      
# [94] ggrepel_0.9.7          htmlwidgets_1.6.4      farver_2.1.2          
# [97] htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8            
# [100] mime_0.13              MASS_7.3-65