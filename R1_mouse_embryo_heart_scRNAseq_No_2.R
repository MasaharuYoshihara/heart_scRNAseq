# embryonic heart at 12.5/14.5/16.5 days, Illumina NovaSeq 6000
# Ren et al., SciData, 2023, doi: 10.1038/s41597-023-02333-6


library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(1)

dat_merge <- readRDS("20260424_mouse_embryo_heart_scRNAseq_merge_before_QC.rds")

# dat_merge
# An object of class Seurat 
# 18378 features across 45135 samples within 1 assay 
# Active assay: RNA (18378 features, 0 variable features)
# 6 layers present: counts.16.5_2, counts.12.5_1, counts.12.5_2, counts.14.5_1, counts.14.5_2, counts.16.5_1

dat_merge[["percent.mt"]] <- PercentageFeatureSet(dat_merge, pattern = "^mt-")

VlnPlot(dat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dat_merge <- subset(dat_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

dat_merge <- NormalizeData(dat_merge, normalization.method = "LogNormalize", scale.factor = 10000)

dat_merge <- FindVariableFeatures(dat_merge, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(dat_merge)
dat_merge <- ScaleData(dat_merge, features = all.genes)

dat_merge <- RunPCA(dat_merge, features = VariableFeatures(object = dat_merge))

# dat_merge
# An object of class Seurat 
# 18378 features across 42447 samples within 1 assay 
# Active assay: RNA (18378 features, 2000 variable features)
# 13 layers present: counts.16.5_2, counts.12.5_1, counts.12.5_2, counts.14.5_1, counts.14.5_2, counts.16.5_1, data.16.5_2, data.12.5_1, data.12.5_2, data.14.5_1, data.14.5_2, data.16.5_1, scale.data
# 1 dimensional reduction calculated: pca

saveRDS(dat_merge, "20260425_mouse_embryo_heart_scRNAseq_after_pca.rds")

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
#   [1] ggplot2_4.0.2      dplyr_1.2.0        Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.1.7           
# [5] magrittr_2.0.4         RcppAnnoy_0.0.23       otel_0.2.0             matrixStats_1.5.0     
# [9] ggridges_0.5.7         compiler_4.5.1         spatstat.geom_3.7-0    png_0.1-8             
# [13] vctrs_0.7.1            reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3       
# [17] fastmap_1.2.0          labeling_0.4.3         promises_1.5.0         purrr_1.2.1           
# [21] jsonlite_2.0.0         goftest_1.2-3          later_1.4.7            spatstat.utils_3.2-1  
# [25] irlba_2.3.7            parallel_4.5.1         cluster_2.1.8.2        R6_2.6.1              
# [29] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-9   
# [33] reticulate_1.45.0      parallelly_1.46.1      spatstat.univar_3.1-6  lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2   
# [41] zoo_1.8-15             sctransform_0.4.3      httpuv_1.6.16          Matrix_1.7-4          
# [45] splines_4.5.1          igraph_2.2.2           tidyselect_1.2.1       abind_1.4-8           
# [49] codetools_0.2-20       spatstat.random_3.4-4  miniUI_0.1.2           spatstat.explore_3.7-0
# [53] listenv_0.10.0         lattice_0.22-9         tibble_3.3.1           plyr_1.8.9            
# [57] withr_3.0.2            shiny_1.13.0           S7_0.2.1               ROCR_1.0-12           
# [61] Rtsne_0.17             future_1.69.0          fastDummies_1.7.5      survival_3.8-6        
# [65] polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1          KernSmooth_2.23-26    
# [69] plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0          
# [73] globals_0.19.0         xtable_1.8-8           glue_1.8.0             lazyeval_0.2.2        
# [77] tools_4.5.1            data.table_1.18.2.1    RSpectra_0.16-2        RANN_2.6.2            
# [81] dotCall64_1.2          cowplot_1.2.0          grid_4.5.1             tidyr_1.3.2           
# [85] nlme_3.1-168           patchwork_1.3.2        cli_3.6.5              spatstat.sparse_3.1-0 
# [89] spam_2.11-3            viridisLite_0.4.3      uwot_0.2.4             gtable_0.3.6          
# [93] digest_0.6.39          progressr_0.18.0       ggrepel_0.9.7          htmlwidgets_1.6.4     
# [97] farver_2.1.2           htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8            
# [101] mime_0.13              MASS_7.3-65