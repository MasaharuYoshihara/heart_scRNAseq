# embryonic heart at 12.5/14.5/16.5 days, Illumina NovaSeq 6000
# Ren et al., SciData, 2023, doi: 10.1038/s41597-023-02333-6

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(1)

count_12.5_1 <- ReadMtx(mtx = "GSM7226270_E12_5_1_matrix.mtx.gz",
                        features = "GSM7226270_E12_5_1_features.tsv.gz",
                        cells = "GSM7226270_E12_5_1_barcodes.tsv.gz")

seurat_12.5_1 <- CreateSeuratObject(count_12.5_1, project = "12.5_1",
                                    min.cells = 3, min.features = 200)

# seurat_12.5_1
# An object of class Seurat 
# 17092 features across 6962 samples within 1 assay 
# Active assay: RNA (17092 features, 0 variable features)
# 1 layer present: counts

count_12.5_2 <- ReadMtx(mtx = "GSM7226271_E12_5_2_matrix.mtx.gz",
                        features = "GSM7226271_E12_5_2_features.tsv.gz",
                        cells = "GSM7226271_E12_5_2_barcodes.tsv.gz")

seurat_12.5_2 <- CreateSeuratObject(count_12.5_2, project = "12.5_2",
                                    min.cells = 3, min.features = 200)

# seurat_12.5_2
# An object of class Seurat 
# 17381 features across 9846 samples within 1 assay 
# Active assay: RNA (17381 features, 0 variable features)
# 1 layer present: counts

count_14.5_1 <- ReadMtx(mtx = "GSM7226272_E14_5_1_matrix.mtx.gz",
                        features = "GSM7226272_E14_5_1_features.tsv.gz",
                        cells = "GSM7226272_E14_5_1_barcodes.tsv.gz")

seurat_14.5_1 <- CreateSeuratObject(count_14.5_1, project = "14.5_1",
                                    min.cells = 3, min.features = 200)

# seurat_14.5_1
# An object of class Seurat 
# 17014 features across 6983 samples within 1 assay 
# Active assay: RNA (17014 features, 0 variable features)
# 1 layer present: counts

count_14.5_2 <- ReadMtx(mtx = "GSM7226273_E14_5_2_matrix.mtx.gz",
                        features = "GSM7226273_E14_5_2_features.tsv.gz",
                        cells = "GSM7226273_E14_5_2_barcodes.tsv.gz")

seurat_14.5_2 <- CreateSeuratObject(count_14.5_2, project = "14.5_2",
                                    min.cells = 3, min.features = 200)

# seurat_14.5_2
# An object of class Seurat 
# 16602 features across 6673 samples within 1 assay 
# Active assay: RNA (16602 features, 0 variable features)
# 1 layer present: counts

count_16.5_1 <- ReadMtx(mtx = "GSM7226274_E16_5_1_matrix.mtx.gz",
                        features = "GSM7226274_E16_5_1_features.tsv.gz",
                        cells = "GSM7226274_E16_5_1_barcodes.tsv.gz")

seurat_16.5_1 <- CreateSeuratObject(count_16.5_1, project = "16.5_1",
                                    min.cells = 3, min.features = 200)

# seurat_16.5_1
# An object of class Seurat 
# 16920 features across 7321 samples within 1 assay 
# Active assay: RNA (16920 features, 0 variable features)
# 1 layer present: counts

count_16.5_2 <- ReadMtx(mtx = "GSM7226276_E16_5_2_matrix.mtx.gz",
                        features = "GSM7226276_E16_5_2_features.tsv.gz",
                        cells = "GSM7226276_E16_5_2_barcodes.tsv.gz")

seurat_16.5_2 <- CreateSeuratObject(count_16.5_2, project = "16.5_2",
                                    min.cells = 3, min.features = 200)

# seurat_16.5_2
# An object of class Seurat 
# 16595 features across 7350 samples within 1 assay 
# Active assay: RNA (16595 features, 0 variable features)
# 1 layer present: counts

merge <- merge(x = seurat_16.5_2,
               y = c(seurat_12.5_1, seurat_12.5_2,
                     seurat_14.5_1, seurat_14.5_2,
                     seurat_16.5_1),
               add.cell.ids = c("16.5_2",
                                "12.5_1", "12.5_2",
                                "14.5_1", "14.5_2",
                                "16.5_1"))

saveRDS(merge, "20260424_mouse_embryo_heart_scRNAseq_merge_before_QC.rds")

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
# [17] fastmap_1.2.0          promises_1.5.0         purrr_1.2.1            jsonlite_2.0.0        
# [21] goftest_1.2-3          later_1.4.7            spatstat.utils_3.2-1   irlba_2.3.7           
# [25] parallel_4.5.1         cluster_2.1.8.2        R6_2.6.1               ica_1.0-3             
# [29] stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-9    reticulate_1.45.0     
# [33] parallelly_1.46.1      spatstat.univar_3.1-6  lmtest_0.9-40          scattermore_1.2       
# [37] Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15            
# [41] sctransform_0.4.3      httpuv_1.6.16          Matrix_1.7-4           splines_4.5.1         
# [45] igraph_2.2.2           tidyselect_1.2.1       abind_1.4-8            codetools_0.2-20      
# [49] spatstat.random_3.4-4  miniUI_0.1.2           spatstat.explore_3.7-0 listenv_0.10.0        
# [53] lattice_0.22-9         tibble_3.3.1           plyr_1.8.9             withr_3.0.2           
# [57] shiny_1.13.0           S7_0.2.1               ROCR_1.0-12            Rtsne_0.17            
# [61] future_1.69.0          fastDummies_1.7.5      survival_3.8-6         polyclip_1.10-7       
# [65] fitdistrplus_1.2-6     pillar_1.11.1          KernSmooth_2.23-26     plotly_4.12.0         
# [69] generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0           globals_0.19.0        
# [73] xtable_1.8-8           glue_1.8.0             lazyeval_0.2.2         tools_4.5.1           
# [77] data.table_1.18.2.1    RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2         
# [81] cowplot_1.2.0          grid_4.5.1             tidyr_1.3.2            nlme_3.1-168          
# [85] patchwork_1.3.2        cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-3           
# [89] viridisLite_0.4.3      uwot_0.2.4             gtable_0.3.6           digest_0.6.39         
# [93] progressr_0.18.0       ggrepel_0.9.7          htmlwidgets_1.6.4      farver_2.1.2          
# [97] htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8             mime_0.13             
# [101] MASS_7.3-65
