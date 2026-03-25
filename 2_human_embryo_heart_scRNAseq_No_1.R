# GSE181346

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

RNAseq <- readRDS("human_6_8_12and19_merged_final_cleaned.rds")
# RNAseq2 <- readRDS("human_6_8_12and19_merged_final.rds")

head(RNAseq)
# head(RNAseq2)

# The following line helps run
RNAseq@graphs <- list()

RNAseq <- UpdateSeuratObject(object = RNAseq)
# DimPlot(RNAseq, reduction = "umap", label = TRUE, group.by = "new_manual_annotation", repel = TRUE)

# RNAseq <- FindNeighbors(RNAseq, dims = 1:15)

library(readxl)
library(tibble)
metadata <- read_excel("NIHMS1858872-supplement-7.xlsx", sheet = "E")
head(metadata)
metadata_to_add <- metadata %>%
  select(Sample, `Cluster names`, Barcode) %>%
  column_to_rownames("Barcode")

matching_barcodes <- rownames(metadata_to_add) %in% rownames(RNAseq@meta.data)
summary(matching_barcodes)

metadata_to_add <- metadata[, c("Sample", "Cluster names")]
rownames(metadata_to_add) <- metadata$Barcode
head(metadata_to_add)
metadata_to_add <- metadata_to_add[match(rownames(RNAseq@meta.data), rownames(metadata_to_add)), ]
RNAseq <- AddMetaData(RNAseq, metadata = metadata_to_add)
head(RNAseq)

colnames(RNAseq@meta.data)[which(colnames(RNAseq@meta.data) == "Cluster names")] <- "Cluster_names"
Idents(RNAseq) <- RNAseq$`Cluster_names`

DimPlot(RNAseq, reduction = "umap", label = TRUE, group.by = "Cluster_names", repel = TRUE)
DimPlot(RNAseq, reduction = "umap", label = TRUE, group.by = "Sample", repel = TRUE)

c_features <- c("TTN", "MYH7", "ANKRD1", "TNNT2", "NPPA", "MYL2", "DSP",
                "PKP2", "HAND1", "MYOCD", "KCNH2", "TNNI3", "ALPK2",
                "SCN5A", "MYL4", "MYL3", "GJA5", "TBX5", "MYL7", "MYBPC3",
                "KCNJ5", "TCAP", "MYH6",
                "PECAM1", "CDH5", "APLNR", "SOX17", "KDR", "DLL4", "FLT1",
                "GATA2", "SELE", "EMCN", "MEOX2", "TEK", "NDNF", "JCAD",
                "ROBO4", "MYDGF", "COL1A1", "MFAP4", "TWIST1","GDNF", "TAGLN",
                "PDGFRB", "MYH11", "ACTA2",
                "DCN", "LUM", "TBX18", "TCF21",
                "WT1", "COL3A1", "VCAM1", "DDR2")
DotPlot(RNAseq, features = c_features) + 
  FontSize(x.title = 30, y.title = 30, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis()

saveRDS(object = RNAseq, "260228_Heart_ATAC_F6_8_12_19_merge_No0.rds")

sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8   
# [3] LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tibble_3.2.1       readxl_1.4.5       ggplot2_4.0.0      patchwork_1.3.2   
# [5] Seurat_5.3.0       SeuratObject_5.0.2 sp_2.2-0           dplyr_1.1.4       
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
# [4] rlang_1.1.5            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
# [10] compiler_4.4.1         spatstat.geom_3.3-5    png_0.1-8             
# [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.6.0         
# [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] utf8_1.2.6             promises_1.5.0         purrr_1.0.4           
# [22] jsonlite_1.9.0         goftest_1.2-3          later_1.4.1           
# [25] spatstat.utils_3.1-2   irlba_2.3.5.1          parallel_4.4.1        
# [28] cluster_2.1.6          R6_2.6.1               ica_1.0-3             
# [31] stringi_1.8.4          RColorBrewer_1.1-3     spatstat.data_3.1-9   
# [34] reticulate_1.41.0      parallelly_1.45.1      spatstat.univar_3.1-1 
# [37] cellranger_1.1.0       lmtest_0.9-40          scattermore_1.2       
# [40] Rcpp_1.0.14            tensor_1.5.1           future.apply_1.20.2   
# [43] zoo_1.8-13             sctransform_0.4.1      httpuv_1.6.15         
# [46] Matrix_1.7-3           splines_4.4.1          igraph_2.1.4          
# [49] tidyselect_1.2.1       rstudioapi_0.18.0      dichromat_2.0-0.1     
# [52] abind_1.4-8            spatstat.random_3.3-2  codetools_0.2-20      
# [55] miniUI_0.1.2           spatstat.explore_3.3-4 listenv_0.10.0        
# [58] lattice_0.22-6         plyr_1.8.9             withr_3.0.2           
# [61] shiny_1.13.0           S7_0.2.0               ROCR_1.0-12           
# [64] Rtsne_0.17             future_1.69.0          fastDummies_1.7.5     
# [67] survival_3.6-4         polyclip_1.10-7        fitdistrplus_1.2-6    
# [70] pillar_1.11.1          KernSmooth_2.23-24     plotly_4.12.0         
# [73] generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0          
# [76] globals_0.19.0         xtable_1.8-8           glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.4.1            data.table_1.17.0     
# [82] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2         
# [85] cowplot_1.2.0          grid_4.4.1             tidyr_1.3.1           
# [88] nlme_3.1-164           cli_3.6.4              spatstat.sparse_3.1-0 
# [91] spam_2.11-1            viridisLite_0.4.3      uwot_0.2.3            
# [94] gtable_0.3.6           digest_0.6.37          progressr_0.18.0      
# [97] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2          
# [100] htmltools_0.5.8.1      lifecycle_1.0.5        httr_1.4.8            
# [103] mime_0.12              MASS_7.3-60.2