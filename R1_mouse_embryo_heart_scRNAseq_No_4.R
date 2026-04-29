# embryonic heart at 12.5/14.5/16.5 days, Illumina NovaSeq 6000
# Ren et al., SciData, 2023, doi: 10.1038/s41597-023-02333-6

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(1)

dat_merge <- readRDS("20260425_mouse_embryo_heart_scRNAseq_after_UMAP.rds")

dat_merge$orig.ident <-paste0("E", dat_merge$orig.ident)

Figure2A <- DimPlot(
  object = dat_merge,
  group.by = 'orig.ident',
  label = FALSE,
  label.size = 6) + ggtitle('Sample') + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))+
  theme(plot.title = element_text(hjust = 0.5))
Figure2A # W700 H630

# DimPlot(dat_merge, reduction = "umap", group.by = "orig.ident")

DimPlot(dat_merge, reduction = "umap", label = TRUE)

FeaturePlot(dat_merge, "Cdh5")
FeaturePlot(dat_merge, "Myh7")
FeaturePlot(dat_merge, "Myh6")
FeaturePlot(dat_merge, "Nppa")
FeaturePlot(dat_merge, "Gja5")
FeaturePlot(dat_merge, "Acta1")
FeaturePlot(dat_merge, "Vim")
FeaturePlot(dat_merge, "Etv2")
FeaturePlot(dat_merge, "Sox17")
FeaturePlot(dat_merge, "Rgcc")
FeaturePlot(dat_merge, "Nrp2")
FeaturePlot(dat_merge, "Gata1")
FeaturePlot(dat_merge, "Klf1")
FeaturePlot(dat_merge, "Runx1")
FeaturePlot(dat_merge, "Tal1")
FeaturePlot(dat_merge, "Cldn6")
FeaturePlot(dat_merge, "Cldn7")
FeaturePlot(dat_merge, "Epcam")
FeaturePlot(dat_merge, "Rgs5")
FeaturePlot(dat_merge, "Hcn4")
FeaturePlot(dat_merge, "Tbx3")
FeaturePlot(dat_merge, "Ephb3")
FeaturePlot(dat_merge, "Cd68") # macrophage
FeaturePlot(dat_merge, "Cdh11") # endocardium
FeaturePlot(dat_merge, "Tnni3")
FeaturePlot(dat_merge, "Hand1")
FeaturePlot(dat_merge, "Tnnt1")
FeaturePlot(dat_merge, "Rspo3")
FeaturePlot(dat_merge, "Myl9")
FeaturePlot(dat_merge, "Kcne1l")
FeaturePlot(dat_merge, "Tgfb1")
FeaturePlot(dat_merge, "Lgr4")
FeaturePlot(dat_merge, "Lrp6")
FeaturePlot(dat_merge, "Col3a1")
FeaturePlot(dat_merge, "Sox10")
FeaturePlot(dat_merge, "Pax7")
FeaturePlot(dat_merge, "Sat1")
FeaturePlot(dat_merge, "Efnb2")
FeaturePlot(dat_merge, "Fabp4") # https://doi.org/10.1111/jcmm.12415 He et al., JCMM, 2014
FeaturePlot(dat_merge, "Hbb-bs")
FeaturePlot(dat_merge, "Hbb-bt")
FeaturePlot(dat_merge, "Ckm")
FeaturePlot(dat_merge, "Cdh2")
FeaturePlot(dat_merge, "Hbb-y")
FeaturePlot(dat_merge, "Tbx5")
FeaturePlot(dat_merge, "Kdr")
FeaturePlot(dat_merge, "Irx6")
FeaturePlot(dat_merge, "Myl4")
FeaturePlot(dat_merge, "Fabp5")

# dat_merge <- JoinLayers(dat_merge)
# temp <- FindMarkers(object = dat_merge, ident.1 = 4, ident.2 = 0)
# temp <- temp[order(-temp$avg_log2FC), ]
# temp <- temp[order(-temp$pct.1), ]
# temp <- temp[order(temp$p_val_adj), ]
# write.csv(temp, "cluster_7.csv")

markers_temp <- c("Myl4", "Nppa", "Myh7", "Myh6", "Ckm","Rspo3", "Kcne1l", 
                  "Cdh5", "Pecam1", "Fabp4", "Fabp5", "Nfatc1", "Cdh2", 
                  "Cdh11", "Hbb-bt", "Cd68")
DotPlot(dat_merge, features = markers_temp)

new.cluster.ids <- c("Atrial Myocyte 1", "Atrial Myocyte 2", 
                     "Ventricular Myocyte 1", "Ventricular Myocyte 2",
                     "Endothelial Cell 1", "Ventricular Myocyte 3", 
                     "Endothelial Cell 2", "Ventricular Myocyte 4",
                     "AVC 1", "Epicardial Cell 1",
                     "Ventricular Myocyte 5", "AVC 2",
                     "Ventricular Myocyte 6", "Atrial Myocyte 3",
                     "Erythroild Cell 1", "Epicardial Cell 2",
                     "Erythroild Cell 2", "Endothelial Cell 3",
                     "Endothelial Cell 4", "Ventricular Myocyte 7",
                     "Macrophage", "Atrial Myocyte 4",
                     "Erythroild Cell 3", "Endothelial Cell 5")
names(new.cluster.ids) <- levels(dat_merge)
dat_merge <- RenameIdents(dat_merge, new.cluster.ids)
dat_merge$celltype <- Idents(dat_merge)

Figure2B <- DimPlot(dat_merge, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle('Cell type') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")  + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure2B # W630 H670

cluster_data <- table(Idents(dat_merge), dat_merge$orig.ident)
head(cluster_data)
write.csv(cluster_data, "R1_mouse_scRNAseq_cluster_data.csv")

Figure2C <- FeaturePlot(dat_merge, features = "Cdh5")+ 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure2C # W630 H670

markers <- c("Myl4", "Nppa", "Myh7", "Myh6", "Ckm","Rspo3", "Kcne1l", 
             "Cdh5", "Pecam1", "Fabp4", "Fabp5", "Nfatc1", "Cdh2", 
             "Cdh11", "Hbb-bt", "Cd68",
             "Notch1", "Notch2", "Notch3", "Notch4", 
             "Hes1", "Hes2", "Hes3", "Hes5", "Hes6", "Hes7", 
             "Hey1", "Hey2", "Heyl"
             )

Figure2D <- DotPlot(dat_merge, features = markers) + 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20)) +
  RotatedAxis() +
  scale_y_discrete(
    limits = c("Atrial Myocyte 1", "Atrial Myocyte 2", 
               "Atrial Myocyte 3", "Atrial Myocyte 4",
               "Ventricular Myocyte 1", "Ventricular Myocyte 2", "Ventricular Myocyte 3", 
               "Ventricular Myocyte 4", "Ventricular Myocyte 5", "Ventricular Myocyte 6", 
               "Ventricular Myocyte 7", 
               "AVC 1", "AVC 2", "Epicardial Cell 1", "Epicardial Cell 2", 
               "Endothelial Cell 1", "Endothelial Cell 2", "Endothelial Cell 3",
               "Endothelial Cell 4", "Endothelial Cell 5",
               "Erythroild Cell 1", "Erythroild Cell 2", "Erythroild Cell 3",
               "Macrophage"
               )
  )
# The following requested variables were not found: Hes2, Hes3 
Figure2D # W1900 H1000

Figure2D_data <- Figure2D$data
head(Figure2D_data)
Figure2D_data_2 <- Figure2D_data[,2:5]
head(Figure2D_data_2)
write.csv(Figure2D_data_2, "R1_mouse_scRNAseq_avgexp_pctexp.csv")

Figure2E <- FeaturePlot(dat_merge, features = "Notch1")+ 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure2E # W500 H530

Figure2F <- FeaturePlot(dat_merge, features = "Hes1")+ 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure2F # W500 H530

Figure2G <- FeaturePlot(dat_merge, features = "Fabp4")+ 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure2G # W500 H530

Figure2H <- FeaturePlot(dat_merge, features = "Nfatc1")+ 
  FontSize(x.title = 20, y.title = 20, x.text = 20, y.text = 20) +
  theme(legend.text = element_text(size = 20), plot.title = element_text(size = 20))
Figure2H # W500 H530

dat_merge <- JoinLayers(dat_merge)
all.markers <- FindAllMarkers(object = dat_merge, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
write.csv(all.markers, "R1_mouse_scRNAseq_all_markers.csv")
top5_all <- all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
write.csv(top5_all, "R1_mouse_scRNAseq_top5_DEGs.csv")

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
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.2.0           
# [5] magrittr_2.0.4         RcppAnnoy_0.0.23       otel_0.2.0             matrixStats_1.5.0     
# [9] ggridges_0.5.7         compiler_4.5.1         spatstat.geom_3.7-0    png_0.1-8             
# [13] vctrs_0.7.1            reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3       
# [17] fastmap_1.2.0          labeling_0.4.3         promises_1.5.0         purrr_1.2.1           
# [21] jsonlite_2.0.0         goftest_1.2-3          later_1.4.7            spatstat.utils_3.2-1  
# [25] irlba_2.3.7            parallel_4.5.1         cluster_2.1.8.2        R6_2.6.1              
# [29] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-9   
# [33] limma_3.66.0           reticulate_1.45.0      parallelly_1.46.1      spatstat.univar_3.1-6 
# [37] lmtest_0.9-40          scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1          
# [41] future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3      httpuv_1.6.16         
# [45] Matrix_1.7-4           splines_4.5.1          igraph_2.2.2           tidyselect_1.2.1      
# [49] rstudioapi_0.18.0      abind_1.4-8            codetools_0.2-20       spatstat.random_3.4-4 
# [53] miniUI_0.1.2           spatstat.explore_3.7-0 listenv_0.10.0         lattice_0.22-9        
# [57] tibble_3.3.1           plyr_1.8.9             withr_3.0.2            shiny_1.13.0          
# [61] S7_0.2.1               ROCR_1.0-12            Rtsne_0.17             future_1.69.0         
# [65] fastDummies_1.7.5      survival_3.8-6         polyclip_1.10-7        fitdistrplus_1.2-6    
# [69] pillar_1.11.1          KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4        
# [73] RcppHNSW_0.6.0         scales_1.4.0           globals_0.19.0         xtable_1.8-8          
# [77] glue_1.8.0             lazyeval_0.2.2         tools_4.5.1            data.table_1.18.2.1   
# [81] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.2.0         
# [85] grid_4.5.1             tidyr_1.3.2            nlme_3.1-168           patchwork_1.3.2       
# [89] cli_3.6.6              spatstat.sparse_3.1-0  spam_2.11-3            viridisLite_0.4.3     
# [93] uwot_0.2.4             gtable_0.3.6           digest_0.6.39          progressr_0.18.0      
# [97] ggrepel_0.9.7          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.9       
# [101] lifecycle_1.0.5        httr_1.4.8             statmod_1.5.1          mime_0.13             
# [105] MASS_7.3-65