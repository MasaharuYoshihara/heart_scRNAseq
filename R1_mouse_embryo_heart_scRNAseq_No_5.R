library(dplyr)

dat <- read.csv("R1_mouse_scRNAseq_all_markers.csv", header = T)

dat <- subset(dat, gene == "Notch1")

dat <- dat[, 2:8]
dat <- dat[order(unique(dat$cluster)),]

write.csv(dat, "R1_mouse_scRNAseq_Notch1.csv")

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
#   [1] dplyr_1.2.0
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1  compiler_4.5.1    magrittr_2.0.4    R6_2.6.1          generics_0.1.4   
# [6] cli_3.6.6         tools_4.5.1       pillar_1.11.1     glue_1.8.0        rstudioapi_0.18.0
# [11] tibble_3.3.1      vctrs_0.7.1       lifecycle_1.0.5   pkgconfig_2.0.3   rlang_1.2.0