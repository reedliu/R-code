remove.packages('Seurat')
pkgs = c( 'mixtools', 'lars', 'dtw', 'doSNOW', 'hdf5r' ) 
#pkgs=c('jackstraw','slingshot')
BiocManager::install(pkgs,ask = F,update = F)
# 以后只需要修改这个版本号即可
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.0.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
library(Seurat)

