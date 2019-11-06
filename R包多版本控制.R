### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-18
### Email: jieandze1314@gmail.com
### Title: R包的多版本控制
### ---------------
# 原本有4个R包路径，现在再添加一个
myPaths <- .libPaths()  
new <- c('~/r-2nd-pkgs')
myPaths <- c(myPaths, new) 
.libPaths(myPaths) 
.libPaths()
# 接下来调整位置，第一个是默认安装路径（位置随意调整）
myPaths <- .libPaths()  
myPaths <- c(myPaths[5], myPaths[2],myPaths[1],myPaths[3],myPaths[4])  
.libPaths(myPaths) 
.libPaths()
# 接下来就能在新的路径下安装其他版本了，例如安装seurat 2.3.4，它还需要其他依赖包
pkgs = c( 'mixtools', 'lars', 'dtw', 'doSNOW', 'hdf5r','fpc','foreach' ) 
BiocManager::install(pkgs,ask = F,update = F)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_2.3.4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

# 最后卸载之前的版本，加载新建版本
detach("package:Seurat", unload=TRUE)
library(Seurat, lib.loc="~/r-2nd-pkgs") 
packageVersion('Seurat')

