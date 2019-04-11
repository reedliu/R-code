###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-04-11
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### CAAS/AGIS/SDAU 
### Update Log:2019-04-1  AnnotationHub & Biostrings
### --------------- 

##################################
# AnnotationHub
##################################
library(AnnotationHub)
options()$BioC_mirror ##查看使用bioconductor的默认镜像
# 检查是否存在新版本
ahub <- AnnotationHub()




##################################
# BioStrings
##################################
# 准备基因组文件
BiocManager::install("BSgenome", version = "3.8")
library(BSgenome)
available.genomes() #总共91个













