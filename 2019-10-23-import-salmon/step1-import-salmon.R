### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: 使用tximport导入salmon数据
### ---------------

# tximport函数主要需要两个参数：定量文件files和转录本与基因名的对应文件tx2gene
rm(list = ls())
options(stringsAsFactors = F)

####################
# 配置files路径
####################
dir <- file.path(getwd(),'quant/')
dir
files <- list.files(pattern = '*sf',dir,recursive = T)
files <- file.path(dir,files)
all(file.exists(files))

####################
# 配置tx2gene
####################
# https://support.bioconductor.org/p/101156/
# BiocManager::install("EnsDb.Mmusculus.v79")
# 如何构建？
if(F){
  library(EnsDb.Mmusculus.v79)
  txdf <- transcripts(EnsDb.Mmusculus.v79, return.type="DataFrame")
  mm10_tx2gene <- as.data.frame(txdf[,c("tx_id", "gene_id")])
  head(mm10_tx2gene)
  write.csv(mm10_tx2gene,file = 'mm10_tx2gene.csv')
}

tx2_gene_file <- 'mm10_tx2gene.csv'
tx2gene <- read.csv(tx2_gene_file,row.names = 1)
head(tx2gene)

####################
# 开始整合
####################
library(tximport)
library(readr)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,ignoreTxVersion = T)
names(txi)
head(txi$length)
head(txi$counts)
# 发现目前counts的列名还没有指定

####################
# 添加列名
####################
# 获取导入文件名称的ID，如SRR7815870
library(stringr)
# 先得到SRRxxxx_quant这一部分
n1 <- sapply(strsplit(files,'\\/'), function(x)x[12])
# 或者
n1 <- sapply(strsplit(files,'\\/'), '[',12)
# 再或者
n1 <- stringr::str_split(files,'\\/', simplify = T)[,12]

# 再得到SRRxxxx这一部分
n2 <- sapply(strsplit(n1,'_'), function(x)x[1])
# 或者
n2 <- sapply(strsplit(n1,'_'), '[',1)
# 再或者
n2 <- stringr::str_split(n1,'_',simplify = T)[,1]

colnames(txi$counts) <- n2
head(txi$counts)

####################
# 操作新得到的表达矩阵
####################
salmon_expr <- txi$counts
# 表达量取整
salmon_expr <- apply(salmon_expr, 2, as.integer)
head(salmon_expr)
# 添加行名
rownames(salmon_expr) <- rownames(txi$counts)
dim(salmon_expr)

save(salmon_expr,file = 'salmon-aggr.Rdata')



