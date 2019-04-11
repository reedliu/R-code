### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-03-24
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### CAAS/AGIS/SDAU 
### Update Log: 2019-03-24  Bioconductor-Part2
### From: https://bioconductor.github.io/BiocWorkshops/r-and-bioconductor-for-everyone-an-introduction.html
### ---------------

# 示例数据eset.RData下载：https://github.com/jmacdon/Bioc2018Anno/blob/master/inst/extdata/eset.Rdata
rm(list=ls())
options(stringsAsFactors = F)
load('eset.Rdata')
eset

head(exprs(eset))
head(pData(phenoData(eset)))
head(pData(featureData(eset)))

# Simple example
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
BiocManager::install("hugene20sttranscriptcluster.db", version = "3.8")
library(hugene20sttranscriptcluster.db)
set.seed(12345)
ids <- featureNames(eset)[sample(1:25000, 5)]
ids
# 注释
AnnotationDbi::select(hugene20sttranscriptcluster.db, ids, "SYMBOL")
# 换个包
if (! require(hgu95av2.db,character.only=T) ) {
  BiocManager::install(hgu95av2.db,ask = F,update = F)
}
columns(hgu95av2.db)
keys <- head( keys(hgu95av2.db) )
keys
AnnotationDbi::select(hgu95av2.db, keys=keys, columns = c("SYMBOL","ENTREZID"))

# 另外一个例子
ids <- c('16737401','16657436' ,'16678303')
AnnotationDbi::select(hugene20sttranscriptcluster.db, ids, c("SYMBOL","MAP"))

mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID")

mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID",multiVals = "list")
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID",multiVals = "CharacterList")
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID",multiVals = "filter")
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID",multiVals = "asNA")

# all keytypes
columns(hugene20sttranscriptcluster.db)
# What gene symbol corresponds to Entrez Gene ID 1000?
mapIds(org.Hs.eg.db, "1000", "SYMBOL", "ENTREZID")
# What is the Ensembl Gene ID for PPARG?
mapIds(org.Hs.eg.db, "PPARG", "ENTREZID", "SYMBOL")
# What is the UniProt ID for GAPDH?
mapIds(org.Hs.eg.db, "GAPDH", "UNIPROT", "SYMBOL")

# How many of the probesets from the ExpressionSet (eset) we loaded map to a single gene? How many don’t map to a gene at all?
mpid <- mapIds(hugene20sttranscriptcluster.db, rownames(eset), "ENTREZID", "PROBEID") 
length(mpid)
# only once
mpid[table(mpid)==1] %>% length()


tmp <- as.data.frame(mpid)
sum(is.na(tmp$mpid))


require("TxDb.Hsapiens.UCSC.hg19.knownGene")
AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene, c("1","10"),
       c("TXNAME","TXCHROM","TXSTART","TXEND"), "GENEID")

library("EnsDb.Hsapiens.v86")
AnnotationDbi::select(EnsDb.Hsapiens.v79, c("1", "10"),
       c("GENEID","GENENAME","SEQNAME","GENESEQSTART","GENESEQEND"), "ENTREZID")








