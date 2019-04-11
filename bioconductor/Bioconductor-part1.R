### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-03-24
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### CAAS/AGIS/SDAU 
### Update Log: 2019-03-24  Bioconductor-Part1
### From: https://bioconductor.github.io/BiocWorkshops/r-and-bioconductor-for-everyone-an-introduction.html
### ---------------

####################################
## 关于R
####################################
# csv下载地址：https://raw.githubusercontent.com/Bioconductor/BiocWorkshops/master/100_Morgan_RBiocForAll/ALL-phenoData.csv

# 找到某个文件位置
fname <- file.choose()
fname
# 检查文件是否存在
file.exists(fname)
#读入文件
pdata <- read.csv(fname)
# 检查
dim(pdata)
summary(pdata)
head(pdata)
# 提前指定每列的类型
pdata <- read.csv(
  fname,
  colClasses = c("character", "factor", "integer", "factor")
)
summary(pdata)
# 检查A中是否包含B的信息
pdata$mol.biol %in% c("BCR/ABL", "NEG")
# 取子集
subset(pdata, sex == "F" & age > 50)
bcrabl <- subset(pdata, mol.biol %in% c("BCR/ABL", "NEG"))
dim( bcrabl )
# 然后table看一些(table的意思就是 tabular summary)
table(bcrabl$mol.biol)
# 因为我们只选了BCR/ABL和NEG两项，因此其他的都是空白，想要只看选出来的这两项信息，可以
factor(bcrabl$mol.biol)

# 更新一下bcrabl的mol.biol列（看前后Levels的变化）
bcrabl$mol.biol <- factor(bcrabl$mol.biol)

# formula使用
boxplot(age ~ mol.biol, bcrabl)
t.test(age ~ mol.biol, bcrabl)

# ggplot 
library(ggplot2)
ggplot(bcrabl, aes(x = mol.biol, y = age)) + geom_boxplot()

####################################
## 关于Bioconductor
####################################

# 示例bed下载：https://raw.githubusercontent.com/Bioconductor/BiocWorkshops/master/100_Morgan_RBiocForAll/CpGislands.Hsapiens.hg38.UCSC.bed

library(rtracklayer)
fname <- file.choose()
cpg <- import(fname)
cpg

# Subsetting operation that returns only the 'standard' chromosomes.
## Use pruning.mode="coarse" to drop list elements with mixed seqlevels:
cpg <- keepStandardChromosomes(cpg, pruning.mode = "coarse")
cpg

head( start(cpg) )
head( cpg$name )
hist(log10(width(cpg)))
subset(cpg, seqnames %in% c("chr1", "chr2"))

# 加载注释信息
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx

tx <- keepStandardChromosomes(tx, pruning.mode="coarse")
tx

olaps <- countOverlaps(tx, cpg)
length(olaps)     # 1 count per transcript
table(olaps)

tx$cpgOverlaps <- countOverlaps(tx, cpg)
tx

subset(tx, cpgOverlaps>10)

## SummarizedExperiment
fname <- file.choose()  # airway_colData.csv
colData <- read.csv(fname, row.names = 1)
colData
head(colData)

counts <- read.csv('airway_counts.csv', row.names=1)
counts <- as.matrix(counts)
dim(counts)
head(counts)

library("SummarizedExperiment")

se <- SummarizedExperiment(assay = counts, colData = colData)
se

subset(se, , dex == "trt")

assay(se) %>% head

colSums(assay(se))

# 添加进se metadata
se$lib.size <- colSums(assay(se))
colData(se)

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds
dds <- DESeq(dds)
head(results(dds))

# volcano
tmp <- as.data.frame(results(dds))

tmp$significant <- as.factor(tmp$pvalue<0.05 & abs(tmp$log2FoldChange) > 1)
ggplot(data=tmp, aes(x=log2FoldChange, y =-log10(pvalue),color =significant)) +
  geom_point() +
  scale_color_manual(values =c("black","red"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)



