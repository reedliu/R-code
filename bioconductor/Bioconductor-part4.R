### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-04-08
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### CAAS/AGIS/SDAU 
### Update Log:2019-04-08  FindOverlaps and GRanges
### ---------------
rm(list=ls())
options(stringsAsFactors = F)

library(IRanges)
qry <- IRanges(start=c(1, 26, 19, 11, 21, 7), 
               end=c(16, 30, 19, 15, 24, 8), 
               names=letters[1:6])

sbj <- IRanges(start=c(1, 19, 10), 
               end=c(5, 29, 16), 
               names=letters[24:26])

# 默认any比对模式
hts <- findOverlaps(qry, sbj)
hts
names(qry)[queryHits(hts)]
names(sbj)[subjectHits(hts)]
# 使用within模式
hts_within <- findOverlaps(qry, sbj, type="within")
hts_within

# 使用select参数
findOverlaps(qry, sbj,select = "first")
findOverlaps(qry, sbj,select = "last")
findOverlaps(qry, sbj,select = "arbitrary")

# 创建GRanges
library(GenomicRanges)
gr1 <- GRanges(seq = c("chr1","chr2","chr3","chr1"),
               ranges = IRanges(start = c(5:8),width = 10),
               strand = "+")
gr1

gr2 <- GRanges(seq = c("chr1","chr2","chr3","chr1"),
               ranges = IRanges(start = c(5:8),width = 10),
               strand = "+",
               gc = seq(10,70,length.out = 4),
               seqlengths = c(chr1=150,chr2=200,chr3=250))
gr2

# 访问seqname
seqnames(gr2)
# ranges
ranges(gr2)
start(gr2)
end(gr2)
# 这里的width能做什么？=》统计全基因组范围内的基因长度分布
width(gr2)
# 看一下总共有多少行
length(gr2)

# strand
strand(gr2)

length(gr2)
names(gr2) <- letters[1:length(gr2)]
gr2

# 取子集
gr2[start(gr2)>7]

# 看看metadata columns
mcols(gr2)
gr2$gc

# 计算1号染色体的GC content
mean(gr2[seqnames(gr2) == "chr1"]$gc)
sort(gr2) #实现排序
gr2[order(gr2$gc)] #按某一列进行排序(从小到大)
gr2[order(gr2$gc,decreasing = T)] #按某一列进行排序(从大到小)

gr5 <- c(gr2,gr2)
length(gr2)
length(gr5)
gr5
unique(gr5)

# 拆分GRanges
gr_split <- split(gr2,seqnames(gr2))
names(gr_split)
gr_split[[1]]

## GRangesList(注意Rle的方便之处)
gr3 <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(2, 3)),
               ranges = IRanges(1:5, end = 6:10),
               strand = Rle(strand(c("-", "+", "+","-")), c(1,1,2,1)),
               score = 1:5, GC = seq(1, 0, length = 5))
gr4 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
grl <- GRangesList(gr3, gr4)
grl
unlist(grl)

dbgrl <- c(grl,grl) #make double
dbgrl 
length(dbgrl)

#循环
sapply(gr_split, function(x) min(start(x)))
sapply(gr_split, length)

##练习：统计human exon total length
# GRCh36 (hg18): ENSEMBL release_52.
# GRCh37 (hg19): ENSEMBL release_59/61/64/68/69/75.
# GRCh38 (hg38): ENSEMBL  release_76/77/78/80/81/82...

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
BiocManager::install("EnsDb.Hsapiens.v86", version = "3.8")
library(EnsDb.Hsapiens.v86)

ensembl.hg38 <- EnsDb.Hsapiens.v86
ensembl.hg38.exon <- exons(ensembl.hg38)
# 看下长度
length(ensembl.hg38.exon)
# 长度分布-v1
hist(width(ensembl.hg38.exon))
# 长度分布-v2(上面的图中受极大值影响，大部分数据展示不出来)
hist(log2(width(ensembl.hg38.exon)+1)) #大部分在2^7~2^8左右，也就是128-256bp之间

# 数据探索
mean(width(ensembl.hg38.exon)) #因此human exon 平均长度是310
median(width(ensembl.hg38.exon))

sum(width(ensembl.hg38.exon)) #错误做法
sum(width(GenomicRanges::reduce(ensembl.hg38.exon))) #正确做法










