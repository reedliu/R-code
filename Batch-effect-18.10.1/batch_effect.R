# batch effect 也就是批次效应来源：
# 比如芯片数据，每次都要用机器读取，那么光照时间和强度每次都可能不一样， 极有可能出现批次效应
# 实验的三个重复，时间有间隔，也会有批次效应，如Western blot的三个重复
# 10年Nature综述：不同平台的数据，同一平台的不同时期的数据，同一个样品不同试剂的数据，以及同一个样品不同时间的数据等等都会产生一种batch effect。比如先测了疾病组，后测了对照组，分析得到的差异基因是和研究的因素有关还是和时间有关呢？
# 解决批次问题，不仅可以让实验更可靠，做多个芯片的联合分析；一般使用sva的中combat包来校正批次效应
source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("sva")
biocLite("bladderbatch")
library(sva)
library(bladderbatch) #膀胱癌样本

data(bladderdata)

#bladder的属性是EsetExpressionSet，所以可以用pData和exprs方法

pheno <- pData(bladderEset) # 注释信息

edata <- exprs(bladderEset) # 表达矩阵

#做一个组，用于批次效应中排除项
pheno$hasCancer <- as.numeric(pheno$cancer == "Cancer")

#先做一个聚类分析看一下
dist_mat <- dist(t(edata))
clustering <- hclust(dist_mat, method = "complete")
plot(clustering, labels = pheno$batch)
plot(clustering, labels = pheno$cancer)

#model 可有可无，有的话就是告诉combat这些分组本来就有差别
model <- model.matrix(~hasCancer, data = pheno)
combat_edata <- ComBat(dat = edata, batch = pheno$batch, mod = model)

dist_mat_combat <- dist(t(combat_edata))
clustering_combat <- hclust(dist_mat_combat,method = "complete")
plot(clustering_combat, labels = pheno$batch)
plot(clustering_combat, labels = pheno$cancer)

