### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: 和featureCount结果比较
### ---------------
rm(list = ls())
options(stringsAsFactors = F)
load('salmon-aggr.Rdata')
dim(salmon_expr)
colnames(salmon_expr)
salmon_expr[1:3,1:3]

# 为了方便比较，我们选第一个样本SRR7815790
salmon_SRR7815790 <- salmon_expr[,1]
names(salmon_SRR7815790) <- rownames(salmon_expr)

###########################
# 然后找featureCounts结果
###########################
ft <- read.table('../hisat/hisat2_counts.txt',row.names = 1,comment.char = '#',header = T)
ft[1:4,1:4]
colnames(ft)
# 两种取子集的方法，grep返回结果的序号；grepl返回逻辑值
ft_SRR7815790 <-  ft[,grep(colnames(salmon_expr)[1], colnames(ft))]
ft_SRR7815790 <-  ft[,grepl(colnames(salmon_expr)[1], colnames(ft))]

names(ft_SRR7815790) <- rownames(ft)

###########################
# 取salmon和featureCounts公共基因
###########################
# salmon使用的Ensembl基因，而featureCounts得到的是Symbol基因
# 先进行基因名转换
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
gene_tr <- clusterProfiler::bitr(names(salmon_SRR7815790),
                      "ENSEMBL","SYMBOL",
                      org.Mm.eg.db)
# 要找salmon在gene_tr中的对应位置，然后取gene_tr的第二列symbol，因此match中把salmon放第一个位置
names(salmon_SRR7815790) <- gene_tr[match(names(salmon_SRR7815790),gene_tr$ENSEMBL),2]

# 找共有基因
uni_gene <- intersect(names(salmon_SRR7815790),names(ft_SRR7815790))
length(uni_gene)

plot(log1p(salmon_SRR7815790[uni_gene]),
     log1p(ft_SRR7815790[uni_gene]))
# 有些基因用不同的流程检测效果是不一样的，使用不同的参考数据得到的结果也是有区别的，但整体上一致











