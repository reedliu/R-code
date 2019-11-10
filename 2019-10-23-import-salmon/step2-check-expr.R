### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: 检查得到的表达矩阵
### ---------------
rm(list = ls())
options(stringsAsFactors = F)
load('salmon-aggr.Rdata')
dim(salmon_expr)
colnames(salmon_expr)
salmon_expr[1:3,1:3]
library(pheatmap)
# 选全部基因
pheatmap(cor(salmon_expr),
         filename = 'salmon_cor_raw.png')
pheatmap(cor(log2(edgeR::cpm(salmon_expr)+1)),
         filename = 'salmon_cor_cpm.png')
dev.off()

# 根据在多少细胞中表达量大于1，选部分基因
pt_salmon_expr <- salmon_expr[apply(salmon_expr,1,
                                    function(x) sum(x>1)>10),]
dim(pt_salmon_expr)
pheatmap(cor(log2(edgeR::cpm(pt_salmon_expr)+1)))

# 根据mad选择前1000
pt_salmon_expr <- log2(edgeR::cpm(salmon_expr)+1)
pt_salmon_expr <- pt_salmon_expr[names(head(sort(
  apply(pt_salmon_expr, 1, mad),decreasing = T),1000)),]
dim(pt_salmon_expr)
pheatmap(cor(pt_salmon_expr))
dev.off()








