########################
# 基因ID转换大法
########################
# no.1 =》mapIds转换
suppressMessages(library(org.Hs.eg.db))
symb <- mapIds(org.Hs.eg.db, keys=rownames(df), 
               keytype="ENSEMBL", column="SYMBOL")
length(rownames(df))
length(na.omit(symb))


# no.2 =》bitr转换
library(clusterProfiler)
library(org.Hs.eg.db)
symb_2 <- bitr(rownames(df), fromType = "ENSEMBL",
                toType = c("SYMBOL"),
                OrgDb = org.Hs.eg.db)
length(na.omit(symb_2$SYMBOL))

# 那些假装没有对应Symbol的Ensembl ID
not_mapped_ensembl <- names(symb[is.na(symb)])
length(not_mapped_ensembl)
head(not_mapped_ensembl,3)

# no.3 =》biomaRt转换
# 当前两者转换不全时，可以考虑用这种
library("biomaRt")
# 用useMart函数链接到人类的数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 除以以外还有很多：使用 listDatasets(ensembl) 查看

attributes <- listAttributes(ensembl)
View(attributes)

value <- not_mapped_ensembl
attr <- c("ensembl_gene_id","hgnc_symbol")
ids <- getBM(attributes = attr,
             filters = "ensembl_gene_id",
             values = value,
             mart = ensembl)







