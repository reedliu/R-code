---
title: "Limma-voom"
author: "Reedliu"
date: "2019-03-11"
---
#############################
## 配置信息 
#############################
library(edgeR)
counts <- read.delim("all_counts.txt", row.names = 1)
head(counts[1:3,1:3])
dim(counts)
# 构建DGEList对象，将counts和sample信息包含进去
d0 <- DGEList(counts)

#############################
## 预处理
#############################
# 计算标准化因子
d0 <- calcNormFactors(d0)
d0
# 过滤低表达基因
cutoff <- 1
cut <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-cut,] 
dim(d) # 剩余基因

# 根据列名提取样本信息(sample name)
spname <- colnames(counts) 
spname
# 分离出分组信息
strain <- substr(spname, 1, nchar(spname) - 2) 
time <- substr(spname, nchar(spname) - 1, nchar(spname) - 1)
strain
time
# 再将这两部分整合进group分组信息中
group <- interaction(strain, time)
group

# Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(group))

suppressMessages(library(RColorBrewer))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1") 
col.group <- as.character(col.group)
plotMDS(d, labels=group, col=col.group) 
title(main="A. Sample groups")


#############################
## Voom转换及方差权重计算
#############################
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
# 不好的结果如下（主要看x轴0附近有没有迅速上升，有的话表示没有
# 过滤low counts）
tmp <- voom(d0, mm, plot = T)

#############################
## Voom转换及方差权重计算
#############################
fit <- lmFit(y, mm)
head(coef(fit),3)

contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# 找差异基因
top.table <- topTable(tmp, sort.by = "P", n = Inf)
DEG <- na.omit(top.table)
head(DEG, 5)
length(which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC)>2 ))

# 换个组比较
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
DEG <- na.omit(top.table)
head(DEG, 5)
length(which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC)>2 ))

#############################
## P.S. 换个思路
#############################
# 构建新的model matrix
mm <- model.matrix(~strain*time)
colnames(mm)

y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
head(coef(fit),3)

# 想比较品系I5和品系C在time6的差异
tmp <- contrasts.fit(fit, coef = 2) 
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
DEG <- na.omit(top.table)
head(DEG, 5)
length(which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC)>2 ))

# 想看time9与品系I5的差异结果
head(coef(fit),3)
# cultivarI5:time9
tmp <- contrasts.fit(fit, coef = 5) 
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
DEG <- na.omit(top.table)
#head(DEG, 5)
length(which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC)>2 ))


#############################
## 更为复杂的模型
#############################
# 考虑批次效应
batch <- factor(rep(rep(1:2, each = 2), 6))
batch

mm <- model.matrix(~0 + group + batch)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
DEG <- na.omit(top.table)
#head(DEG, 5)
length(which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC)>2 ))

# 有一个连续型变量rate，它可能是pH、光照等等对研究材料的影响值
# Generate example rate data[行数要与count矩阵的列数相等]
set.seed(10)
rate <- rnorm(n = 24, mean = 5, sd = 1.7)
rate
# 指定矩阵模型
mm <- model.matrix(~rate)
head(mm)

y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
tmp <- contrasts.fit(fit, coef = 2) # test "rate" coefficient
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
DEG <- na.omit(top.table)
#head(DEG, 5)
length(which(DEG$adj.P.Val < 0.05 & abs(DEG$logFC)>2 ))

AT1G01060 <- y$E["AT1G01060",]
plot(AT1G01060 ~ rate, ylim = c(6, 12))
intercept <- coef(fit)["AT1G01060", "(Intercept)"]
slope <- coef(fit)["AT1G01060", "rate"]
abline(a = intercept, b = slope)






