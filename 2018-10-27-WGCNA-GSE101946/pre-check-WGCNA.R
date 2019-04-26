samples <- read.csv('sample_metadata.csv')
raw_counts <- read.csv('count_table.csv', row.names=1)
datTraits <- samples

# low-count filtering
low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
sprintf("Removing %d low-count genes (%d remaining).",sum(low_count_mask), sum(!low_count_mask))

filt_raw_counts <- raw_counts[!low_count_mask,]

# log-transforming data
log_counts <- log2(filt_raw_counts + 1)
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]



RNAseq_voom <- log_counts

# check WGCNA power
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad)[1:3300], decreasing = T),])
datExpr <- WGCNA_matrix

# 检查表达矩阵是否有太多的缺失值
library(WGCNA)
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK # 返回TRUE则继续【这里显示有一些基因没有通过】
# （可选）如果存在太多的缺失值
if (!gsg$allOK){
    # 把含有缺失值的基因或样本打印出来
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(colnames(datExpr)[!gsg$goodGenes], 
                                                  collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples],
                                                    collapse = ", ")));
    # 去掉那些缺失值
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 【图】画出软阈值分布:
if(T){
    pdf(file = "Soft threshold.pdf", width = 18, height = 10)
    par(mfrow = c(1,2));
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.85,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
}
#sft$powerEstimate

# net
if(T){
    net = blockwiseModules(
        datExpr,
        power = 4,
        maxBlockSize = 4000, #设置值大于datExpr基因数
        minModuleSize = 30, #每个模块基因数最少是30
        TOMType = "unsigned", #标准TOM模式
        reassignThreshold = 0, #p值为0的基因才可以重新在模块间分配
        mergeCutHeight = 0.25,#聚类的模块高度低于0.25的需要融合
        numericLabels = TRUE, pamRespectsDendro = FALSE,
        saveTOMs = FALSE,
        saveTOMFileBase = "Dm-genetype-TOM",
        verbose = 3
    )
}


# Convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
table(moduleColors)
# 【图】合并模块前的模块聚类图
if(T){
    pdf(file = "DendroAndColors.pdf", width = 18, height = 10)
    #plotDendroAndColors接受一个聚类的对象，以及该对象里面包含的所有个体所对应的颜色
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main= "Gene dendrogram and module colors")
    dev.off()
}

