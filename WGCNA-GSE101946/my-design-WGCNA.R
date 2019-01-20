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

# soft power
if(F){
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
}

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

# labeled Heatmap
if(T){
    # 先看模块
    moduleColors <- labels2colors(net$colors)
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0) #根据模块中基因与表型的相关性进行排序
    # 再看表型
    # 针对离散型的表型数据（这里的subtype），构建差异比较矩阵
    design=model.matrix(~0+ make.names(datTraits$tissue))
    colnames(design)=c("Control", "EcR-DN",
                       "E75-RNAi","Sox14-RNAi")
    
    # 计算二者相关性
    moduleTraitCor = cor(MEs, design , use = "p")#关于use的解释见cor {stats}
    nSamples <- nrow(datExpr)
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    # 【图】模块、表型热图 + 【表】模块表型相关性、p值
    if(T){
        pdf(file = "labeledHeatmap.pdf", width = 6, height = 4)
        # 文字展示相关系数和p值
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        par(mar = c(8, 8.5, 3, 3))
        par(cex = 0.6)
        # Display the correlation values within a heatmap plot
        labeledHeatmap(Matrix = moduleTraitCor,
                       xLabels = colnames(design),
                       yLabels = names(MEs),
                       ySymbols = names(MEs),
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = 0.5,
                       cex.lab.x = 0.8,
                       zlim = c(-1,1),
                       main = paste("Module-trait relationships"))
        dev.off()
        options(digits=3)
        cor_module <- "EcR-DN" #选择要计算相关性的表型
        ModSmpCor <- cbind(moduleTraitCor[,cor_module],
                           moduleTraitPvalue[,cor_module])
        colnames(ModSmpCor) <- c("Correlation", "P-value")
        write.table(ModSmpCor,"Module-sample-cor.xls", sep = "\t", quote = F)
    }
}

# merged module cluster
if(T){
    geneTree = net$dendrograms[[1]]
    library(flashClust)
    #Calculate dissimilarity of module eigenegenes
    MEDiss= 1-cor(MEs)
    #Cluster module eigengenes
    METree= flashClust(as.dist(MEDiss), method= "average")
    #plot(METree, main= "Clustering of module eigengenes", 
    #xlab= "", sub= "")
    #abline(h=MEDissThres, col="red")
    MEDissThres= 0.42 #相关性系数大于0.58的模块合并掉
    merge= mergeCloseModules(datExpr, moduleColors, 
                             cutHeight= MEDissThres, verbose =3)
    
    mergedColors= merge$colors
    #mergedMEs= merge$newMEs
    pdf(file = "Merged-DendroAndColors.pdf", width = 18, height = 10)
    plotDendroAndColors(geneTree, 
                        cbind(moduleColors, mergedColors), 
                        c("Dynamic Tree Cut", "Merged dynamic"), 
                        dendroLabels= FALSE, hang=0.03, 
                        addGuide= TRUE, guideHang=0.05
    )
    dev.off()
}

# 【图】各个模块基因显著性
if(T){
    GS=as.numeric(cor(as.data.frame(design)[,2],datExpr, use="p"))
    GeneSignificance=abs(GS)
    ModuleSignificance=tapply(GeneSignificance,
                              moduleColors, mean, na.rm=T)
    pdf(file = "Gene-significance.pdf", width = 18, height = 10)
    par(mar=c(11,5,2,2))
    plotModuleSignificance(GeneSignificance,moduleColors,
                           cex=1.5, cex.main=2, cex.lab=1.5, cex.axis=1.5,
                           cex.sub=1.5,axis.lty=1, las=3, ylim=c(0,0.7))
    dev.off()
    maxModule <- which.max(ModuleSignificance[names(ModuleSignificance != "grey")])
    cat("The most significant module is:", names(maxModule), maxModule)
}

# 某个模块的散点图
if(T){
    #接下来选某个模块看基因相关性，看散点图
    ## 1.计算模块与基因的相关性矩阵（MM: Module Membership）
    if(T){
        # 提取module（color）的名字
        modNames = substring(names(MEs), 3)
        # 得到矩阵
        geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
        # 矩阵t检验
        MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
        # 修改列名
        names(geneModuleMembership) = paste("MM", modNames, sep="")
        names(MMPvalue) = paste("p.MM", modNames, sep="")
    }
    
    ## 2.计算表型与基因的相关性矩阵（GS: Gene Significance）
    # 看模块热图决定=》这里以Luminal为例
    if(T){
        E75 = as.data.frame(design[,2])
        names(E75) = "E75"
        # 得到矩阵
        geneTraitSignificance <-  as.data.frame(cor(datExpr, E75, use = "p"))
        # 矩阵t检验
        GSPvalue <-  as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
        # 修改列名
        names(geneTraitSignificance) = paste("GS.", names(E75), sep="")
        names(GSPvalue) = paste("p.GS.", names(E75), sep="")
    }
    
    ## 【图】3.合并两个相关性矩阵,找到是否存在模块、表型都高度相关的基因
    if(T){
        module <-  "tan"
        pdf(file = paste0(module,"-MM-GS-scatterplot.pdf"), width = 10, height = 10)
        column <-  match(module, modNames) #找到目标模块所在列
        moduleGenes <-  moduleColors==module #找到模块基因所在行
        par(mfrow = c(1,1))
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = "Gene significance for body weight",
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        dev.off()
    }
}

# GO
module = "tan"
# GO-BP
if(T){
    # Select module probes
    probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
    inModule = (moduleColors==module)
    modProbes = probes[inModule]
    library(clusterProfiler)
    library(org.Dm.eg.db)
    library(ggplot2)
    #write.table(modProbes,"Dm-ID.txt", sep = "\t", quote = F)
    #transID <- read.csv("FlyBase_Fields_download.txt", sep = "\t")
    goID <- bitr(modProbes, 'SYMBOL', "ENTREZID", org.Dm.eg.db)
    ego.BP <- enrichGO(gene          = goID$ENTREZID,
                       OrgDb         = org.Dm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
}
dotplot(ego.BP, showCategory=10,title = "GO biological process (BP) enrichment in tan module")
barplot(ego.BP)
emapplot(ego.BP)
plotGOgraph(ego.BP)

GO_BP <- as.data.frame(ego.BP)
write.table(GO_BP, paste0("GO.BP-module-",module,".xls"),sep = "\t", quote = F)

#GO-CC
if(T){
    ego.CC <- enrichGO(gene          = goID$ENTREZID,
                       OrgDb         = org.Dm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
}
dotplot(ego.CC, showCategory=10, title = "GO cellular components (CC)enrichment in tan module")
barplot(ego.CC)
emapplot(ego.CC)
plotGOgraph(ego.CC)
GO_CC <- as.data.frame(ego.CC)
write.table(GO_CC, paste0("GO.CC-module-",module,".xls"),sep = "\t", quote = F)

#GO-MF
if(T){
    ego.MF <- enrichGO(gene          = goID$ENTREZID,
                       OrgDb         = org.Dm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
    dotplot(ego.MF, showCategory=10, title = "GO enrichment in tan module")
}

#KEGG
if(T){
    keggID <- bitr(modProbes, 'SYMBOL', c("ENTREZID", "UNIPROT"),org.Dm.eg.db)
    kk <- enrichKEGG(gene   = keggID$UNIPROT,
                     keyType = "uniprot",
                     organism     = 'dme',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
    head(kk)
    KEGG <- as.data.frame(kk)
    write.table(KEGG, paste0("KEGG-module",module),sep = "\t", quote = F)
    browseKEGG(kk, 'dme00534')
}


# 找hub基因
if(T){
    TOM = TOMsimilarityFromExpr(datExpr, power = 4)
    # 选择模块
    module = "tan"
    # 指定模块中的基因
    probes = colnames(as.data.frame(datExpr)) ## probe就是基因名
    inModule = (moduleColors==module)
    modProbes = probes[inModule]
    
    # 指定模块的TOM矩阵
    modTOM = TOM[inModule,inModule]
    dimnames(modTOM) = list(modProbes, modProbes)
    nTop = 30
    modConn = softConnectivity(datExpr[, modProbes])
    top = (rank(-modConn) <= nTop)
    
    # 转换hubgenes的基因名（=》entrez ID）
    subTOM <- modTOM[top, top]
    subGenes <- bitr(colnames(subTOM),fromType = "SYMBOL", toType = c("ENTREZID","UNIPROT"),
                     OrgDb = org.Dm.eg.db)
    write.table(subGenes,"subgenes.xls",sep = "\t")
    
}

## cytoscape
if(T){
    cyt = exportNetworkToCytoscape(
        modTOM,
        edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
        weighted = TRUE,
        threshold = 0.02,
        nodeNames = modProbes, 
        nodeAttr = moduleColors[inModule]
    )
    
}
