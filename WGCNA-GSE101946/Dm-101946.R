# GSE101946-Drosophila WGCNA
#############################
## step0: 获得表达矩阵datExpr
##############################
a=read.table('./GSE101946_Supplemental_Table_5f.txt',
             sep = '\t',stringsAsFactors = F) 
colnames(a) <- a[1,]
a <- a[-1,]
rownames(a) <- a[,1]
a <- a[,-1]
# character to numeric
b <- as.data.frame(sapply(a, as.numeric))
fpkm <- b

library(magrittr)
RNAseq_voom <- fpkm[,c(1:43)]
rownames(RNAseq_voom) <- rownames(a)

WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad)[1:7500], decreasing = T),])
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

## 对样本进行聚类,检测异常值
if(T){
    sampleTree = hclust(dist(cor(datExpr)), method = "average")
    pdf(file = "pre-sampleClustering.pdf", width = 12, height = 9)
    par(cex = 0.6)
    par(mar = c(0,4,2,0))
    plot(sampleTree, main = "Sample clustering to detect outliers", 
         sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
    #abline(h = 42000, col = "red") #先画一条辅助线
    dev.off()
}

########################
## 获得表型矩阵datTraits
########################
library(GEOquery)
a=getGEO('GSE101946') 
b <- pData(a[[1]]) #看一下全部表型信息
metadata <- b[,c(2,13,12)] #选择2,8,13列作为表型信息原型
metadata[,3] %>% unique() 

datTraits = data.frame(gsm=metadata[,1],
                       state=trimws(sapply(as.character(metadata$characteristics_ch1.3),
                                           function(x) strsplit(x,":")[[1]][2])),
                       tissue=trimws(sapply(as.character(metadata$characteristics_ch1.2),
                                              function(x) strsplit(x,":")[[1]][2]))
)
datTraits <- datTraits[65:107,]
## 保证表型与样本名字对应
sampleNames = rownames(datExpr)
rownames(datTraits) = sampleNames




######################正式开始####################################
## step1: 确定最佳β值 =》软阀值（soft thresholding power)
#################################################################
library(WGCNA)
load('GSE111057-Dm-input.RData')
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
powers = c(c(1:10), seq(from = 12, to=16, by=2))
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
# 最佳的beta值就是sft$powerEstimate
# 【图】检验选定的β值下记忆网络是否逼近 scale free
# 【一般要求k与p(k)的相关性达到0.85时的power作为β值】
if(T){
    softPower <- sft$powerEstimate
    k <- softConnectivity(datE=datExpr,power=softPower) 
    pdf(file = "checkBeta.pdf", width = 10, height = 8)
    par(mfrow=c(1,2))
    hist(k)
    scaleFreePlot(k,main="Check Scale free topology\n")
    dev.off()
}

##########################################
##step2: 一步法构建共表达矩阵【关键一步】
##########################################
# 计算基因邻接性=》基因相似性=》基因相异性系数=》基因间的系统聚类树
# =》动态剪切法确定基因模块=》模块的特征向量值=》模块聚类=》合并距离较近的模块
if(T){
    net = blockwiseModules(
        datExpr,
        power = sft$powerEstimate,
        maxBlockSize = 6000, #设置值大于datExpr基因数
        minModuleSize = 30, #每个模块基因数最少是30
        TOMType = "unsigned", #标准TOM模式
        reassignThreshold = 0, #p值为0的基因才可以重新在模块间分配
        mergeCutHeight = 0.25,#聚类的模块高度低于0.25的需要融合
        numericLabels = TRUE, pamRespectsDendro = FALSE,
        saveTOMs = TRUE,
        saveTOMFileBase = "BRCA-cellline-TOM",
        verbose = 3
    )
}
table(net$colors) 

save(fpkm,RNAseq_voom,net,datTraits,datExpr,file = 'GSE101946-Dm-input.RData')
##########################################
##step3:众多基因=》模块：聚类
# 无法归类的基因为灰色
#【灰色过多表示前期筛选表达矩阵不合格】
##########################################
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
## 结果就是所有的基因被分到十几个模块中，并且模块做了聚类
## 【图】合并模块后的模块聚类图【MEs在第四部分】
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

#################################################
##step4:模块=》表型：热图（研究哪种表型就选哪个模块）
# 重点就是创建moduleTraitCor =》模块和表型相关矩阵
#################################################
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
    pdf(file = "labeledHeatmap.pdf", width = 5, height = 8)
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
    cor_module <- "Isolated_neurons_E75_RNAi" #选择要计算相关性的表型
    ModSmpCor <- cbind(moduleTraitCor[,cor_module],
                       moduleTraitPvalue[,cor_module])
    colnames(ModSmpCor) <- c("Correlation", "P-value")
    write.table(ModSmpCor,"Module-sample-cor.xls", sep = "\t", quote = F)
}


## 看各个模块
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
E75 = as.data.frame(design[,2])
names(E75) = "E75 RNAi"
MET = orderMEs(cbind(MEs, E75)) #MET=MEs+Trait
#【图】合并模块聚类、热图
if(T){
    pdf(file = "Eigengene-dengro-heatmap.pdf", width = 10, height = 10)
    par(cex = 1)
    plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                          marHeatmap = c(6,6,1,2), 
                          cex.lab = 0.8, xLabelsAngle = 90)
    dev.off()
}
#【图】分开模块聚类、热图
if(T){
    pdf(file = "Eigengene-dendrogram.pdf", width = 10, height = 10)
    par(cex = 1.0)
    plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                          plotHeatmaps = FALSE)
    dev.off()
    pdf(file = "Eigengene-heatmap.pdf", width = 10, height = 10)
    par(cex = 1.0)
    plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(6,6,2,2),
                          plotDendrograms = FALSE, xLabelsAngle = 90)
    dev.off()
}
##【图】各模块热图加相应样本柱状图
if(T){
    module="blue"#replace with module of interest
    pdf(file = paste0(module,"-Eigen-heatmap-ME-boxplot.pdf"), width = 7, height = 10)
    ME=MEs[, paste("ME",module, sep="")]
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
    plotMat(t(scale(datExpr[,moduleColors==module ])),
            nrgcols=30,rlabels=F,rcols=module,
            main=module, cex.main=2)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME, col=module, main="",  cex.names=1, cex.main=2,
            #legend.text = sapply(rownames(datExpr) ,function(x) gsub("FPKM_", "", x)),
            ylab="eigengene expression",xlab="sample")
    dev.off()
    cat(paste0("The sample of max eigengene in module", module,
               "is:",rownames(datExpr)[which.max(ME)],"\n"))
    cat(paste0("The sample of min eigengene in module", module,
               "is:",rownames(datExpr)[which.min(ME)]))
}

#################################################
##step5:选出来模块=》模块内基因
# 同时看模块、表型与基因的关系
#################################################
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
# 【图】组合各个模块基因相关性的散点图
#merged前的情况
if(T){
    whichTrait="Isolated.γ.neurons.E75.RNAi" #Replace this with the trait of interest
    nGenes = ncol(datExpr);
    nSamples = nrow(datExpr);
    selTrait = as.data.frame(design[,whichTrait]);
    names(selTrait) = whichTrait
    modNames = substring(names(MEs), 3)
    geneModuleMembership = as.data.frame(signedKME(datExpr, MEs));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    names(geneModuleMembership) = paste("MM", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    geneTraitSignificance = as.data.frame(cor(datExpr, selTrait, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
    names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
    pdf(file = "All-GS-MM-scatter.pdf", width = 20, height = 20)
    par(mfrow=c(6,5))
    for(module in modNames[1:length(modNames)]){
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste(module,"module membership"),
                           ylab = paste("GS for", whichTrait),
                           col = module,mgp=c(2.3,1,0))
    }
    dev.off()
}
#merged后的情况
if(T){
    mergedMEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
    mergedMEs = orderMEs(mergedMEs0)
    whichTrait="Luminal" #Replace this with the trait of interest
    nGenes = ncol(datExpr);
    nSamples = nrow(datExpr);
    selTrait = as.data.frame(design[,whichTrait]);
    names(selTrait) = whichTrait
    modNames = substring(names(mergedMEs), 3)
    geneModuleMembership = as.data.frame(signedKME(datExpr, mergedMEs));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    names(geneModuleMembership) = paste("MM", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    geneTraitSignificance = as.data.frame(cor(datExpr, selTrait, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
    names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
    pdf(file = "Merged-All-GS-MM-scatter.pdf", width = 20, height = 20)
    par(mfrow=c(4,5))
    for(module in modNames[1:length(modNames)]){
        column = match(module, modNames);
        moduleGenes = mergedColors==module;
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste(module,"module membership"),
                           ylab = paste("GS for", whichTrait),
                           col = module,mgp=c(2.3,1,0))
    }
    dev.off()
}

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
    E75 = as.data.frame(design[,3])
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
# 看结果的p值和cor值，进一步佐证模块中基因的价值


#################################################
##step6: 全部/随机部分基因热图
#################################################
## 【图】选择400个基因画图
if(T){
    nSelect = 400
    set.seed(10)
    dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6) #基因之间的相异度
    #上面一步法等价于：（1）得到beta值：softPower <- sft$powerEstimate
    #（2）邻接矩阵：adjacency = adjacency(multiExpr[[1]]$data, power = softPower)
    #（3）邻接转TOM矩阵：TOM = TOMsimilarity(adjacency)
    #（4）dissTOM = 1-TOM
    nGenes <- ncol(datExpr)
    select = sample(nGenes, size = nSelect)
    selectTOM = dissTOM[select, select]
    # 再计算基因之间的距离树(对于基因的子集，需要重新聚类)
    selectTree = hclust(as.dist(selectTOM), method = "average")
    selectColors = moduleColors[select]
    
    pdf(file = paste0("Sub400-netheatmap.pdf"), width = 10, height = 10)
    plotDiss = selectTOM^7
    diag(plotDiss) = NA #将对角线设成NA，在图形中显示为白色的点，更清晰显示趋势
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
    dev.off()
}
## 【图】选择全部基因画图(耗时较久，生成的文件很大)
# 结果可以看到各个区块的颜色差异,沿着对角线的深色区块就是模块Module
if(T){
    pdf(file = "All-netheatmap.pdf", width = 20, height = 20)
    plotTOM = dissTOM^7
    diag(plotTOM) = NA
    geneTree = net$dendrograms[[1]]
    TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
    dev.off()
}

###############################
##step7: 指定模块GO、KEGG
###############################
# Select module
module = "black"
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module)
modProbes = probes[inModule]
library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)
#write.table(modProbes,"Dm-ID.txt", sep = "\t", quote = F)
#transID <- read.csv("FlyBase_Fields_download.txt", sep = "\t")
goID <- bitr(rownames(RNAseq_voom), 'SYMBOL', "ENTREZID", org.Dm.eg.db)
# GO 分析：
if(T){
    ego.BP <- enrichGO(gene          = goID$ENTREZID,
                       OrgDb         = org.Dm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1)
    dotplot(ego.BP, showCategory=30)
    barplot(ego.BP)
    emapplot(ego.BP)
    plotGOgraph(ego.BP)
    
    # 保存数据
    GO_BP <- as.data.frame(ego.BP)
    write.table(GO_BP, paste0("GO.BP-module",module,".xls"),sep = "\t", quote = F)
}
#KEGG分析
if(T){
    keggID <- bitr(rownames(RNAseq_voom), 'SYMBOL', c("ENTREZID", "UNIPROT"),org.Dm.eg.db)
    kk <- enrichKEGG(gene   = keggID,
                     keyType = "uniprot",
                     organism     = 'dme',
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
    head(kk)
    KEGG <- as.data.frame(kk)
    write.table(KEGG, paste0("KEGG-module",module),sep = "\t", quote = F)
    browseKEGG(kk, 'hsa04520')
}


###############################
##step8: 找模块hub.genes
###############################
# 重新计算TOM
allowWGCNAThreads()
TOM = TOMsimilarityFromExpr(datExpr, power = 6)
# 选择模块
module = "black"
# 指定模块中的基因
probes = colnames(as.data.frame(datExpr)) ## probe就是基因名
inModule = (moduleColors==module)
modProbes = probes[inModule]

# 指定模块的TOM矩阵
modTOM = TOM[inModule,inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# 选出模块中连接度最高的基因【这里以30为例】
nTop = 30
modConn = softConnectivity(datExpr[, modProbes])
top = (rank(-modConn) <= nTop)

# 转换hubgenes的基因名（=》entrez ID）
subTOM <- modTOM[top, top]
subGenes <- bitr(colnames(subTOM),fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Dm.eg.db)

# 【图】30个hub基因的热图：
if(T){
    pdf("hubgenes-heatmap.pdf", width = 15, height = 15)
    par(mar=c(10,2,2,6))
    plotNetworkHeatmap(data.frame(datExpr),
                       plotGenes = subGenes$ENTREZID,
                       networkType = "unsigned",
                       useTOM = TRUE,
                       power=6,
                       main=paste0("Top",nTop,"  ","hubgenes","  ","correlations"))
    dev.off()
}

###############################
##step9: 导出模块hub.genes
###############################
#可以选择导出全部modTOM或者部分subTOM
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
## VisANT
if(T){
    vis = exportNetworkToVisANT(modTOM,
                                file = paste("VisANTInput-", module, ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0)
}






