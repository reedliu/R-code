## WGCNA处理转录组数据
# 思路：
#1.两个必须文件=》表达矩阵（rpkm或其他归一化表达量）
#2.看全局哪个模块颜色最深或者看要研究表型的哪个模块最深
#3.有了模块看基因，毕竟一个模块中包括了上百个基因，从中选hub genes
#4.有了核心基因（hub genes），GO、KEGG富集，网络分析

#############################
## step0: 获得表达矩阵datExpr
##############################
#wget -c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55839/suppl/GSE55839_RAW.tar
#tar -xf GSE48213_RAW.tar
#gzip -d *.gz
# 将该文件夹下所有文件一行行提取出来，并且每个文件都去掉EnsEMBL_Gene_ID这一行
# shell 命令：awk '{print FILENAME"\t"$0}' * |grep -v EnsEMBL_Gene_ID >tmp.txt
a=read.table('GSE76150_EdgeR_normalised_FB_counts.txt',sep = '\t',stringsAsFactors = F, header = T)
a <- a[,-c(2,3)]
rownames(a) <- a[,1]
a <- a[,-1]
raw_counts <- a
library(reshape2)
library(magrittr)
fpkm <- dcast(a,formula = V2~V1) #长数据框变宽数据框，指定V2是行，V1是列
rownames(fpkm)=fpkm[,1] #把第一列当作行名
fpkm=fpkm[,-1] # 去掉第一列
#只取列名前半部分（以_分割的），sapply得到的是列表，我们取第一个列表的第一个字符串
colnames(fpkm)=sapply(colnames(fpkm),function(x) strsplit(x,"_")[[1]][1])
RNAseq_voom <- fpkm[c(1:56)]#选择自己要研究的样本（这里全选）

#可以看到有36000多个基因，一次没必要分析这么多，选出5000个变化大的基因，但怎么选是个问题
#【用绝对中位差mad可以估计方差，相当于掐头去尾，不考虑极端值，
#用处于中间位置的数据与中位数的距离来反映原向量的数据波动情况】
# 因此可以降序排列，选前5000个mad值大的，也就是波动大的数据更有可能是差异基因
# 选出来之后，因为要对基因进行聚类，所以需要转置；（p.s.针对样本的聚类分析用hclust即可)
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr <- WGCNA_matrix  ## top 5000 mad genes
# 自己分析的数据，就是从总基因中提取DEGs来转换得到datExpr

# 检查表达矩阵是否有太多的缺失值
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK # 返回TRUE则继续
## 对样本进行聚类,检测异常值
if(F){
    datExprTree = hclust(dist(datExpr), method = "average")
    par(cex = 0.6)
    par(mar = c(0,4,2,0))
    plot(datExprTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
}
## 如果有，则去除异常值，得到过滤后的表达矩阵
if(F){
    abline(h = 15, col = "red") #画一条辅助线,h的值根据上图自定义
    # 把高于15的切除
    clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
    table(clust) # 0代表切除的，1代表保留的
    keepSamples = (clust==1)
    datExpr = datExpr0[keepSamples, ]
}

########################
## 获得表型矩阵datTraits
########################
library(GEOquery)
b=getGEO('GSE76150')
c <- pData(b[[1]])
metadata=pData(b[[1]])[,c(2,8)] 
metadata$characteristics_ch1.2 %>% unique() 
#可以看到这56个乳腺癌细胞系样本被分成了5组，即5个亚型

datTraits = data.frame(gsm=metadata[,1],
                       source=as.character(metadata$source_name_ch1))

# 将metadata的2、3列分别取出：后面的内容
# 其中trimws是trim white space，将取出的内容去除开头的空格

## 保证表型与样本名字对应
sampleNames = rownames(datExpr)
rownames(datTraits) = datTraits[sampleNames%in%datTraits$gsm, 1] 

save(fpkm,datTraits,datExpr,file = 'GSE48213-wgcna-input.RData')

######################正式开始####################################
## step1: 确定最佳β值 =》软阀值（soft thresholding power)
#################################################################
library(WGCNA)
load('GSE48213-wgcna-input.RData')
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
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
    abline(h=0.90,col="red")
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
    pdf(file = "checkBeta.pdf", width = 15, height = 10)
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

## 【图】合并模块后的模块聚类图
if(T){
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

# 【不需要运行】试着对样本进行聚类分析，结果可以看出：即使聚类关系
# 很近的细胞系，属于的样本也不同，因此仅仅根据分组进行差异分析，不全面
if(F){
    datExprTree<-hclust(dist(datExpr), method = "average")
    par(mar = c(0,5,2,0))
    plot(datExprTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
         cex.axis = 1, cex.main = 1,cex.lab=1)
    #没有表型的话，做到上面得到样本聚类图就好了；
    #有表型数据，可以加进来，看看是否聚类合理
    sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
                                    colors = c("white","blue","red","green","orange"),signed = FALSE)
    par(mar = c(1,4,3,1),cex=0.8)
    plotDendroAndColors(datExprTree, sample_colors,
                        groupLabels = colnames(sample_colors),
                        cex.dendroLabels = 0.8,
                        marAll = c(1, 4, 3, 1),
                        cex.rowText = 0.01,
                        main = "Sample dendrogram and trait heatmap")
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
design=model.matrix(~0+ datTraits$subtype) 
colnames(design)=levels(datTraits$subtype)
# 计算二者相关性
moduleTraitCor = cor(MEs, design , use = "p")#关于use的解释见cor {stats}
nSamples <- nrow(datExpr)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 【图】模块、表型热图 + 【表】模块表型相关性、p值
if(T){
    pdf(file = "labeledHeatmap.pdf", width = 10, height = 15)
    # 文字展示相关系数和p值
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3))
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
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    dev.off()
    options(digits=3)
    cor_module <- "Luminal" #选择要计算相关性的表型
    ModSmpCor <- cbind(moduleTraitCor[,cor_module],
                       moduleTraitPvalue[,cor_module])
    colnames(ModSmpCor) <- c("Correlation", "P-value")
    write.table(ModSmpCor,"Module-sample-cor.xls", sep = "\t", quote = F)
}
#每一种乳腺癌都有跟它强烈相关的模块

## 看各个模块
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
Luminal = as.data.frame(design[,3])
names(Luminal) = "Luminal"
MET = orderMEs(cbind(MEs, Luminal)) #MET=MEs+Trait
#【图】合并模块聚类、热图
if(T){
    pdf(file = "Eigengene-dengro-heatmap.pdf", width = 10, height = 15)
    par(cex = 0.9)
    plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,4,1,2), cex.lab = 0.8, xLabelsAngle
                          = 90)
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
    plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(4,5,2,2),
                          plotDendrograms = FALSE, xLabelsAngle = 90)
    dev.off()
}

##【图】各模块热图加相应样本柱状图
if(T){
    module="brown" #replace with module of interest
    pdf(file = paste0(module,"-Eigen-heatmap-ME-boxplot.pdf"), width = 10, height = 15)
    ME=MEs[, paste("ME",module, sep="")]
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
    plotMat(t(scale(datExpr[,moduleColors==module ])),
            nrgcols=30,rlabels=F,rcols=module,
            main=module, cex.main=2)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME, col=module, main="",  cex.names=1, cex.main=2,
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
    GS=as.numeric(cor(as.data.frame(design)$Luminal,datExpr, use="p"))
    GeneSignificance=abs(GS)
    ModuleSignificance=tapply(GeneSignificance,
                              moduleColors, mean, na.rm=T)
    pdf(file = "Gene-significance.pdf", width = 18, height = 10)
    par(mar=c(8.3,5,2,2))
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
    whichTrait="Luminal" #Replace this with the trait of interest
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
    luminal = as.data.frame(design[,3])
    names(luminal) = "luminal"
    # 得到矩阵
    geneTraitSignificance <-  as.data.frame(cor(datExpr, luminal, use = "p"))
    # 矩阵t检验
    GSPvalue <-  as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    # 修改列名
    names(geneTraitSignificance) = paste("GS.", names(luminal), sep="")
    names(GSPvalue) = paste("p.GS.", names(luminal), sep="")
}

## 【图】3.合并两个相关性矩阵,找到是否存在模块、表型都高度相关的基因
if(T){
    module <-  "brown"
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
module = "brown";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module)
modProbes = probes[inModule]
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# GO 分析：
if(T){
    ego.CC <- enrichGO(gene          = modProbes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.8,
                       qvalueCutoff  = 0.5)
    dotplot(ego.CC, showCategory=30)
    barplot(ego.CC)
    emapplot(ego.CC)
    plotGOgraph(ego.CC)
    
    # 保存数据
    GO_CC <- as.data.frame(ego.CC)
    write.table(GO_CC, paste0("GO.cc-module",module),sep = "\t", quote = F)
}
#KEGG分析
if(T){
    kk.genes <- bitr(modProbes, fromType = "ENSEMBL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
    kk <- enrichKEGG(gene   = kk.genes$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = 0.8,
                     qvalueCutoff = 0.8)
    head(kk)
    KEGG <- as.data.frame(kk)
    write.table(KEGG, paste0("KEGG-module",module),sep = "\t", quote = F)
    browseKEGG(kk, 'hsa04520')
}

#对于非模式物种
# 情况一：可以在AnnotationHub上在线抓取Org.Db的非模式生物
if(F){
    # 抓取 OrgDb
    require(AnnotationHub)
    hub <- AnnotationHub()
    query(hub, "Helicoverpa")
    Helicoverpa <- hub[['AH59375']]
    # 查看 OrgDb 里的Gene数量
    length(keys(Helicoverpa))
    # 查看 OrgDb 里的Gene ID类型
    columns(Helicoverpa)
    # 自己的数据需要转换ID，例如：
    require(clusterProfiler)
    bitr(keys(Helicoverpa)[1], 'ENTREZID', c("REFSEQ", "GO", "ONTOLOGY"), Helicoverpa)
    # ID转换成和org.db一致的样式后，就能进行enrichGO【这里直接用org.db中的部分基因】
    sample_genes <- keys(Helicoverpa)[1:10000]
    head(sample_genes)
    ego.CC <- enrichGO(gene        = sample_genes,
                       OrgDb         = Helicoverpa,
                       keyType       = 'ENTREZID',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1)
    dotplot(ego.CC, showCategory=30)
    barplot(ego.CC)
    
}
# 情况二：不能在线抓取Org.Db的非模式生物
# 1.下载注释文件
# 2.准备所要富集【基因列表】
# 3.基因与GO号的对应关系【term2gene】
# 4.GO号与其对应注释的文件【term2name】
if(F){
    # require将会根据包的存在与否返回true或者false
    require("clusterProfiler") 
    # 导入基因列表
    gene <- read.csv(file,header = T,sep=",")
    gene <- as.factor(gene$V1)
    # 导入注释文件
    term2gene <- read.csv(file,header=F,sep=",")
    term2name <- read.csv(file,header=F,sep=",")
    # 富集分析
    x <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
    # 设置结果输出文件
    ouf <- paste(out_file,sep ="\t")
    # 输出结果
    write.csv(x,ouf)
    # 绘制条形图
    barplot(x)
    # 绘制气泡图
    dotplot(x)
}
# 如果要做【多列基因富集对比】，
# 还是先准备基因列表（多列多个header对应多个表型）
if(F){
    # 导入文件
    gene <- read.csv(file,header = T,sep=",")
    # 富集分析
    x <- compareCluster(gene, fun='enricher',TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
    # 绘制气泡图
    dotplot(x, showCategory=10,includeAll=TRUE)
}


###############################
##step8: 找模块hub.genes
###############################
# 重新计算TOM
allowWGCNAThreads()
TOM = TOMsimilarityFromExpr(datExpr, power = 6)
# 选择模块
module = "brown"
# 指定模块中的基因
probes = colnames(datExpr) ## probe就是基因名
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
subGenes <- bitr(colnames(subTOM),fromType = "ENSEMBL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# 【图】30个hub基因的热图：
if(T){
    pdf("./plots/hubgenes-heatmap.pdf", width = 15, height = 10)
    par(mar=c(1,2,2,6))
    plotNetworkHeatmap(datExpr,
                       plotGenes = subGenes$ENSEMBL,
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


