### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-02-25
### Email: yunzeliu@um.edu.mo
### Title: GO result from Web via R to PPT
### ---------------

# 描述：一般从GeneOntology上传基因名，然后会用panther计算，结果是一个Excel表格（如果有多个样本，可以汇总到一个Excel中的不同sheet，最后批量操作、作图，最后导出为PPT）

rm(list=ls())
options(stringsAsFactors = F)

# install.packages("rio")
library(rio)
library(stringr)
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
library(dplyr)
library(patchwork)

####################################
# 首先导入Excel的各个sheet为数据框
####################################
# 多个sheet用import_list，读进来是个列表
tmp <- import_list("go 1.5 fc graph.xlsx")
# names(tmp)

# mcf7_up <- tmp$`mcf7 goup fc  1.5`
# mcf7_dn <- tmp$`mcf7 godown fc  1.5`
# mb231_up <- tmp$`231 go UP fc  1.5`
# mb231_dn <- tmp$`231 godown fc  1.5`

p <- list()
ppt <- list()
for(i in 1:length(tmp)){
    dat = tmp[[i]]
    dat = dat[,1:8]
    # 然后更改列名
    colnames(dat)[1] <- c("Description")
    colnames(dat)[length(colnames(dat))-1] <- c("pvalue")
    colnames(dat)[length(colnames(dat))] <- c("P.adjust")
    colnames(dat)
    dat$fold <- dat$`upload_1 (fold Enrichment)`
    dat <- dat[,-5]
    # 添加GO id列：先提取Description，然后以（）为分界线，提取其中的GO ID
    dat$ID <- str_split(str_split(dat[,1],"\\(",simplify = T)[,2],
                       "\\)",simplify = T)[,1]
    # 添加Count列
    dat$Count <- dat[,3]
    # 重新排序
    dat <- dat[,c(9,1:7,10,8)]
    dat <- dat[,-6]
    # 把dat的Description首字母大写
    dat$Description <- paste(toupper(substring(dat$Description, 1,1)), 
                             substring(dat$Description, 2), sep="")
    ####################################
    # 然后把数据框处理成作图的格式
    ####################################
    test <- dat
    rownames(test) <- test$ID
    test=arrange(test,fold)
    test$Description = factor(test$Description,levels = test$Description,ordered = T)
    
    ####################################
    # 作图
    ####################################
    library(ggplot2)
    p[[i]] <- ggplot(test,aes(x = fold,y = Description))+
        geom_point(aes(color = P.adjust,
                       size = Count))+
        scale_color_continuous(low="red", high="blue", name = 'P.adjust', guide=guide_colorbar(reverse=TRUE),
                               labels = function(x)format(x,scientific = F))+
        xlab("Fold Enrichment")+
        xlim(0,ceiling(max(test$fold)))+
        theme_bw()
    
    ####################################
    # 组合
    ####################################
    setting=element_text(size = 6,color="black",family="Arial")
    
    ppt[[i]]=p[[i]]+theme(axis.title.x = setting,
                 axis.title.y = element_blank())+
        theme(axis.text.x = setting,
              axis.text.y = setting,
              legend.text = setting,
              legend.title = setting)
}

####################################
# 排两列，加a-d标记，导出为PPT
####################################
final_ppt = ppt[[1]] + ppt[[2]] + ppt[[3]] + ppt[[4]] + plot_layout(ncol = 2) + 
    plot_annotation(tag_levels = 'a')

# Final_ppt
# 这个长宽直接影响到结果的图的比例
# orient设置竖向
graph2ppt(final_ppt,"final_ppt.pptx",width = 10, height = 8, orient="portrait")




















