### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-02-25
### Email: yunzeliu@um.edu.mo
### Title: 多张图组合后修改字体+导出为PPT
### ---------------

rm(list=ls())
options(stringsAsFactors = F)
load('grid_plots.Rdata')

#原来每张图是这样
p1=p_mb231_up
p2=p_mb231_down
p3=p_mcf7_up
p4=p_mcf7_down

# 用patchwork可以组合在一起，并且大小比例已经调整好
library(patchwork)
pp = p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
pp


## 之后分别调整theme字体（设成6号大小），再组合
# 先设置好setting，后面的可以随便做了
# 设置element_blank()可以不显示（比如Y轴title）
setting=element_text(size = 6,color="black",family="Arial")

p1x=p1+theme(axis.title.x = setting,
             axis.title.y = element_blank())+
    theme(axis.text.x = setting,
          axis.text.y = setting,
          legend.text = setting,
          legend.title = setting)

p2x=p2+theme(axis.title.x = setting,
             axis.title.y = element_blank())+
    theme(axis.text.x = setting,
          axis.text.y = setting,
          legend.text = setting,
          legend.title = setting)

p3x=p3+theme(axis.title.x = setting,
             axis.title.y = element_blank())+
    theme(axis.text.x = setting,
          axis.text.y = setting,
          legend.text = setting,
          legend.title = setting)

p4x=p4+theme(axis.title.x = setting,
             axis.title.y = element_blank())+
    theme(axis.text.x = setting,
          axis.text.y = setting,
          legend.text = setting,
          legend.title = setting)

# 排两列，加a-d标记，导出为PPT
pp2x = p1x + p2x + p3x + p4x + plot_layout(ncol = 2) + 
    plot_annotation(tag_levels = 'a')
# pp2x
# 这个长宽直接影响到结果的图的比例
graph2ppt(pp2x,"pp2x.pptx",width = 10, height = 8, orient="portrait")


