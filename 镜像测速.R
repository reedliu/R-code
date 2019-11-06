remotes::install_github("GuangchuangYu/mirrorselect")
library(mirrorselect)
# 测速CRAN
n <- get_mirror()
cn_cran <- n[grepl("cn",n)]
x <- mirrorselect(cn_cran)
knitr::kable(x)
# 测速Bioconductor
m  <- get_mirror('BioC')
cn_bioc <- m[grepl("cn",m)]
y <- mirrorselect(cn_bioc)
knitr::kable(y)