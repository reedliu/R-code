### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-11-06
### Email: jieandze1314@gmail.com
### R获取芯片平台与注释
### ---------------

## 配置几个包
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
BiocManager::install("GEOmetadb")
install.packages("RSQLite")
devtools::install_github("ggrothendieck/sqldf")
library(GEOmetadb)
library(RSQLite)
library(sqldf)

## 使用SQLite
getSQLiteFile() 
# 这个函数也是直接下载到当前目录下
# trying URL 'http://starbuck1.s3.amazonaws.com/sradb/GEOmetadb.sqlite.gz'
# Content type 'binary/octet-stream' length 528141472 bytes (503.7 MB)
# 网速可以的话，先自己下载到本地

################################
## 如果你的下载出现了问题，直接跳过上面的操作
################################
# 我也把这个文件放在微云：链接：https://share.weiyun.com/5njbq8z 
# 密码：3b8pkg
# 然后解压缩GEOmetadb.sqlite.gz，解压缩后是8G
con <- dbConnect(RSQLite::SQLite(),'GEOmetadb.sqlite')
gplToBioc <- dbGetQuery(con,'select gpl,bioc_package,title from gpl where bioc_package is not null')






