### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-11-08
### Email: jieandze1314@gmail.com
### Title: R包批量发邮件
### ---------------

# 一系列的报错及试错----------------------------------
#这个过程难在安装依赖包，需要java环境
install.packages('rJava',dependencies = T)
library(rJava)
# 出现报错
# Unable to find any JVMs matching version "(null)".
# No Java runtime present, try --request to install.
# Library not loaded: /Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk/Contents/Home/lib/server/libjvm.dylib

# 可能需要安装JDK11版本，目前mac上只有JDK 13，估计是版本太高
# JDK的版本命名：jdk-interim.update.patch.jdk
# JDK 历史版本在：https://jdk.java.net/archive/
# 下载的tar.gz文件解压后放在/Library/Java/JavaVirtualMachines目录
# 使用java --version查看版本

# 删除旧版本，比如之前安装过JDK 13，需要删掉
# sudo rm -rf jdk-13.0.1.jdk
# 然后在terminal输入：sudo R CMD javareconf
# 结果就能看到版本变化了：Java version: 11.0.1

# 参考了：https://blog.csdn.net/weixin_38986122/article/details/80931223
# 中的”在R中安装rJava及建立连接“

# 尝试安装----------------------------------
install.packages('mailR',dependencies = T)
library(mailR)
# 可见只要按照报错的要求操作，缺啥补啥，就没问题

# 尝试使用----------------------------------
sender <- "bioinfoplanet520@yeah.net"
recipients <- c("jieandze1314@gmail.com")
send.mail(from = sender,
          to = recipients,
          subject = "Program Done.",
          body = "My program is finished.",
          smtp = list(host.name = "smtp.yeah.net", port = 465,
                      user.name = "bioinfoplanet520@yeah.net",
                      passwd = "bioinfo520", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

# 继续解决报错----------------------------------
# 结果又发现报错，显示Error: NoClassDefFoundError (Java): javax/activation/DataHandler

# 搜索一次：https://github.com/rpremraj/mailR/issues/77
# "wush978 "说让把两个.jar文件放到下面👇的目录
system.file("java", package = "mailR")
# 这次还是报错，先检查下载的文件，结果发现下载的两个命名有问题，都是1.2.0这样的名称，而不是.jar后缀

# 搜索第二次：这次搜索这两个文件
# https://code.google.com/archive/p/javamail-android/downloads

# 这次的下载命名是对的，下载的路径也是对的
# 结果成功了："Java-Object{org.apache.commons.mail.SimpleEmail@6e509ffa}"







