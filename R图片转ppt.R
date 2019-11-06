if(!require(export))install.packages("export")

library(export)
plot(iris$Sepal.Length)
graph2ppt(file="tryplot.pptx", width=7, height=5)
