library(CNPBayes)
library(CNPAric)
lrr <- readRDS("vi_lrr.rds")
vigr <- readRDS("VI-granges.rds")
data(se_aric)

# 63772765, 42909764, 8755522


rmat <- lrr[[223]][,colnames(aricCNPs)]
rmat <- rmat[-which(is.na(rowSums(rmat))),, drop=FALSE]
y = prcomp(t(rmat), center=TRUE, scale.=FALSE)
y <- y$x[,1]
hist(y, breaks=200, col="gray")
hist(colMeans(rmat)/100, breaks=200, col="gray", border="gray")
