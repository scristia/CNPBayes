lrr.list <- readRDS("~scristia/MixtureModels/data/LRRList_adj.rds")

pdf("~/MixtureModels/snp-cn_hists/snp-cn_hist%03d.pdf", onefile=FALSE)
for(i in 1:length(lrr.list)) {
    cn.ind <- grep("CN", rownames(lrr.list[[i]]))
    snp.ind <- grep("SNP", rownames(lrr.list[[i]]))
    if(length(cn.ind) >= 10 & length(snp.ind) >= 10) {
        cn.means <- colMeans(lrr.list[[i]][cn.ind,], na.rm=TRUE)/100
        snp.means <- colMeans(lrr.list[[i]][snp.ind,], na.rm=TRUE)/100
        hist(snp.means, breaks=200, col=rgb(0.7,0.7,0.7,0.4), border=rgb(0.7,0.7,0.7,0.4),main=print(paste("Region", i)), freq=FALSE)
        hist(cn.means, breaks=200, col=rgb(0,0.8,1,0.4), border=rgb(0,0.8,1,0.4), freq=FALSE, add=TRUE)
        legend("topleft", legend=c("CN", "SNP"), col=c(rgb(0,0.8,1), rgb(0.8, 0.8, 0.8)), lty=1, bty="n")
    }
}
dev.off()


    cn.ind <- grep("CN", rownames(lrr.list0[[7]]))
    snp.ind <- grep("SNP", rownames(lrr.list0[[7]]))
    pdf("snpvscn.pdf")
    hist(colMeans(lrr.list0[[3]][snp.ind,], na.rm=TRUE)/100, breaks=200, col=rgb(0.8,0.8,0.8,0.4), border=rgb(0.8,0.8,0.8,0.5),main="Region 3", freq=FALSE)
    hist(colMeans(lrr.list0[[3]][cn.ind,], na.rm=TRUE)/100, breaks=200, col=rgb(0,0.8,1,0.4), border=rgb(0,0.5,1,0.5), freq=FALSE, add=TRUE)
    legend("topleft", legend=c("CN", "SNP"), col=c(rgb(0,0.8,1,0.8), rgb(0.8, 0.8, 0.8)), lty=1, bty="n")
    dev.off()
