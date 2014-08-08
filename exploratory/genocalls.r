library(snpStats)

path <- "/thumper/ctsa/ARIC_active/ForRob/White/W_ARIC_gwas_frz3v2_chr"

for(i in 1:23) {

    pedfile <- read.pedfile(file=paste0(path, i, ".ped"),
                            snps=paste0(path, i, ".map"), which=1, split=" ")

    genodata <- pedfile$genotypes@.Data
    colnames(genodata) <- pedfile$map[,2]

    saveRDS(genodata, file=paste0(dirname(path),"/W_genotypes_chr",i,".rds"))

}
