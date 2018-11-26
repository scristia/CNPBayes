library(tidyverse)
library(PancCnvs)
hapmap.dir <- "PancCnvsData2/inst/extdata/hapmap"
hap.list <- readRDS(file.path(hapmap.dir, "cnp_57.rds"))
experiment <- readRDS(file.path(hapmap.dir, "simulation_params.rds")) %>%
  mutate(id=as.numeric(id)) %>%
  filter(id == 156)
full.data <- hapmapData(experiment[1, ], hap.list)
saveRDS(full.data,
        file=file.path("CNPBayes/inst/extdata", "hapmap.rds"))
if(FALSE){
  ## batch 3 has no homozygous deletions
  A <- ggplot(full.data, aes(medians, ..density..)) +
         geom_histogram(bins=100,
                        color="gray70",
                        fill="gray70",
                        alpha=0.1) +
    geom_density(adjust=1, alpha=0.4, size=0.75, color="gray30") +
    facet_wrap(~batch)
}

