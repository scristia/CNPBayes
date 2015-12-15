label_components <- function(model) {
    all <- readRDS("~/Code/CNPBayesPaper/output/mb.RDS")
    model <- all$model[["CNP candidate 10"]]

    # default to all four copy number gain/deletion
    labels <- c("Homozygous Deletion",
                "Hemizygous Deletion",
                "Normal",
                "Copy Number Gain")

    # get overall means of components
    means <- mu(model)

    # get overall variances of components
    # N.B. that this accessor may need to change
    vars <- colMeans(sigma2(model))

    return(labels)
}
