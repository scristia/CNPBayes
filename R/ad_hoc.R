label_components <- function(model) {
    all <- readRDS("~/Code/CNPBayesPaper/output/mb.RDS")
    model <- all$model[["CNP candidate 10"]]

    # default to all four copy number gain/deletion
    labels <- c("Homozygous Deletion",
                "Hemizygous Deletion",
                "Normal",
                "Copy Number Gain")

    # get data and components
    y <- y(model)
    z <- z(model)
    k <- k(model)

    # get overall means of components
    means <- mu(model)

    # get overall variances of components
    # N.B. that this accessor may need to change
    vars <- colMeans(sigma2(model))

    # define "big variance"
    big.var <- 0.5

    # if the overall variance for a component is bigger than "big"
    # then note the component as outlier
    outlier.component <- var[var > big.var]

    # if there are any outlier components then do something about it

    return(labels)
}
