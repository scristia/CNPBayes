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

    # assume no copy number gain/deletion if k==1 (valid assumption?)
    if (k == 1) {
        return(labels[3])
    }

    # get overall means of components
    means <- mu(model)

    # set unknown copy number and sort
    names(means) <- rep(NA, k)
    means <- sort(means)

    # see if there is a normal and hemizygous deletion
    normal <- min(means[means >= 0])
    hemi <- max(means[means < 0])

    if (hemi > -1 & normal < 0.5) {
        names(means[means == normal]) <- labels[3]
        names(means[means == hemi]) <- labels[2]
    }


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
