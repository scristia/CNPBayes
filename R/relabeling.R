# return fraction that are not relabeled
# this works for marginal
# the batch method will work by running this method on each batch
relabeling <- function(model, merge=TRUE) {
    # handle case of merge == FALSE
    if (merge) {
        merged <- DensityModel(model, merge=TRUE)
        comp_map <- clusters(merged)
    }

    comp_orig <- as.integer(names(comp_map))
    comp_merged <- names(comp_map) != comp_map

    thetas <- theta(chains(model))
    ordered <- apply(thetas, 1, function(x) paste(order(x), collapse=""))



    for (i in comp_orig) {
        if (comp_merged[i]) {
            ordered <- gsub(as.character(i), as.character(comp_map[i]), ordered)
        }
    }
    `
    table_ordered <- sort(table(ordered), decreasing=TRUE)
    proportion <- table_ordered[1] / sum(table_ordered)

    return(proportion)
}
