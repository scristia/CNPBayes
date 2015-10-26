# return fraction that are not relabeled
# this works for marginal
# the batch method will work by running this method on each batch
relabeling <- function(thetas, comp_map) {
    comp_orig <- as.integer(names(comp_map))
    comp_merged <- names(comp_map) != comp_map

    ordered <- apply(thetas, 1, function(x) paste(order(x), collapse=""))

    for (i in comp_orig) {
        if (comp_merged[i]) {
            ordered <- gsub(i, comp_map[i], ordered)
        }
    }
    
    table_ordered <- sort(table(ordered), decreasing=TRUE)
    proportion <- table_ordered[1] / sum(table_ordered)

    return(proportion)
}
