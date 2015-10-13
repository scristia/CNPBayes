test_that("test_consensusCNP", {
    set.seed(100)
    starts <- rpois(1000, 100) + 10000000L
    ends <- rpois(1000, 100) + 10100000L
    cnv1 <- GRanges("chr1", IRanges(starts, ends))
    cnv1$id <- paste0("sample", seq_along(cnv1))
    starts <- rpois(500, 1000) + 101000000L
    ends <- rpois(500, 1000) + 101400000L
    cnv2 <- GRanges("chr5", IRanges(starts, ends))
    cnv2$id <- paste0("sample", seq_along(cnv2))
    cnv3 <- GRanges("chr1", IRanges(9000000L, 15000000L), id = "sample1400")
    starts <- seq(5000000L, 200000000L, 10000000L)
    ends <- starts + rpois(length(starts), 25000L)
    cnv4 <- GRanges("chr1", IRanges(starts, ends), id = paste0("sample",
        sample(1000:1500, length(starts))))
    all_cnvs <- suppressWarnings(c(cnv1, cnv2, cnv3, cnv4))
    grl <- split(all_cnvs, all_cnvs$id)
    cnps <- consensusCNP(grl)
    truth <- GRanges("chr1", IRanges(10000100L, 10100100L))
    seqinfo(truth) <- seqinfo(grl)
    expect_identical(truth, cnps)
    cnps <- consensusCNP(grl, max.width = 5e+05)
    truth <- GRanges(c("chr1", "chr5"), IRanges(c(10000100L,
        101000999L), c(10100100L, 101400999L)))
    seqlevels(truth, force = TRUE) <- seqlevels(grl)
    seqinfo(truth) <- seqinfo(grl)
    expect_identical(truth, cnps)
})
