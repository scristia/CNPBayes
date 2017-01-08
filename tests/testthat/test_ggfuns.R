context("gg functions")

test_that("ggfun", {
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  ggSingleBatch(truth)


  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.02, 0.02, 0.02))
  ## hard to see marginal when it looks the same as the truth
  ggSingleBatch(truth)


  model <- truth
  colors <- c("#999999", "#56B4E9", "#E69F00", "#0072B2",
              "#D55E00", "#CC79A7",  "#009E73")
  df <- singleBatchDensities(model)
  df2 <- df
  df2$d <- 0
  df3 <- rbind(df, df2)



  tmp$y <- 0

  ggplot(df3, aes(x, d)) +
    geom_histogram(data=df.observed,
                   aes(y, ..density..), bins=bins,
                   inherit.aes=FALSE) +
    geom_polygon(aes(fill=name, group=name), alpha=0.4)

    stat_function(fun=dnorm,
                  args=list(mean=theta(model)[1],
                            sd=sigma(model)[1]),
                  aes(fill="blue"))

    geom_area(stat="identity", aes(color=name, fill=name),
              alpha=0.4) +
    xlab("quantiles") + ylab("density") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    guides(fill=guide_legend(""), color=guide_legend(""))


})
