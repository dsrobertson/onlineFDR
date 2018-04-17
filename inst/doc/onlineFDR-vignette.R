## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ------------------------------------------------------------------------
library(onlineFDR)

sample.df <- data.frame(
    id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
        'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
        'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
    date = as.Date(c(rep("2014-12-01",3),
                    rep("2015-09-21",5),
                    rep("2016-05-19",2),
                    "2016-11-12",
                    rep("2017-03-27",4))),
    pval = c(2.90e-14, 0.06743, 0.01514, 0.08174, 0.00171,
            3.61e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
            0.69274, 0.30443, 0.000487, 0.72342, 0.54757))

results <- LORD(sample.df)
results


## ------------------------------------------------------------------------
results.indep <- LOND(sample.df)    # for independent p-values
results.dep <- LOND(sample.df, dep=TRUE)   # for dependent p-values

# compare adjusted test levels
cbind(independent = results.indep$alphai, dependent = results.dep$alphai)


## ------------------------------------------------------------------------
results.LORD1 <- LORD(sample.df, version=1)
results.LORD2 <- LORD(sample.df, version=2)
results.LORD3 <- LORD(sample.df)   # default version

# compare adjusted test levels
cbind(LORD1 = results.LORD1$alphai,
    LORD2 = results.LORD2$alphai,
    LORD3 = results.LORD3$alphai)


