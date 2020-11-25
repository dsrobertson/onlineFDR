## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

library(onlineFDR)
library(Rcpp)

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

set.seed(1)

## -----------------------------------------------------------------------------
set.seed(1); results.indep <- LOND(sample.df)    # for independent p-values
set.seed(1); results.dep <- LOND(sample.df, dep=TRUE)   # for dependent p-values

# compare adjusted significance thresholds
cbind(independent = results.indep$alphai, dependent = results.dep$alphai)


## -----------------------------------------------------------------------------
set.seed(1); results.LORD.plus <- LORD(sample.df)
set.seed(1); results.LORD3 <- LORD(sample.df, version=3)
set.seed(1); results.LORD.discard <- LORD(sample.df, version='discard')
set.seed(1); results.LORD.dep <- LORD(sample.df, version='dep') 

# compare adjusted significance thresholds
cbind(LORD.plus = results.LORD.plus$alphai,
    LORD3 = results.LORD3$alphai,
    LORD.discard  = results.LORD.discard$alphai,
    LORD.dep = results.LORD.dep$alphai)


