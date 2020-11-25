## ----setup, include = FALSE---------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

knitr::include_graphics("stream-diagram.png")

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
library(onlineFDR)
library(Rcpp)

set.seed(1)
LOND_results <- LOND(sample.df)
LOND_results

## -----------------------------------------------------------------------------
LOND_results %>%
  filter(R == 1) %>%
  nrow()

## -----------------------------------------------------------------------------
set.seed(1)
LORD_results <- LORD(sample.df)

cbind(LOND_results, LORD_results$alphai) %>%
      mutate(index = row_number(),
             LOND = log(alphai),
             LORD = log(LORD_results$alphai),
             Bonferroni = log(0.05/index),
             Unadjusted = rep(log(0.05), nrow(.))) %>%
      pivot_longer(cols = c(LOND, LORD, Bonferroni, Unadjusted),
                   names_to = "adjustment",
                   values_to = "alpha") %>%
  ggplot(aes(x = index, y = alpha, col = adjustment)) + 
    geom_line()

## -----------------------------------------------------------------------------
sample.df <- data.frame(
    id = c('A15432', 'B90969', 'C18705'),
    date = as.Date(c(rep("2014-12-01",3))),
    pval = c(2.90e-14, 0.06743, 0.01514))

set.seed(1)
LOND_results <- LOND(sample.df)

## -----------------------------------------------------------------------------
#after you've completed more experiments
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
LOND_results <- LOND(sample.df)

