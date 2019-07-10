test.pval <- c(1e-07, 3e-04, 0.1, 6e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)

test_that("Correct rejections", {
    expect_identical(ADDIS(test.pval)$R, c(1,1,0,1))
    
    expect_identical(ADDIS(test.pval,
                           async=TRUE, decision.times=seq_len(4)+3)$R,
                     c(1,1,0,0))
})

test_that("Check that ADDIS with async=FALSE is a special case of async=TRUE", {
              expect_equal(ADDIS(test.pval)$alphai,
                           ADDIS(test.pval, async=TRUE,
                                 decision.times=seq_len(4))$alphai)
})

test_that("Check that ADDIS gives same results as SAFFRON and SAFFRONstar
          algorithms with discard=TRUE", {
              expect_equal(ADDIS(test.pval)$alphai,
                           SAFFRON(test.df, discard=TRUE)$alphai)
              
              expect_equal(ADDIS(test.pval, async=TRUE,
                                 decision.times=seq_len(4)+2)$alphai,
                           SAFFRONstar(test.pval, version='async', discard=TRUE,
                                       decision.times=seq_len(4)+2)$alphai)
})
