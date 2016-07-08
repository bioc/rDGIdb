library(rDGIdb)

test_that("Wrong input arguments", {
    expect_error(plotInteractionsBySource(NULL))
    expect_error(plotInteractionsBySource("a"))
    expect_error(plotInteractionsBySource(queryDGIdb('XYZA')))
})
