library(PharmacoGx)

context("Testing .Hill function")

test_that("Returns right maths",{

    expect_equal(.Hill(0, c(1, 0, -Inf)), 0)
    expect_equal(.Hill(0, c(0, 0, 0)), 1/2)
    expect_equal(.Hill(0, c(1, 0, Inf)), 1)
    expect_equal(.Hill(-Inf, c(1, 0.2, 1)), 1)
    expect_equal(.Hill(Inf, c(1, 0.2, 1)), 0.2)
})
