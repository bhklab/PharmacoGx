test_that("signature thread cap falls back to one when CPU count is unknown", {
  expect_equal(.capSignatureThreads(1L, available = NA_integer_), 1L)
  expect_equal(.capSignatureThreads(4L, available = NA_integer_), 1L)
  expect_equal(.capSignatureThreads(4L, available = 2L), 2L)
  expect_equal(.capSignatureThreads(2L, available = 4L), 2L)
})
