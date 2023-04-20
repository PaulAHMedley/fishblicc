# Test table functions
library("tibble")


test_that("Prior table tibble exists", {
  expect_true(is_tibble(table_prior(eg_ld)))
})
