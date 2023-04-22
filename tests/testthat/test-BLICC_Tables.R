# Test table functions
library("tibble")


test_that("Prior table tibble exists", {
  expect_true(is_tibble(table_prior(eg_ld)))
})


test_that("Prior table tibble exists", {
  expect_true(is_tibble(table_prior(set_Mk(eg_ld, ref_length=25))))
})
