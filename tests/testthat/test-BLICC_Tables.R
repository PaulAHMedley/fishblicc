# Test table functions


test_that("Prior table tibble exists", {
  expect_true(tibble::is_tibble(blicc_prior(eg_ld)))
})


test_that("Prior table tibble exists", {
  expect_true(tibble::is_tibble(blicc_prior(blip_Mk(eg_ld, ref_length=25))))
})


test_that("Results table tibble exists", {
  expect_true(tibble::is_tibble(blicc_results(eg_rp)))
})
