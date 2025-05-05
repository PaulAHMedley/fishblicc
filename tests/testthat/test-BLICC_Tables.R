# Test table functions


test_that("Prior table tibble exists", {
  expect_true(tibble::is_tibble(blicc_prior(gillnet_ld)))
})


test_that("Prior table tibble exists", {
  expect_true(tibble::is_tibble(blicc_prior(blip_Mk(trgl_ld, ref_length=25))))
})


test_that("Results table tibble exists", {
  expect_true(tibble::is_tibble(blicc_results(trgl_rp)))
})
