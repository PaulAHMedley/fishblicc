test_that("blicc_expected_catches produces tibble ",{
  expect_true(tibble::is_tibble(blicc_expected_catches(trgl_rp)))
})

gl_fit <- blicc_mpd(gillnet_ld)
gl_rp <- blicc_ref_pts(gl_fit, gillnet_ld)

test_that("blicc reference points produces rp_df tibble ",{
  expect_true(tibble::is_tibble(gl_rp$rp_df))
})
  