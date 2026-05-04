test_that("blicc_expected_catches produces tibble ",{
  expect_true(tibble::is_tibble(blicc_expected_catches(trgl_rp)))
})


test_that("blicc reference points produces rp_df tibble ",{
  gl_fit <- blicc_mpd(gillnet_ld)
  gl_rp <- blicc_ref_pts(gl_fit, gillnet_ld)
  expect_true(tibble::is_tibble(gl_rp$rp_df) & tibble::is_tibble(gl_rp$dr_df))
})
  


test_that("blicc reference points produces rp_df tibble ",{
  gl_rp <- blicc_ref_pts(trgl_slim, trgl_ld, vdir = c(1,0))
  expect_true(tibble::is_tibble(gl_rp$rp_df) & tibble::is_tibble(gl_rp$dr_df) & is.list(gl_rp$scenario))
})
