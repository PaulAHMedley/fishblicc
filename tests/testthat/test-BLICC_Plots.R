library("ggplot2")
library("testthat")

# Check works with MPD results
ld <- blicc_dat(LLB = 25:40, fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                 Linf=c(35, 3))
mpd_fit <- blicc_mpd(ld)
mpdrp_df <- blicc_ref_pts(mpd_fit, ld)
mpdlx_df <- blicc_expect_len(mpdrp_df, ld)


test_that("Plot expected frequency exist",{
  p <- plot_expected_frequency(mpdrp_df, mpdlx_df, ld)
  expect_true(is.ggplot(p))
})

test_that("Plot selectivity exists",{
  p <- plot_selectivity(mpdlx_df, ld)
  expect_true(is.ggplot(p))
})

test_that("Plot SPR density rejected",{
  p <- plot_SPR_density(mpdrp_df)
  expect_identical(p, "To obtain a density, you will need to obtain sufficient values (>100) from MCMC.")
})

test_that("Plot FkF40 density rejected",{
  p <- plot_FkF40_density(mpdrp_df)
  expect_identical(p, "To obtain a density, you will need to obtain sufficient values (>100) from MCMC.")
})

test_that("Plot expected frequencies for F reference points exists",{
  p <- plot_efq_FRP(mpdrp_df, mpdlx_df, ld)
  expect_true(is.ggplot(p))
})

test_that("Plot expected frequencies for selectivity reference points exists",{
  p <- plot_efq_SRP(mpdrp_df, mpdlx_df, ld)
  expect_true(is.ggplot(p))
})

test_that("Plot SPR contour plot exists",{
  p <- plot_SPR_contour(mpdrp_df, ld)
  expect_true(is.ggplot(p))
})

test_that("Plot YPR contour plot exists",{
  p <- plot_YPR_contour(mpdrp_df, ld)
  expect_true(is.ggplot(p))
})



