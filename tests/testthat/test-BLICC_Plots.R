library("ggplot2")
library("testthat")

# The following carries out basic checks that plot routines run and produce
# the expected output. If they do not, it might imply something wrong
# in other routines.
# Subtleties in output are not tested. This requires analytical evaluation
# rather software tests.


# Check works with MPD results
suppressWarnings(
  ld <- blicc_dat(model_name = "Test Graphs",
                  LLB = 25:40,
                  fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                  sel_fun=4L,
                  gear_names = "Otter Trawl",
                  Linf=c(35, 3))
)
mpd_fit <- blicc_mpd(ld)
mpd_rp <- blicc_ref_pts(mpd_fit, ld)

# plot_expected_frequency

test_that("Plot mpd expected frequency exist",{
  p <- plot_expected_frequency(mpd_rp)
  expect_true(ggplot2::is.ggplot(p))
})

test_that("Plot mcmc expected frequency exist",{
  p <- plot_expected_frequency(eg_rp)
  expect_true(ggplot2::is.ggplot(p))
})

# plot_residuals

test_that("Plot mpd residuals exist",{
  p <- plot_residuals(mpd_rp)
  expect_true(ggplot2::is.ggplot(p))
})

test_that("Plot mcmc residuals exist",{
  p <- plot_residuals(eg_rp)
  expect_true(ggplot2::is.ggplot(p))
})

# plot_selectivity

test_that("Plot mpd selectivity exists",{
  p <- plot_selectivity(mpd_rp)
  expect_true(ggplot2::is.ggplot(p))
})


test_that("Plot mcmc selectivity exists",{
  p <- plot_selectivity(eg_rp)
  expect_true(ggplot2::is.ggplot(p))
})

# plot_SPR_density
test_that("Plot mpd SPR density rejected",{
  expect_error(plot_SPR_density(mpd_rp))
})

test_that("Plot mcmc SPR density exists",{
  p <- plot_SPR_density(eg_rp)
  expect_true(ggplot2::is.ggplot(p))
})

# plot_FkF40_density

test_that("Plot mpd FkF40 density rejected",{
  expect_error(plot_FkF40_density(mpd_rp))
})

test_that("Plot mcmc FkF40 density exists",{
  p <- plot_FkF40_density(eg_rp)
  expect_true(ggplot2::is.ggplot(p))
})

# plot_efq_FRP

test_that("Plot mpd expected frequencies for F reference points exists",{
  p <- plot_efq_FRP(mpd_rp)
  expect_true(is.ggplot(p))
})

test_that("Plot mcmc expected frequencies for F reference points exists",{
  p <- plot_efq_FRP(eg_rp)
  expect_true(is.ggplot(p))
})

# plot_efq_SRP

test_that("Plot mpd expected frequencies for selectivity reference points exists",{
  p <- plot_efq_SRP(mpd_rp)
  expect_true(is.ggplot(p))
})

test_that("Plot mcmc expected frequencies for selectivity reference points exists",{
  p <- plot_efq_SRP(eg_rp)
  expect_true(is.ggplot(p))
})

# plot_SPR_contour

test_that("Plot mpd SPR contour plot exists",{
  p <- plot_SPR_contour(mpd_rp)
  expect_true(is.ggplot(p))
})

test_that("Plot mcmc SPR contour plot exists",{
  p <- plot_SPR_contour(eg_rp)
  expect_true(is.ggplot(p))
})

# plot_YPR_contour

test_that("Plot mpd YPR contour plot exists",{
  p <- plot_YPR_contour(mpd_rp)
  expect_true(is.ggplot(p))
})


test_that("Plot mcmc YPR contour plot exists",{
  p <- plot_YPR_contour(eg_rp)
  expect_true(is.ggplot(p))
})

