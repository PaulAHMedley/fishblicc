
test_that("Mpd produces tibble for 1 gear",{
  ld1 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                   sel_fun=4L,
                   gear_names = "Otter Trawl",
                   Linf=c(35, 3),
                   L50=20,
                   L95=21,
                   Mk = 1.5)
  expect_true(tibble::is_tibble(blicc_mpd(ld1)))
})


test_that("Fit produces stanfit object.",{
  ld1 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                   sel_fun=4L,
                   gear_names = "Otter Trawl",
                   Linf=c(35, 3),
                   L50=20,
                   L95=21,
                   Mk = 1.5)
  suppressWarnings(
  res <- blicc_fit(ld1, ntarget=100, nwarmup = 100, nchain=1)
  )
  expect_true(class(res)[1]=="stanfit")
})


test_that("Mpd produces tibble for 1 gear",{
  ld1 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                   sel_fun=4L,
                   gear_names = "Otter Trawl",
                   Linf=c(35, 3),
                   L50=20,
                   L95=21,
                   Mk = 1.5)
  expect_true(tibble::is_tibble(blicc_mpd(ld1)))
})

test_that("Mpd produces tibble for 1 gear",{
  ld1 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                   sel_fun=4L,
                   gear_names = "Otter Trawl",
                   Linf=c(35, 3),
                   L50=20,
                   L95=21,
                   Mk = 1.5)
  expect_true(tibble::is_tibble(blicc_mpd(ld1)))
})


test_that("Mpd produces tibble for 2 gears ",{
  ld2 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=list(c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                           c(0,1,1,5,6,12,24,34,31,30,20,15,10,5,3,0)),
                   Catch = c(0.6, 0.4),
                   sel_fun = c(4L, 4L),
                   gear_names = c("Otter Trawl", "Handline"),
                   Linf=c(35, 3),
                   Mk = 1.5)
  expect_true(tibble::is_tibble(blicc_mpd(ld2)))
})

test_that("Change selectivity functions ",{
  ld2 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=list(c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                           c(0,1,1,5,6,12,24,34,31,30,20,15,10,5,3,0)),
                   Catch = c(0.6, 0.4),
                   sel_fun = c(4L, 4L),
                   gear_names = c("Otter Trawl", "Handline"),
                   Linf=c(35, 3),
                   Mk = 1.5)
  ld2 <- blicc_selfun(ld2, sel_fun = c("logistic", "normal"),
                      sel_indx = c(1, 2),
                      model_name = "Test")
  expect_true(is.list(ld2))
})


test_that("Change selectivity mixtures ",{
  ld2 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=list(c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                           c(0,1,1,5,6,12,24,34,31,30,20,15,10,5,3,0)),
                   Catch = c(0.6, 0.4),
                   sel_fun = c("logistic", "normal"),
                   gear_names = c("Otter Trawl", "Handline"),
                   Linf=c(35, 3),
                   Mk = 1.5)
  ld2 <- blicc_gear_sel(ld2, gear_sel = list("Otter Trawl" = c(1, 2),
                                             "Handline" = c(2)),
                        gear = c("Otter Trawl", "Handline"),
                        model_name = "Test")
  expect_true(is.list(ld2))
})

test_that("Remove population ",{
  ld2 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=list(c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                           c(0,1,1,5,6,12,24,34,31,30,20,15,10,5,3,0)),
                   period_freq = c(1, 2),
                   period_names = c("1995", "2016"),
                   sel_fun = "logistic",
                   gear_names = "Otter Trawl",
                   Linf=c(35, 3),
                   Mk = 1.5)
  ld2 <- blicc_population_filter(ld2, population=2)
  expect_true(is.list(ld2) & ld2$NN==1L)
})

  

