library("tibble")
suppressWarnings(
  ld1 <- blicc_dat(model_name = "Test Graphs",
                  LLB = 25:40,
                  fq=c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                  sel_fun=4L,
                  gear_names = "Otter Trawl",
                  Linf=c(35, 3))
)
suppressWarnings(
  ld2 <- blicc_dat(model_name = "Test Graphs",
                   LLB = 25:40,
                   fq=list(c(0,1,2,26,72,66,40,36,31,25,19,12,10,8,2,0),
                           c(0,1,1,5,6,12,24,34,31,30,20,15,10,5,3,0)),
                   Catch = c(0.6, 0.4),
                   sel_fun = c(4L, 1L),
                   gear_names = c("Otter Trawl", "Handline"),
                   Linf=c(35, 3))
)


test_that("Mpd produces tibble for 1 gear",{
  expect_true(is_tibble(blicc_mpd(ld1)))
})
test_that("Mpd produces tibble for 2 gears ",{
  expect_true(is_tibble(blicc_mpd(ld2)))
})
