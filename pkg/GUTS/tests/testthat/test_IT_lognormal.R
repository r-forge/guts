context("IT lognormal")

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "lognormal",
  model = "IT"
)

para <- c(hb = 0, kd = 1.3, t1 = 3, t2 = 2)

#print(guts)
test_that("new and old versions are similar (up to tolerance 1e-3)", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = para), 
    c(1.00000000, 0.63369535, 0.40483971, 0.15860946, 0.09072692),
    tolerance = 1e-3
    )
  expect_equal(
    guts_calc_loglikelihood(guts, par = para), 
    -12.59616,
    tolerance = 1e-3
  )
  expect_equal(
    guts_report_sppe(guts), 
    -9.072692,
    tolerance = 1e-3
  )
  expect_equal(
    guts_report_squares(guts), 
    16.49783,
    tolerance = 1e-3
  )
})


data(diazinon)
gts.lognormal <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "lognormal", model = "IT")

sigma2 <- log( 1 + 6.495^2 / 19.099^2)
mu <- log(19.099) - 0.5 * sigma2
lognormal.thresholds <- sort(rlnorm(10000, meanlog = mu, sdlog = sqrt(sigma2)))
gts.external <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "external", model = "IT",
  study = "Test IT external",
  Clevel = "arbitrary")

test_that("lognormal external distribution and internal distribution are similar for diazinon (up to tolerance 1e-2 -- might fail due to random evaluation)", {
  expect_equal(
    guts_calc_loglikelihood(gts.external, c(0.051, 0.126), lognormal.thresholds), 
    guts_calc_loglikelihood(gts.lognormal, c(0.051, 0.126, 19.099, 6.495)), 
    tolerance = 1e-2
  )
  expect_equal(
    guts_calc_survivalprobs(gts.external, c(0.051, 0.126), lognormal.thresholds), 
    guts_calc_survivalprobs(gts.lognormal, c(0.051, 0.126, 19.099, 6.495)), 
    tolerance = 1e-2
  )
})

test_that("Wrong parameter settings throw error", {
  expect_error(
    guts_calc_loglikelihood(gts.external, c(0.051, 0.126, 4, 7), lognormal.thresholds)
  )
  expect_error(
    guts_calc_loglikelihood(gts.external, c(0.051, 0.126))
  )
  expect_error(
    guts_calc_survivalprobs(gts.external, c(0.051, 0.126, 4, 7), lognormal.thresholds)
  )
  expect_error(
    guts_calc_survivalprobs(gts.lognormal, c(0.051, 0.126, 19.099, 6.495, 17))
  )
  expect_error(
    guts_calc_loglikelihood(gts.lognormal, 1)
  )
})  
