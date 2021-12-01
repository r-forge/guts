context("IT Proper")

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "loglogistic",
  model = "Proper",
  N = 10000,
  M = 10000,
  study = "Test proper loglogistic",
  Clevel = "arbitrary"
)

para <- c(hb = 0, kd = 1.3, kk = 0.07, alpha = 3, beta = 2)

#print(guts)
test_that("new and old versions are similar (up to tolerance 1e-5)", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = para), 
    c(1.0000000, 0.9910818, 0.9678977, 0.9052989, 0.7982619),
    tolerance = 1e-5
    )
  expect_equal(
    guts_calc_loglikelihood(guts, par = para), 
    -41.80747,
    tolerance = 1e-5
  )
  expect_equal(
    guts_report_sppe(guts), 
    -79.82619,
    tolerance = 1e-5
  )
  expect_equal(
    guts_report_squares(guts), 
    235.2989,
    tolerance = 1e-5
  )
})

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "lognormal",
  model = "Proper",
  N = 10000,
  M = 10000,
  study = "Test lognormal",
  Clevel = "arbitrary"
)

para <- c(hb = 0, kd = 1.3, kk = 0.07, mn = 3, sd = 2)
#print(guts)
test_that("new and old versions are similar (up to tolerance 1e-5)", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = para), 
    c(1.0000000, 0.9923859, 0.9683345, 0.8941076, 0.7645970),
    tolerance = 1e-5
  )
  expect_equal(
    guts_calc_loglikelihood(guts, par = para), 
    -42.51648,
    tolerance = 1e-5
  )
  expect_equal(
    guts_report_sppe(guts), 
    -76.4597,
    tolerance = 1e-5
  )
  expect_equal(
    guts_report_squares(guts), 
    228.4952,
    tolerance = 1e-5
  )
})

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "delta",
  model = "Proper",
  N = 10000,
  M = 10000,
  study = "Test proper delta",
  Clevel = "arbitrary"
)

para <- c(hb = 1e-5, kd = 1.3, kk = 0.1, t1 = 3)

#print(guts)
test_that("new and old versions are similar (up to tolerance 1e-5)", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = para), 
    c(1.0000000, 0.9999900, 0.9999800, 0.9319453, 0.7475945),
    tolerance = 1e-5
  )
  expect_equal(
    guts_calc_loglikelihood(guts, par = para), 
    -96.48211,
    tolerance = 1e-5
  )
  expect_equal(
    guts_report_sppe(guts), 
    -74.75945,
    tolerance = 1e-5
  )
  expect_equal(
    guts_report_squares(guts), 
    238.0985,
    tolerance = 1e-5
  )
})

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "external",
  model = "Proper",
  N = 10000,
  M = 10000,
  study = "Test proper external",
  Clevel = "arbitrary"
)
para <- c(0.051, 0.126, 1.618, 19.099, 6.495)
sigma2 <- log( 1 + 6.495^2 / 19.099^2)
mu <- log(19.099) - 0.5 * sigma2
lognormal.thresholds <- sort(rlnorm(100, meanlog = mu, sdlog = sqrt(sigma2)))
lognormal.thresholds.orig <- lognormal.thresholds

gts.lognormal <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1, 
  dist = "lognormal", model = "Proper")

test_that("lognormal external distribution and internal distribution are similar (up to tolerance 1e-2 -- might fail due to random evaluation)", {
  expect_equal(
    guts_calc_loglikelihood(guts, par = para[1:3], external_dist = lognormal.thresholds), 
    guts_calc_loglikelihood(gts.lognormal, c(0.051, 0.126, 1.618, 19.099, 6.495)), 
    tolerance = 1e-2
  )
  expect_equal(
    guts_calc_survivalprobs(guts, par = para[1:3], external_dist = lognormal.thresholds), 
    guts_calc_survivalprobs(gts.lognormal, c(0.051, 0.126, 1.618, 19.099, 6.495)), 
    tolerance = 1e-2
  )
})

test_that("lognormal external distribution has not changed during manipulation (Rcpp correctly clones the distribution vector)", {
  expect_equal(
    lognormal.thresholds, 
    lognormal.thresholds.orig
  )
})

data(diazinon)
gts.lognormal <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  M = 10000,
  dist = "lognormal", model = "Proper")

sigma2 <- log( 1 + 6.495^2 / 19.099^2)
mu <- log(19.099) - 0.5 * sigma2
lognormal.thresholds <- sort(rlnorm(10000, meanlog = mu, sdlog = sqrt(sigma2)))
gts.external <- guts_setup(
  C = diazinon$C1, Ct = diazinon$Ct1,
  y = diazinon$y1, yt = diazinon$yt1, 
  dist = "external", model = "Proper",
  N = 1000,
  M = 10000,
  study = "Test proper external",
  Clevel = "arbitrary")

test_that("lognormal external distribution and internal distribution are similar for diazinon (up to tolerance 1e-2 -- might fail due to random evaluation)", {
  expect_equal(
    guts_calc_loglikelihood(gts.external, c(0.051, 0.126, 1.618), lognormal.thresholds), 
    guts_calc_loglikelihood(gts.lognormal, c(0.051, 0.126, 1.618, 19.099, 6.495)), 
    tolerance = 1e-2
  )
  expect_equal(
    guts_calc_survivalprobs(gts.external, c(0.051, 0.126, 1.618), lognormal.thresholds), 
    guts_calc_survivalprobs(gts.lognormal, c(0.051, 0.126, 1.618, 19.099, 6.495)), 
    tolerance = 1e-2
  )
})
