context("IT loglogistic")

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "loglogistic",
  model = "IT",
  N = NA,
  M = NA,
  SVR = 1,
  study = "Test loglogistic",
  Clevel = "arbitrary"
)

para <- c(hb = 0, kd = 1.3, t1 = 3, t2 = 2)

#print(guts)
test_that("new and old versions are similar (up to tolerance 1e-3) - used M = 500000, N = 100000 in old version", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = para), 
    c(1.0000000, 0.6861127, 0.5189911, 0.3004810, 0.2221791),
    tolerance = 1e-3
    )
  expect_equal(
    guts_calc_loglikelihood(guts, par = para), 
    -13.96819,
    tolerance = 1e-3
  )
  expect_equal(
    guts_report_sppe(guts), 
    -22.21791,
    tolerance = 1e-3
  )
  expect_equal(
    guts_report_squares(guts), 
    34.03945,
    tolerance = 1e-3
  )
})

guts_SVR_0.2 <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "loglogistic",
  model = "IT",
  N = NA,
  M = NA,
  SVR = 0.2,
  study = "Test loglogistic",
  Clevel = "arbitrary"
)

test_that("SVR correctly scales", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = c(para[1], 0.2 * para[2], para[3], para[4])), 
    guts_calc_survivalprobs(guts_SVR_0.2, par = para)
  )
})

para <- c(hb = 0.03, kd = 0.7, t1 = 2, t2 = 1.3)
test_that("new and old versions are similar (including hb) (up to tolerance 1e-3) - used M = 500000, N = 100000 in old version", {
  expect_equal(
    guts_calc_survivalprobs(guts, par = para), 
    c(1.0000000, 0.5847193, 0.4295508, 0.2825590, 0.2124839),
    tolerance = 1e-3
  )
  expect_equal(
    guts_calc_loglikelihood(guts, par = para), 
    -12.59041,
    tolerance = 1e-3
  )
  expect_equal(
    guts_report_sppe(guts), 
    -21.24839,
    tolerance = 1e-3
  )
  expect_equal(
    guts_report_squares(guts), 
    21.22359,
    tolerance = 1e-3
  )
})

guts1 <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "loglogistic",
  model = "IT",
  N = 10000,
  M = 10000,
  study = "Test loglogistic",
  Clevel = "arbitrary"
)

guts2 <- guts_setup(
  C = c(4, 5, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "loglogistic",
  model = "IT",
  N = 10000,
  M = 10000,
  study = "Test loglogistic",
  Clevel = "arbitrary"
)

guts_lst <- list(guts1, guts2)

para <- c(hb = 0.01, kd = 1.3, t1 = 3, t2 = 2)

test_that("multiple guts-objects in loops keep their values in calls by reference", {
  expect_equal(
    lapply(guts_lst, guts_calc_survivalprobs, par = para), 
    list(
      guts_calc_survivalprobs(guts1, para),
      guts_calc_survivalprobs(guts2, para)
    ),
    tolerance = 1e-6
  )
  expect_equal(
    sapply(guts_lst, guts_calc_loglikelihood, par = para), 
    c(
      guts_calc_loglikelihood(guts1, para),
      guts_calc_loglikelihood(guts2, para)
    )
  )
  expect_equal(
    guts_lst, list(guts1, guts2)
  )
})
