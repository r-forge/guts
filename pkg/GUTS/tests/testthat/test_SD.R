context("SD")

guts <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "",
  model = "SD",
  N = 10000,
  M = 10000,
  study = "Test loglogistic",
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
