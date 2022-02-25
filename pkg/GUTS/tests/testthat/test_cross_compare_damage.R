context("cross compare damage")

guts_IT <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "loglogistic",
  model = "IT",
  N = NA,
  M = NA,
  SVR = 1,
  study = "IT",
  Clevel = "arbitrary"
)

guts_SD <- guts_setup(
  C = c(4, 2, 4, 6, 6),
  Ct = seq_len(5) - 1,
  y = c(10,3,2,1,0),
  yt = seq_len(5) - 1,
  dist = "delta",
  model = "SD",
  M = 5000,
  N = NA,
  SVR = 1,
  study = "SD",
  Clevel = "arbitrary"
)

guts_calc_survivalprobs(guts_IT, par = c(hb = 0, kd = 1.3, t1 = 3, t2 = 2))
guts_calc_survivalprobs(guts_SD, par = c(hb = 0, kd = 1.3, kk = 0.1, t1 = 3))
damage_IT <- guts_report_damage(guts_IT)
damage_SD <- guts_report_damage(guts_SD)

#plot(damage_IT)
#points(damage_SD, col = "red")
find_close_timestep <- function(tim, tim_tab) {
  ind <- sapply(
    tim,
    function(x, x_tab) max(which(x_tab <= x)),
    x_tab = tim_tab
  )
}

test_that("damage estimates of SD and IT are similar (up to tolerance 1e-3)", {
  expect_equal(
    damage_IT$damage,
    damage_SD$damage[find_close_timestep(damage_IT$time, damage_SD$time)],
    tolerance = 1e-3
    )
})

#expos <- data.frame(
#	time = seq(0, 50) * 0.1,
#	conc = pmax(0,runif(51, min = -50, max = 1000))
#)
#dput(expos)

expos <- data.frame(
	time = c(
		0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
		0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 
		2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 
		3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 
		4.7, 4.8, 4.9, 5
	), 
	conc = c(
		555, 914, 785, 104, 20.8, 294, 912, 
		271, 140, 317, 764, 620, 92.3, 604, 889, 306, 231, 791, 943, 
		71.8, 399, 129, 41.3, 72.3, 904, 667, 0, 681, 702, 489, 825, 
		955, 214, 41.5, 647, 556, 241, 50.2, 95.4, 772, 627, 264, 977, 
		862, 0, 164, 197, 812, 297, 570, 158
	)
)

yt <- sort(unique(c(0, max(expos$time), runif(200, 0, max(expos$time)))))

guts_IT <- guts_setup(
	C = expos$conc,
	Ct = expos$time,
	y = c(100, rep(0, length(yt)-1)),
	yt = yt,
	dist = "loglogistic",
	model = "IT",
	N = NA,
	M = NA,
	SVR = 1,
	study = "IT",
	Clevel = "arbitrary"
)

guts_SD <- guts_setup(
	C = expos$conc,
	Ct = expos$time,
	y = c(100, rep(0, length(yt)-1)),
	yt = yt,
	dist = "delta",
	model = "SD",
	M = 100 * length(expos$time),
	N = NA,
	SVR = 1,
	study = "SD",
	Clevel = "arbitrary"
)

surv_IT <- guts_calc_survivalprobs(guts_IT, par = c(hb = 0, kd = 1.3, t1 = 3, t2 = 2))
surv_SD <- guts_calc_survivalprobs(guts_SD, par = c(hb = 0, kd = 1.3, kk = 0.1, t1 = 3))
damage_IT <- guts_report_damage(guts_IT)
damage_SD <- guts_report_damage(guts_SD)

#par(mfrow = c(2,1), mar = c(4,4,0.5,0.5))
#plot(guts_IT$Ct, guts_IT$C, type = "l")
#plot(damage_IT, type = "l")
#lines(damage_SD, col = "red")
#par(mfrow = c(1,1))
test_that("damage estimates of SD and IT for complex exposure are similar (up to tolerance 1e-3)", {
	expect_equal(
		damage_IT$damage,
		damage_SD$damage[find_close_timestep(damage_IT$time, damage_SD$time)],
		tolerance = 1e-3
	)
})
