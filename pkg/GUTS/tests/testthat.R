path <- file.path("test", "testthat")
test_files <- dir(path, pattern = "test_")

sapply(file.path(path, test_files), source)
