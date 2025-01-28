context("iNEXT")
test_that("iNEXT for abundance-based data", {
  # Test input by a demo data
  data(spider)
  out <- iNEXT(spider, q=0, datatype="abundance")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "n")
  expect_equal(nrow(out$DataInfo), length(spider))
  
  # Test input by a vector
  x <- spider$Girdled
  out <- iNEXT(x, q=0, datatype="abundance")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "n")
  expect_equal(nrow(out$DataInfo), 1)
  
  # Test input by a data.frame
  x <- data.frame(a=c(10,20,30,40,50,0,0), b=c(11,22,0,0,33,44,55))
  out <- iNEXT(x, q=0, datatype="abundance")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "n")
  expect_equal(nrow(out$DataInfo), ncol(x))
  
})

test_that("iNEXT for sampling-unit-based incidence frequencies data", {
  # Test input by a demo data
  data(ant)
  out <- iNEXT(ant, q=0, datatype="incidence_freq")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), length(ant))
  
  # Test input by a vector
  out <- iNEXT(ant$h50m, q=0, datatype="incidence_freq")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), 1)
})


test_that("iNEXT for species by sampling-units incidence matrix", {
  # Test input by a demo data
  data(ciliates)
  options(warn=-1)
  out <- iNEXT(ciliates, q=0, datatype="incidence_raw")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), length(ciliates))
  
  # Test input by a data.frame
  x <- ciliates$EtoshaPan
  expect_equal(class(x)[1], "matrix")
  out <- iNEXT(x, q=0, datatype="incidence_raw")
  expect_is(out, "iNEXT")
  expect_output(str(out), "List of 3")
  expect_equal(names(out$DataInfo)[2], "T")
  expect_equal(nrow(out$DataInfo), 1)
})


