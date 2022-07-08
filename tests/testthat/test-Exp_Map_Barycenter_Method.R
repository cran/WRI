test_that("function Exp_Map_Barycenter_Method works for standard normal distribution", {
  
  # cdf of standard normal: by R code
  dSup <- seq(-3, 3, length.out=100)
  qSup <- seq(0, 1, length.out=length(dSup))
  qf <- qnorm(qSup) # Quantile function of standard normal N(0,1)
  
  cdf <- Exp_Map_Barycenter_Method(dSup, qf, qSup)

  # cdf of standard normal: by R code
  cdf_r <- pnorm(dSup)

  # Test
  expect_equal(cdf, cdf_r, tolerance=0.01)
})