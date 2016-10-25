y0 <- c(1., 2., 4., 3.)
sfun0  <- stepfun(1:3, y0, f = 0)

curve(sfun0, from = -1, to = 4)
integrate(sfun0, lower = 1, upper = 3)
