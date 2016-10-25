library(pracma)

l2.inner <- function(f, g, min, max, grid.eps = 1e-4) {
	# x.seq <- seq(0, 2*pi, by = grid.eps)
	# x.seq <- seq(min, max, by = grid.eps)
	x.seq <- seq(min, max, length.out = 1e4)

	# f <- function(x) rep(1, length(x))
	# g <- function(x) rep(1, length(x))

	func.val <- f(x.seq) * g(x.seq)
	# trapz(x.seq, func.val)/pi
	trapz(x.seq, func.val)
}


play1 <- function(x) rep(1, length(x))
play2 <- function(x) rep(1, length(x))
# l2.inner(play1, play2, min = 0, max = 2*pi, grid.eps = 1e-4)
l2.inner(play1, play2, min = 0, max = 2*pi)
# ================================================================================
# library(fda.usc)
# library(mgcv)