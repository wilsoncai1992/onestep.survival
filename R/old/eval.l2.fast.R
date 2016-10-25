l2.inner <- function(f, g, min, max) {
	# f.times.g <- function(f, g) {
		out.func <- function(x) {f(x) * g(x)}
		# return(out.func)
	# }
	
	integrate(out.func, min, max, subdivisions = 1e5L)$value
}



# haha <- l2.inner(sin, cos, -pi, pi)
# haha$value
