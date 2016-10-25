compute.Qn <- function(qn.current) {
	# cumsum(rev(qn.current[1,]))
	qn.curr.func <- gen.step.func(c(qn.current[1,], 0), T.uniq)
	for (it1 in 1:length(T.uniq)) {
		T.it <- T.uniq[it1]
		integrate(qn.curr.func, lower = T.it, upper = T.max, subdivisions = 1e5)$value
	}
}