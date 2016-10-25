# compute l2 inner product of two step functions: f and g
# input: two step-function pdf with shared jump points (T.grid)
# can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
# output: scalar

l2.inner.step <- function(f.step, g.step, T.grid) {
	source('./compute.step.cdf_vec.R')	
	if(is.vector(f.step) & is.vector(g.step)){
		# both f and g are one sample
		f.times.g <- f.step * g.step
	}
	if(!is.vector(f.step) & is.vector(g.step)){
		# f: multi-sample
		# g: one-sample
		f.times.g <- sweep(f.step,MARGIN=2,g.step,`*`) # multiply g to each row of f.
	}
	if(is.vector(f.step) & !is.vector(g.step)){
		# f: one-sample
		# g: multi-sample
		f.times.g <- sweep(g.step,MARGIN=2,f.step,`*`) # multiply f to each row of g.
	}
	if(!is.vector(f.step) & !is.vector(g.step)){
		# both f and g are multi-sample of same sample size
		if(nrow(f.step) != nrow(g.step)) stop('f and g have different sample size!')
		f.times.g <- f.step * g.step
	}
	# ------------------------------------------------------------------------------------
	result <- compute.step.cdf(f.times.g, T.grid)
	if(!is.vector(f.step) | !is.vector(g.step)){
		# there is multi-sample
		result <- apply(result, 1, function(obj) tail(obj, 1))
	}else{
		# both f and g are one-sample
		result <- tail(result, 1)
	}
	return(result)
}



# haha <- l2.inner(sin, cos, -pi, pi)
# haha$value
# haha <- l2.inner(Pn.D1.t.func, Pn.D1.t.func, min = 0, max = T.max)

# play <- l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq)
# 
# library(rbenchmark)
# benchmark(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq), 
# 					l2.inner(Pn.D1.t.func, Pn.D1.t.func, min = 0, max = T.max), 
# 					replications = 1e4)
# ----------------------------------------------------------------------------------------------------
# play <- l2.inner(rep(c(1,0), length(T.uniq)/2), rep(1, length(T.uniq)), T.uniq)
# play <- l2.inner(rep(c(1,0), length(T.uniq)/2), rep(c(0,1), length(T.uniq)/2), T.uniq)
# 0

# haha <- l2.inner(gen.step.func(c(0,rep(c(500,0), length(T.uniq)/2)), T.uniq), 
								 # gen.step.func(c(0,rep(c(0,0.8), length(T.uniq)/2)), T.uniq), min = 0, max = T.max)
