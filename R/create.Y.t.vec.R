# compute I{T > t}, loop over t.vec
# output:
# a binary vector, of length = t.vec

# 2016-05-24
create.Y.t.vec <- function(Time, t.vec) {
	out.vec <- (Time > t.vec) + 0
	return(out.vec)
}