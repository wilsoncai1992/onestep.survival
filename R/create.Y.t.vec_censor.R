# compute I{T.tilde >= t, Delta = 1}, loop over t.vec
# compute I{T.tilde >= t}, loop over t.vec
# output:
# a binary vector, of length = t.vec

# 2016-10-23
create.Y.t.vec <- function(Time, t.vec) {
	out.vec <- (Time >= t.vec) + 0
	return(out.vec)
}

create.Y.t.vec_Delta1 <- function(Time, Delta, t.vec) {
    out.vec <- ((Time == t.vec) & (Delta == 1)) + 0
    return(out.vec)
}

