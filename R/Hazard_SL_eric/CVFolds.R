# CVFolds creates the list of row number in each of the V folds.  special cases are when shuffling row number, id cluster identification and if stratify by the outcome to maintain (near) balance in each fold.

CVFolds <- function(V, N.all, shuffle, id, stratifyCV, Y) {
	if(!stratifyCV) {
		if(shuffle) {
			if(is.null(id)) {
				DATA.split <- split(sample(1:N.all), rep(1:V, length=N.all))
			} else {
				n.id <- length(unique(id))
				id.split <- split(sample(1:n.id), rep(1:V, length=n.id))
				DATA.split <- vector("list", V)
				for(v in seq(V)) {
					DATA.split[[v]] <- which(id %in% unique(id)[id.split[[v]]])
				}
			}
			
		} else {
			if(is.null(id)) {
				DATA.split <- split(1:N.all, rep(1:V, length=N.all))
			} else {
				n.id <- length(unique(id))
				id.split <- split(1:n.id, rep(1:V, length=n.id))
				DATA.split <- vector("list", V)
				for(v in seq(V)) {
					DATA.split[[v]] <- which(id %in% unique(id)[id.split[[v]]])
				}
			}
		}
	} else {
		if(length(unique(Y)) != 2) {
			stop("stratifyCV only implemented for binary Y")
		}
		if(sum(Y) < V | sum(!Y) < V) {
			stop("number of (Y=1) or (Y=0) is less than the number of folds")
		}
		if(shuffle) {
			if(is.null(id)) {
				within.split <- suppressWarnings(tapply(1:N.all, INDEX = Y, FUN = split, rep(1:V)))
				DATA.split <- vector("list", length = V)
				names(DATA.split) <- paste(seq(V))
				for(vv in seq(V)) {
					DATA.split[[vv]] <- c(within.split[[1]][[vv]], within.split[[2]][[vv]])
				}
			} else {
				stop("stratified sampling with id not currently implemented")
			}
		} else {
			if(is.null(id)) {
				within.split <- suppressWarnings(tapply(1:N.all, INDEX = Y, FUN = split, rep(1:V)))
				DATA.split <- vector("list", length = V)
				names(DATA.split) <- paste(seq(V))
				for(vv in seq(V)) {
					DATA.split[[vv]] <- c(within.split[[1]][[vv]], within.split[[2]][[vv]])
				}
			} else {
				stop("stratified sampling with id not currently implemented")
			}
		}
	}	
	invisible(DATA.split)
}


# # testing
# N.all <- 200
# Y <- rbinom(N.all, 1, 0.2)
# V <- 10
# shuffle <- TRUE
# id <- NULL
# stratifyCV <- TRUE
# 
# within.split <- suppressWarnings(tapply(1:N.all, INDEX = Y, FUN = split, rep(1:V)))
# DATA.split <- vector("list", length = V)
# 
# for(vv in seq(V)) {
# 	DATA.split[[vv]] <- c(within.split[[1]][[vv]], within.split[[2]][[vv]])
# }
# 
# # check
# for(i in seq(V)) {
# 	print(mean(Y[DATA.split[[i]]]))
# }
# 
# fooS <- CVFolds(V=V, N.all = N.all, shuffle = FALSE, id = NULL, stratifyCV = TRUE, Y = Y)
# fooN <- CVFolds(V=V, N.all = N.all, shuffle = FALSE, id = NULL, stratifyCV = FALSE, Y = Y)
# 
# for(i in seq(V)) {
# 	print(mean(Y[fooS[[i]]]))
# }
# for(i in seq(V)) {
# 	print(mean(Y[fooN[[i]]]))
# }