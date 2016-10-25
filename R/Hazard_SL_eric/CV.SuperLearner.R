# V-fold Cross-Validation wrapper for Convex Super Learner function.  works with both family="gaussian" and "binomial"
# 
#  CV.SuperLearner.R
#  SuperLearner
#  
#  Created by Eric Polley on 2009-03-03.

CV.SuperLearner <- function(Y, X, SL.library, outside.V=20, inside.V=20, shuffle = TRUE, verbose = FALSE, family=gaussian(), method="NNLS", id=NULL, save.fit.library=FALSE, trim.logit=0.001, obsWeights = NULL, stratifyCV = FALSE, ...) {
	# create V-folds
	call <- match.call()
	n <- length(Y)
	folds <- CVFolds(V=outside.V, N.all = n, shuffle = shuffle, id = id, stratifyCV = stratifyCV, Y = Y)
	# if(shuffle) {
	# 	folds <- split(sample(1:n), rep(1:outside.V, length=n))
	# } else {
	# 	folds <- split(1:n, rep(1:outside.V, length=n))
	# }
	# if(!is.null(id)){
	# 	n.id <- length(unique(id))
	# 	id.split <- split(sample(1:n.id), rep(1:outside.V, length=n.id))
	# 	folds <- vector("list", outside.V)
	# 	for(v in seq(outside.V)) {
	# 		folds[[v]] <- which(id %in% unique(id)[id.split[[v]]])
	# 	}
	# }
	if(is.null(obsWeights)) {
		obsWeights <- rep(1, n)
	}
	if(!identical(length(obsWeights), n)) {
		stop("obsWeights vector must have the same dimension as Y")
	}
	# create placeholders for output
	library <- .createLibrary(SL.library)
	.check.SL.library(SL.library = c(library$predAlgorithm$library, library$screenAlgorithm))
	k.cand <- nrow(library$predAlgorithm)
	k.screen <- length(library$screenAlgorithm)
	cand.names <- rep(NA, k.cand)
	for(jj in seq(k.cand)) {
		cand.names[jj] <- paste(library$predAlgorithm[jj, 1], library$screenAlgorithm[library$predAlgorithm[jj, 2]], sep="_")
	}
	
	All.SL <- vector("list", outside.V)
	names(All.SL) <- paste("training", 1:outside.V, sep="")
	pred.SL <- rep(NA, n)
	pred.discreteSL <- rep(NA, n)
	whichDiscreteSL <- rep(NA, outside.V)
	pred.library <- matrix(NA, nrow=n, ncol=k.cand)
	colnames(pred.library) <- cand.names
	coef.SL <- matrix(NA, nrow=outside.V, ncol=k.cand)
	colnames(coef.SL) <- cand.names
	
	for(v in seq(outside.V)){
		valid <- folds[[v]]
		CV.Y <- Y[-valid]
		CV.X <- X[-valid, , drop = FALSE]
		CV.newX <- X[valid, , drop = FALSE]
		CV.id <- id[-valid]
		CV.obsWeights <- obsWeights[-valid]
		 
		fit.SL <- SuperLearner(Y=CV.Y, X=CV.X, newX=CV.newX, SL.library=SL.library, V=inside.V, shuffle=shuffle, verbose=verbose, family=family, method=method, id=CV.id, save.fit.library=save.fit.library, trim.logit=trim.logit, obsWeights = CV.obsWeights, stratifyCV = stratifyCV)
		if(verbose) {
			print(paste("completed CV ", v))
		}
		All.SL[[v]] <- fit.SL
		pred.SL[valid] <- fit.SL$SL.predict
		pred.discreteSL[valid] <- fit.SL$library.predict[, which.min(fit.SL$cv.risk)]
		whichDiscreteSL[v] <- names(which.min(fit.SL$cv.risk))
		pred.library[valid, ] <- fit.SL$library.predict
		coef.SL[v, ] <- fit.SL$coef
	}
	out <- list(CV.fit.SL=All.SL, pred.SL=pred.SL, pred.discreteSL = pred.discreteSL, whichDiscreteSL = whichDiscreteSL, pred.library=pred.library, coef.SL=coef.SL, folds=folds, method = method, call=call, obsWeights = obsWeights, id = id, V = c(outside.V = outside.V, inside.V = inside.V), Y = Y, X = X)
	class(out) <- c("CV.SuperLearner")
	return(out)
}