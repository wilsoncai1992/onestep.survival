
CV.SuperLearner.hazard <- function(N.delta, n.controls=1, time.df=5, time, X, SL.library, outside.V=20, inside.V = 20, shuffle=FALSE, verbose=FALSE, family=binomial(), method="NNLS", id=NULL, save.fit.library=FALSE, trim.logit=0.0001, discreteTime = TRUE, obsWeights = NULL) {
	require(nnls)
	require(gam)
	n <- length(N.delta)
	folds <- CVFolds(V=outside.V, N.all = n, shuffle = shuffle, id = id, stratifyCV = FALSE, Y = N.delta)
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
		CV.N.delta <- N.delta[-valid]
		CV.X <- X[-valid, ]
		CV.time <- time[-valid]
		CV.newX <- X[valid, ]
		CV.time.newX <- time[valid]
		CV.id <- id[-valid]
		CV.obsWeights <- obsWeights[-valid]
		 
		fit.SL <- SuperLearner.hazard(N.delta=CV.N.delta, n.controls = n.controls, time.df = time.df, time = CV.time, X=CV.X, time.newX = CV.time.newX, newX=CV.newX, SL.library=SL.library, V=inside.V, shuffle=shuffle, verbose=verbose, family=family, method=method, id=CV.id, save.fit.library=save.fit.library, trim.logit=trim.logit, obsWeights = CV.obsWeights)
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
	out <- list(CV.fit.SL=All.SL, pred.SL=pred.SL, pred.discreteSL = pred.discreteSL, whichDiscreteSL = whichDiscreteSL, pred.library=pred.library, coef.SL=coef.SL, folds=folds, method = method, call=call, obsWeights = obsWeights, id = id, V = c(outside.V = outside.V, inside.V = inside.V), N.delta = N.delta, time = time, X = X, n.controls = n.controls, time.df = time.df)
	class(out) <- c("CV.SuperLearner")
	return(out)
}
