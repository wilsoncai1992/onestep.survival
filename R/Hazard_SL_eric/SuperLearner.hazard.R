# Convex Super Learner function for Hazard estimation.
# 
#  SuperLearner.R
#  SuperLearner
#  
#  Created by Eric Polley on 2009-07-03.
#  Copyright 2009 
# 

# basically identical to binomial SuperLearner. Additional arguements:
# n.controls -- number of "controls" in incidence sampling
# time -- possibly delta.u from create.discrete()
# expects the input to be in the form from create.discrete()
# Y == N.delta
# X should not include time

SuperLearner.hazard <- function(N.delta, n.controls=1, time.df=5, time, X, newX, time.newX, SL.library, V=20, shuffle=FALSE, verbose=FALSE, family=binomial(), method="NNLS", id=NULL, save.fit.library=TRUE, trim.logit=0.0001, discreteTime = TRUE, obsWeights = NULL) {
	require(nnls)
	require(gam)
	library <- .createLibrary(SL.library)
	.check.SL.library(SL.library = c(library$predAlgorithm$library, library$screenAlgorithm))
	call <- match.call()
	X <- as.data.frame(X)
	var.names <- names(X)
	N.all <- length(N.delta)
	p <- ncol(X)
	k.cand <- nrow(library$predAlgorithm)
	k.screen <- length(library$screenAlgorithm)
	newZ.X <- newZ <- matrix(NA, nrow=N.all, ncol=k.cand)
	
	cand.names <- rep(NA, k.cand)
	for(jj in seq(k.cand)){
		cand.names[jj] <- paste(library$predAlgorithm[jj, 1], library$screenAlgorithm[library$predAlgorithm[jj, 2]], sep="_")
	}
	errorsInCVLibrary <- rep(0, k.cand)
	errorsInLibrary <- rep(0, k.cand)
	
	if(missing(newX)) {
		newX <- X
		time.newX <- time
	}
	# if(missing(time.newX)){time.newX <- time}
	
	fit.gam.library <- fit.library <- vector("list", k.cand)
	names(fit.library) <- cand.names
	
	if(is.character(family))
		family <- get(family, mode="function", envir=parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
	DATA.split <- CVFolds(V=V, N.all = N.all, shuffle = shuffle, id = id, stratifyCV = FALSE, Y = N.delta)
	if(is.null(id)) {
		id <- seq(N.all)
	}
	if(!identical(length(id), N.all)) {
		stop("id vector must have the same dimension as Y")
	}
	if(is.null(obsWeights)) {
		obsWeights <- rep(1, N.all)
	}
	if(!identical(length(obsWeights), N.all)) {
		stop("obsWeights vector must have the same dimension as Y")
	}
	X <- as.data.frame(X)
	newX <- as.data.frame(newX)
	namesX <- names(X)
	if(!identical(colnames(X), colnames(newX))) {
		stop("The variable names and order in newX must be identical to the variable names and order in X")
	}
	# test for character of factor variables
	if(sum(!sapply(X, is.numeric)) > 0) {
		stop("Currently, only numeric variables are allowed.  Please convert any character or factor variables to numeric.")
	}
	# test for missing values
	if(sum(is.na(N.delta)) > 0) {
		stop("missing data is currently not supported")
	}
	if(!is.numeric(N.delta)) {
		stop("the outcome Y must be a numeric vector")
	}
	
	## now for the candidates for the super learner
	vv <- 1
	for(bb in DATA.split){
		tempLearn <- X[-(bb), , drop=FALSE]
		tempOutcome <- N.delta[-(bb)]
		tempTime <- time[-(bb)]
		tempValid <- X[bb, , drop=FALSE]
		tempid <- id[-(bb)]
		tempobsWeights <- obsWeights[-(bb)]
		tempLearn.is <- incidenceSample(long.DATA=data.frame(N.delta=tempOutcome, delta.u=tempTime, tempid = tempid, tempobsWeights = tempobsWeights, tempLearn), n.controls=n.controls)
		
		tempWhichScreen <- matrix(NA, nrow = k.screen, ncol = p)
		for(s in seq(k.screen)){
			testScreen <- try(do.call(library$screenAlgorithm[s], list(Y.temp = tempLearn.is$N.delta, X.temp = tempLearn.is[, -c(1:4)], family = family, id = tempLearn.is$tempid, obsWeights = tempLearn.is$tempobsWeights)))
			if(inherits(testScreen, "try-error")) {
				warning(paste("replacing failed screening,", library$screenAlgorithm[s], ", algorithm with All() in fold", vv ,"\n ")) 
				tempWhichScreen[s, ] <- TRUE
				next
			} else {
				tempWhichScreen[s, ] <- testScreen
			}
		}
		
		for(k in 1:k.cand){
			testAlg <- try(do.call(library$predAlgorithm[k, 1], list(Y.temp = tempLearn.is$N.delta, X.temp = subset(tempLearn.is[, -c(1:4)], select = tempWhichScreen[library$predAlgorithm[k, 2], ], drop=FALSE), newX.temp = subset(X, select = tempWhichScreen[library$predAlgorithm[k, 2], ], drop=FALSE), family = family, id = tempLearn.is$tempid, obsWeights = tempLearn.is$tempobsWeights)))
			if(inherits(testAlg, "try-error")) {
				warning(paste("Error in algorithm", library$predAlgorithm[k, 1], " on fold", vv, "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" )) 
				errorsInCVLibrary[k] <- 1
				next
			} else {
				newZ.X[, k] <- testAlg$out
			}
			fit.gam.data <- data.frame(N.delta=tempOutcome, time=tempTime, newz.x = trimLogit(newZ.X[-(bb), k], trim=0.0001))
			fit.gam.model <- as.formula(paste("N.delta~s(time,", time.df, ") + offset(newz.x)" ))
			fit.gam.t <- gam(fit.gam.model, family=binomial, data=fit.gam.data)
			newZ[bb,k] <- predict(fit.gam.t, newdata=data.frame(time=time[bb], newz.x = trimLogit(newZ.X[bb, k], trim=0.0001)), type="response")
			if(verbose) print(paste("CV", cand.names[k]))
		}#end library
		if(verbose) print(paste("V-fold:", vv))
		vv <- vv + 1
	}#end V-fold


	# now estimte the weights for each algorithm in library
	if(sum(errorsInCVLibrary) > 0) {
		newZ[, as.logical(errorsInCVLibrary)] <- 0 
	}
	if(all(newZ==0)) {
		stop("All algorithms dropped from library")
	}
	getweights <- .computeCoef(newZ=newZ, Y=N.delta, cand.names=cand.names, method=method, trim=trim.logit, verbose=verbose, obsWeights = obsWeights)
	init.coef <- getweights$init.coef
	update.coef <- getweights$update.coef
	names(update.coef) <- cand.names
	names(init.coef) <- cand.names
		
	## now fit all candidate on entire training set and predict on newX
	m <- nrow(newX)
	predY <- matrix(NA, m, k.cand)	
	predY.X <- matrix(NA, N.all, k.cand)
	whichScreen <- matrix(NA, nrow = k.screen, ncol = p)
	
	Train.DATA <- data.frame(N.delta=N.delta, delta.u=time, id=id, obsWeights=obsWeights, X)
	Test.DATA <- data.frame(N.delta=NA, delta.u=time.newX, newX)
	Train.is <- incidenceSample(Train.DATA, n.controls=n.controls)
	
	for(s in seq(k.screen)) {
		testScreen <- try(do.call(library$screenAlgorithm[s], list(Y.temp = Train.is$N.delta, X.temp = Train.is[, -c(1:4)], family = family, id = Train.is$id, obsWeights = Train.is$obsWeights)))
		if(inherits(testScreen, "try-error")) {
			warning(paste("replacing failed screening,", library$screenAlgorithm[s], ", algorithm with All() in full data", "\n ")) 
			whichScreen[s, ] <- TRUE
			next
		} else {
			whichScreen[s, ] <- testScreen
		}
	}
	
	for(k in 1:k.cand){
		fitProb <- try(do.call(library$predAlgorithm[k, 1], list(Y.temp=Train.is$N.delta, X.temp = subset(Train.is[, -c(1:4)], select=whichScreen[library$predAlgorithm[k, 2], ], drop=FALSE), newX.temp=subset(rbind(X, newX), select=whichScreen[library$predAlgorithm[k, 2], ], drop=FALSE), family=family, id = Train.is$id, obsWeights = Train.is$obsWeights)))
		if(inherits(fitProb, "try-error")) {
			warning(paste("Error in algorithm", library$predAlgorithm[k, 1], " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" )) 
			errorsInLibrary[k] <- 1
			next
		}
				
		fit.gam.data <- data.frame(N.delta=N.delta, time=time, predy.x=trimLogit(fitProb$out[1:N.all], trim=0.0001))
		fit.gam.model <- as.formula(paste("N.delta~s(time,", time.df, ") + offset(predy.x)" ))
		fit.gam.time <- gam(fit.gam.model, family=binomial, data=fit.gam.data)
		predY[,k] <- predict(fit.gam.time, newdata=data.frame(time=time.newX, predy.x = trimLogit(fitProb$out[-c(1:nrow(X))], trim=0.0001)), type="response")
		if(save.fit.library) {
			fit.library[[k]] <- fitProb$fit
			fit.gam.library[[k]] <- fit.gam.time
		}
		if(verbose) {
			print(paste("full", cand.names[k]))
		}
	}#end library	
	
	if(sum(errorsInLibrary) > 0) {
		if(sum(update.coef[as.logical(errorsInLibrary)]) > 0) {
			warning(paste("re-running estimation of coefficients removing failed algorithm(s) \n Orignial coefficients are: \n", update.coef, "\n"))
			newZ[, as.logical(errorsInLibrary)] <- 0
			if(all(newZ==0)) {
				stop("All algorithms dropped from library")
			}
			getweights <- .computeCoef(newZ=newZ, Y=N.delta, cand.names=cand.names, method=method, trim=trim.logit, verbose=verbose, obsWeights = obsWeights)
			init.coef <- getweights$init.coef
			update.coef <- getweights$update.coef
			names(update.coef) <- cand.names
			names(init.coef) <- cand.names
		} else {
			warning(" coefficients already 0 for all failed algorithm(s)")
		}
	}
	
	## super learner predictions
	getPred <- .computePred(predY=predY, init.coef=init.coef, update.coef=update.coef, method=method, trim=trim.logit)
	
	final <- list(call=call, cand.names=cand.names, SL.library=library, SL.predict=getPred$SL.pred.update, init.coef=init.coef, coef=update.coef, library.predict=predY, newZ=newZ, cv.risk=getweights$cv.risk, family=family, fit.library=fit.library, id=id, obsWeights = obsWeights, namesX=namesX, DATA.split=DATA.split, fit.gam.library=fit.gam.library, method=method, whichScreen=whichScreen, trim.logit=trim.logit, errorsInCVLibrary = errorsInCVLibrary, errorsInLibrary = errorsInLibrary)
	class(final) <- c("SuperLearner")
	return(final)
}



# test
# test <- SuperLearner.hazard(N.delta=long.DATA$N.delta,time=long.DATA$delta.u,X=long.DATA[,c(4:8)],V=10,verbose=TRUE,SL.library=c("SL.glm"),id=long.DATA$ID)