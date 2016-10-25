# 
#  createLibrary.R
#  SuperLearner
#  
#  Created by Eric Polley on 2009-09-29.
# 

# takes the input from SL.library and sets up the library for SuperLearner
# SL.library may be a character vector or a list
# return list of prediction algorithms and list of screening algorithms.  If screening algorithms are used, must generate the whichScreen matrix and assign the appropriate row number to the prediction algorithm.  whichScreen[1, ] will always be the full data set, even if not used.

.createLibrary <- function(SL.library) {
	if (is.character(SL.library)) { 
		k.cand <- length(SL.library)
		whichScreen <- matrix(1, nrow = 1, ncol = k.cand)
		screenAlgorithm <- "All"
		predAlgorithm <- data.frame(library = SL.library, row = 1, stringsAsFactors=FALSE)
	} else {
		# LibraryNames <- laply(SL.library, .fun = "[", 1)
		LibraryNames <- sapply(SL.library, FUN="[", 1)
		# NumberScreen <- (laply(SL.library, .fun = length) - 1)
		NumberScreen <- (sapply(SL.library, FUN=length) - 1)
		if (sum(NumberScreen == 0) > 0) {
			for(ii in which(NumberScreen == 0)){
				SL.library[[ii]] <- c(SL.library[[ii]], "All")
				NumberScreen[ii] <- 1
			}
		}
		screenAlgorithmFull <- unlist(sapply(SL.library, FUN="[", -1))
		screenAlgorithm <- unique(screenAlgorithmFull)
		
		predAlgorithm <- data.frame(library = rep(LibraryNames, times=NumberScreen), row = match(screenAlgorithmFull, screenAlgorithm), stringsAsFactors = FALSE)
	}
	
	out <- list(predAlgorithm = predAlgorithm, screenAlgorithm = screenAlgorithm)
	return(out)
}


# # examples:
# # just character vector
# SL.library <- c("SL.glm", "SL.foo")
# 
# SL.library <- list("SL.glm", "SL.foo")
# # including screening
# SL.library <- list("SL.glm", c("SL.foo", "screen.lm", "screen.rf"), "SL.rF", c("SL.els", "All", "screen.lm", "screen.rf"))
# 
# SL.library <- list(c("SL1", "All"), c("SL2", "All"), c("SL3", "screen2"))
# 
# 
# # screening template
# All <- function(X.temp, ...) {
# 	return(TRUE)
# }
# 
# # screen functions must return a logical vector of length ncol(X)
# screen.template <-function (Y.temp, X.temp, family = gaussian(), ...) 
# {
#     if (family$family == "gaussian") {
# 	
#     }
#     if (family$family == "binomial") {
# 	
#     }
# 	whichVariable <- rep(TRUE, ncol(X.temp))
#     return(whichVariable)
# }
# 
# 
# # examples:
# n <- 100
# p <- 10
# X <- matrix(rnorm(n*p), n, p)
# X <- as.data.frame(X)
# Y <- rnorm(n)
# whichScreen <- matrix(NA, 2, p)
# whichScreen[1, ] <- TRUE
# whichScreen[2, ] <- rep(c(FALSE,TRUE), 5)
# 
# SL.glm <- function (Y.temp, X.temp, newX.temp, family = gaussian(), ...) 
# {
#     fit.glm <- glm(Y.temp ~ ., data = X.temp, family = family)
#     out <- predict(fit.glm, newdata = newX.temp, type = "response")
#     fit <- list(object = fit.glm)
#     foo <- list(out = out, fit = fit)
#     class(foo) <- c(SL.glm)
#     return(foo)
# }
# 
# SL.randomForest <- function (Y.temp, X.temp, newX.temp, ntree = 1000, family = gaussian(), ...) 
# {
#     tryCatch(require(randomForest), warning = function(...) {
#         stop("you have selected randomForest as a library algorithm but do not have the randomForest package installed")
#     })
#     if (family$family == "gaussian") {
#         fit.rf <- randomForest(Y.temp ~ ., data = X.temp, ntree = ntree, 
#             xtest = newX.temp, keep.forest = TRUE)
#         out <- fit.rf$test$predicted
#         fit <- list(object = fit.rf)
#     }
#     if (family$family == "binomial") {
#         fit.rf <- randomForest(y = as.factor(Y.temp), x = X.temp, 
#             ntree = ntree, xtest = newX.temp, keep.forest = TRUE)
#         out <- fit.rf$test$votes[, 2]
#         fit <- list(object = fit.rf)
#     }
#     foo <- list(out = out, fit = fit)
#     class(foo) <- c("SL.randomForest")
#     return(foo)
# }
#  
# do.call(SL.glm, list(Y.temp=Y, X.temp=X[, whichScreen[2, ]], newX.temp = X[, whichScreen[2, ]]))