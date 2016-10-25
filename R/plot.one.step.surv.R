plot.onestep.surv <- function(onestepfit, col = 'green', lty = 1, add = FALSE) {
	if (add){
		lines(onestepfit[[1]] ~ onestepfit[[2]], col = col, lty = lty)
	}else{
		plot(onestepfit[[1]] ~ onestepfit[[2]], col = col, type = 'l', lty = lty)
	}

}