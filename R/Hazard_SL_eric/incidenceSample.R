incidenceSample <- function(long.DATA, n.controls=1) {
	# takes output from createDiscrete() and implements incidence sampling.  
	# naming structure from createDiscrete() must be intact
	if(n.controls == 0){
		return(long.DATA)
	}
	if(!('N.delta' %in% names(long.DATA))) {
	  stop('N.delta not found, long.DATA is not in the correct format')
	}
  # cases <- subset(long.DATA, N.delta == 1)
  cases <- long.DATA[long.DATA$N.delta == 1, ]
	n.cases <- nrow(cases)
	controls <- new("matrix", nrow = (n.controls*n.cases), ncol = ncol(long.DATA))
	controls <- as.data.frame(controls)
	names(controls) <- names(long.DATA)

	for(ii in seq(n.cases)){
		risk.set <- long.DATA[(long.DATA$delta.upper == cases$delta.upper[ii] & long.DATA$N.delta==0), ]
		n.risk.set <- nrow(risk.set)
		if(n.risk.set < n.controls) { print(paste("insufficient number of controls at case number",ii)) ; next}
		control.rows <- sample(seq(n.risk.set),size=n.controls)
		controls[((ii*n.controls - (n.controls-1)):(ii*n.controls)), ] <- risk.set[control.rows, ]
	}

	combi.sample <- rbind(cases,controls)
	combi.sample <- combi.sample[complete.cases(combi.sample), ]
	return(combi.sample)
}
