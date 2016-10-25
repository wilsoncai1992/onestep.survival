# source("./Hazard_SL/createDiscrete.R")
# source("./Hazard_SL/SuperLearner.hazard.R")
# source("./Hazard_SL/createLibrary.R")
# source("./Hazard_SL/checkSL.R")
# source("./Hazard_SL/CVFolds.R")
# source("./Hazard_SL/incidenceSample.R")

library(SuperLearnerOld)

event <- (dat$T[1] == T.uniq) + 0
createDiscrete(time = T.uniq, event = event, dataX = W)
createDiscrete(time = c(1, 2), event = c(0, 1), dataX = W)


SL.lib <- c("SL.glm", "SL.step", "SL.glm.interaction" )
SuperLearner.hazard(N.delta = event, n.controls = 0, time = T.uniq, X = as.matrix(rep(W[1,], length(event))), SL.library = SL.lib)
# CV.SuperLearner.hazard(N.delta = event, time = T.uniq, X = W, SL.library = SL.lib)
