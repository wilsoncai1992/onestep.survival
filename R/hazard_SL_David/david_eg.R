source('./David_hazard.R')
library(SuperLearner)
# ====================================================================================================
# make data
# ====================================================================================================
n <- 20
set.seed(1234)
dat <- data.frame(
  id = 1:n,
  W1 = runif(n), W2=runif(n),
  T.tilde = rpois(n,10),
  Delta.J = rbinom(n,1,0.9),
  Z = rbinom(n,1,0.5)
)

test <- makeDataList(
  dat=dat, J=1, nZ=2, Z=c(0,1), t0=5, bounds=NULL
)

# check length
length(test)

# look at it
head(test[[1]],10)
head(test[[2]],10)
head(test[[3]],10)

# first df has nrow = T.tilde for each subject
nrow(test[[1]][test[[1]]$id==1,]) == dat$T.tilde[dat$id==1]

# second/third df has nrow = t0
nrow(test[[2]][test[[2]]$id==1,]) == 5
nrow(test[[3]][test[[3]]$id==1,]) == 5

# first df has Z = observed treatment
table(test[[1]]$Z)

# second/third are set to the values of Z given to makeDataList
table(test[[2]]$Z)
table(test[[3]]$Z)

# ====================================================================================================
# SL hazard
# ====================================================================================================
testHaz <- estimateHazards(
  dataList = test, J=1, verbose=FALSE, strata=NULL,
  adjustVars = dat[,c("W1","W2")], # stupid, need to recode this
  SLlibrary.event = c("SL.mean","SL.glm"),
  glmFormula.event = NULL,
  bounds = NULL
)
head(testHaz[[1]],10)


# in case you want to use GLM
testHazGLM <- estimateHazards(
    dataList = test, J=1, verbose=FALSE, strata=NULL,
    adjustVars = dat[,c("W1","W2")], # stupid, need to recode this
    SLlibrary.event = NULL, # use GLM ONLY
    glmFormula.event = "Z + t + W1*W2", # columns Z and t added to adjustVars in fitting procedures
    bounds = NULL
  )