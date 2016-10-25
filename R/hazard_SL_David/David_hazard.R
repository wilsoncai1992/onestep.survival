## makeDataList
# this function takes as input:
# dat = a data.frame with columns named:
#         id = unique subject id
#         T.tilde = (possibly censored) failure time
#         Delta.J= (possible censored) failure type -- can set all to 1 if no censoring 
#         Z = binary treatment indicator
#         W = other columns for baseline variables
# J = the unique positive values of dat$Delta.J (set J=1 if only one type of failure)
# nZ = set to length of unique(dat$Z)
# Z = the treatment values that are of interest (e.g., Z=c(0,1))
# t0 = the time point you eventually want predictions at
# bounds = just leave NULL

# output:
# list[[1]]: orginal data into stacked form
# list[[2]]: data to predict on, set Z to 0
# list[[3]]: data to predict on, set Z to 1

# returns a list of length nZ + 1. first entry is data.frame used by SL
# next nZ entries are data.frames where the column Z has been set to a value 
# given by the unique values of the input Z; for these data.frames each subject
# has t0 rows (compared to T.tilde rows in the first entry). This allows for prediction
# of the survival function at t0.

# each data frame will have columns:
#   - id = unique subject id
#   - T.tilde = same as before
#   - Delta.J = same as before
#   - Z = same as before
#   - Nj = indicator of failure of type j at time t, one for each j \in J
#   - C = indicator of censoring at time t
#   - t = time
#   - lj,uj = for the bounds option -- just ignore

# example usage:
# n <- 20
# set.seed(1234)
# dat <- data.frame(
#   id = 1:n,
#   W1 = runif(n), W2=runif(n),
#   T.tilde = rpois(n,10),
#   Delta.J = rbinom(n,1,0.9),
#   Z = rbinom(n,1,0.5)
# )
# 
# test <- makeDataList(
#   dat=dat, J=1, nZ=2, Z=c(0,1), t0=5, bounds=NULL
# )
# 
# # check length
# length(test)
# 
# # look at it
# head(test[[1]],10)
# head(test[[2]],10)
# head(test[[3]],10)
# 
# # first df has nrow = T.tilde for each subject
# nrow(test[[1]][test[[1]]$id==1,]) == dat$T.tilde[dat$id==1]
# 
# # second/third df has nrow = t0
# nrow(test[[2]][test[[2]]$id==1,]) == 5 
# nrow(test[[3]][test[[3]]$id==1,]) == 5 
# 
# # first df has Z = observed treatment
# table(test[[1]]$Z)
# 
# # second/third are set to the values of Z given to makeDataList
# table(test[[2]]$Z)
# table(test[[3]]$Z)

makeDataList <- function(dat, J, nZ, Z, t0, bounds=NULL){
  n <- nrow(dat)
  dataList <- vector(mode="list",length=nZ+1)
  
  # first element used for estimation
  dataList[[1]] <- dat[rep(1:nrow(dat),dat$T.tilde),]
  for(j in J){
    eval(parse(text=paste("dataList[[1]]$N",j," <- 0",sep="")))
    eval(parse(text=paste("dataList[[1]]$N",j,"[cumsum(dat$T.tilde)] <- as.numeric(dat$Delta.J==j)",sep="")))  
  }
  dataList[[1]]$C <- 0
  dataList[[1]]$C[cumsum(dat$T.tilde)] <- as.numeric(dat$Delta.J==0)
  
  n.row.ii <- nrow(dataList[[1]])
  row.names(dataList[[1]])[row.names(dataList[[1]]) %in% paste(row.names(dat))] <- paste(row.names(dat),".0",sep="")
  dataList[[1]]$t <- as.numeric(paste(unlist(strsplit(row.names(dataList[[1]]),".",fixed=T))[seq(2,n.row.ii*2,2)]))+1
  
  if(!is.null(bounds)){
    boundFormat <- data.frame(t=bounds$t)
    for(j in J){
      if(paste("l",j,sep="") %in% names(bounds)){
        eval(parse(text=paste("boundFormat$l",j," <- bounds$l",j,sep="")))
      }else{
        eval(parse(text=paste("boundFormat$l",j," <- 0",sep="")))
      }
      if(paste("u",j,sep="") %in% names(bounds)){
        eval(parse(text=paste("boundFormat$u",j," <- bounds$u",j,sep="")))
      }else{
        eval(parse(text=paste("boundFormat$u",j," <- 1",sep="")))
      }
    }
    suppressMessages(
      dataList[[1]] <- join(x=dataList[[1]],y=boundFormat,type="left")
    )  
  }else{
    for(j in J){
      eval(parse(text=paste("dataList[[1]]$l",j," <- 0",sep="")))
      eval(parse(text=paste("dataList[[1]]$u",j," <- 1",sep="")))
    }
  }
  
  # subsequent elements used for prediction
  for(i in 1:nZ){
    dataList[[i+1]] <- dat[sort(rep(1:nrow(dat),t0)),]
    dataList[[i+1]]$t <- rep(1:t0,n)
    for(j in J){
      typejEvents <- dat$id[which(dat$Delta.J==j)]
      eval(parse(text=paste("dataList[[i+1]]$N",j," <- 0",sep="")))
      eval(parse(text=paste("dataList[[i+1]]$N",j,"[dataList[[i+1]]$id %in% typejEvents &  dataList[[i+1]]$t >= dataList[[i+1]]$T.tilde] <- 1",sep="")))
    }
    censEvents <- dat$id[which(dat$Delta.J==0)]
    dataList[[i+1]]$C <- 0
    dataList[[i+1]]$C[dataList[[i+1]]$id %in% censEvents & dataList[[i+1]]$t >= dataList[[i+1]]$T.tilde] <- 1
    dataList[[i+1]]$Z <- Z[i]
    dataList[[i+1]]$T.tilde <- t0
    
    if(!is.null(bounds)){
      suppressMessages(
        dataList[[i+1]] <- join(x=dataList[[i+1]],y=boundFormat,type="left")
      ) 
    }else{
      for(j in J){
        eval(parse(text=paste("dataList[[",i,"+1]]$l",j," <- 0",sep="")))
        eval(parse(text=paste("dataList[[",i,"+1]]$u",j," <- 1",sep="")))
      }
    }
  }
  names(dataList) <- c("obs",Z)
  return(dataList)
}

## estimate hazards using Super Learning
# function for computing hazards using Super Learner or GLM
# takes as input:
#   dataList = output of makeDataList
#   J = same as above
#   verbose = print messages as SL runs?
#   strata = would leave as NULL, dont remember exactly how it is coded
#   adjustVars = data.frame of predictors for SL, dont know why I have it coded this way...
#                 just put in e.g. dat[,c("W1","W2")]. by default it will add the Z column
#                 and t column in dataList[[1]] when it fits the SL
#   SLlibrary.event = Super Learner library
#   glmFormula.event = a glm formula if you don't want to use SL
#   superLearnerSummary = doesn't do anything, just leave blank
#   bounds = ignore

# output:
# added column for the initial data.frame
# QjHaz <- result, for each subject, h_t from time 0 to T-1; organized in T_i rows
# QjPseudoHaz <- 
# returns list of length = length(dataList) with added columns name QjPseudoHaz and 
# QjHaz for all j \in J. Ignore the PseudoHaz variables -- they're used when there are 
# multiple event types; if there's only one type of event they'll == Haz values anyway

# example useage (continued from above)
# testHaz <- estimateHazards(
#   dataList = test, J=1, verbose=FALSE, strata=NULL, 
#   adjustVars = dat[,c("W1","W2")], # stupid, need to recode this 
#   SLlibrary.event = c("SL.mean","SL.glm"),
#   glmFormula.event = NULL,
#   bounds = NULL
# )
# head(testHaz[[1]],10)

# in case you want to use GLM
# testHazGLM <- estimateHazards(
#     dataList = test, J=1, verbose=FALSE, strata=NULL,
#     adjustVars = dat[,c("W1","W2")], # stupid, need to recode this
#     SLlibrary.event = NULL,
#     glmFormula.event = "Z + t + W1*W2", # columns Z and t added to adjustVars in fitting procedures
#     bounds = NULL
#   )

estimateHazards <- function(dataList, J, verbose,strata=NULL,adjustVars,
                            SLlibrary.event, glmFormula.event,
                            superLearnerSummary, bounds){
  # check for missing inputs
  if(is.null(SLlibrary.event) & is.null(glmFormula.event) & is.null(strata)){
    warning("Super Learner library, glm formula, and strata for events not specified. Proceeding 
            with empirical estimates")
    glmFormula.event <- "Z*factor(t)"
  }
  
  if(is.null(SLlibrary.event) & is.null(strata)){
    if(is.null(bounds)){
      for(j in J){
        # formula
        Qj.form <- sprintf("%s ~ %s", paste("N",j,sep=""), glmFormula.event)
        
        # add up all events less than current j to see who to include in regression
        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for(i in J[J<j]){
          eval(parse(text=paste("NlessthanJ <- NlessthanJ + dataList[[1]]$N",i,sep="")))
        }
        
        # fit glm
        Qj.mod <- glm(as.formula(Qj.form), data=dataList[[1]][NlessthanJ==0,], family="binomial")
        
        # get predictions back
        dataList <- lapply(dataList, function(x,j){
          eval(parse(text=paste("x$Q",j,"PseudoHaz <- predict(Qj.mod, type='response', newdata=x)",sep="")))
          if(j != min(J)){
            eval(parse(text=paste("x$hazLessThan",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',J[J<j],'Haz')]))",sep="")))
            eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz * (1-x$hazLessThan",j,")",sep="")))
          }else{
            eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz",sep="")))
          }
          x
        }, j=j)
      }
    }else{
      for(j in J){
        Qj.form <- sprintf("%s ~ %s", paste("N",j,sep=""), glmFormula.event)
        X <- model.matrix(as.formula(Qj.form),data=dataList[[1]])
        
        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for(i in J[J<j]){
          eval(parse(text=paste("NlessthanJ <- NlessthanJ + dataList[[1]]$N",i,sep="")))
        }
        
        dataList <- lapply(dataList, function(x,j){
          if(j != min(J)){
            eval(parse(text=paste("x$hazLessThan",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',J[J<j],'Haz')]))",sep="")))
          }else{
            eval(parse(text=paste("x$hazLessThan",j," <- 0",sep="")))
          }
          x
        },j=j)
        
        eval(parse(text=paste("Ytilde <- (dataList[[1]]$N",j,"-dataList[[1]]$l",j,")/(pmin(dataList[[1]]$u",j,", 1 - dataList[[1]]$hazLessThan",j,")  - dataList[[1]]$l",j,")",sep="")))
        fm <- optim(par=rep(0,ncol(X)), fn=LogLikelihood, Y=Ytilde, X=X, 
                    method="BFGS",gr=grad,
                    control=list(reltol=tol/1e4))
        if(fm$convergence!=0){
          return("convergence failure")
        }else{
          beta <- fm$par
          
          dataList <- lapply(dataList, function(x,j){
            newX <- model.matrix(as.formula(Qj.form),data=x)
            eval(parse(text=paste("x$Q",j,"PseudoHaz <- plogis(newX%*%beta)",sep="")))
            eval(parse(text=paste("x$Q",j,"Haz <- (pmin(x$u",j,", 1 - x$hazLessThan",j,")  - x$l",j,") * x$Q",j,"PseudoHaz + x$l",j,sep="")))
            x
          },j=j)
        }
      }
    }
  }else if(is.null(glmFormula.event) & is.null(strata)){
    for(j in J){
      # add up all events less than current j to see who to include in regression
      NlessthanJ <- rep(0, nrow(dataList[[1]]))
      for(i in J[J<j]){
        eval(parse(text=paste("NlessthanJ <- NlessthanJ + dataList[[1]]$N",i,sep="")))
      }
      
      Qj.mod <- eval(parse(text=paste("SuperLearner(Y=dataList[[1]]$N",j,"[NlessthanJ==0],
                             X=dataList[[1]][NlessthanJ==0,c('t', 'Z', names(adjustVars))],
                             id=dataList[[1]]$id[NlessthanJ==0],
                             family=binomial(),
                             SL.library=SLlibrary.event,
                             verbose=verbose)",sep="")))
      
      # get predictions back
      dataList <- lapply(dataList, function(x,j){
        eval(parse(text=paste("x$Q",j,"PseudoHaz <- predict(Qj.mod, onlySL=T, newdata=x)[[1]]",sep="")))
        if(j != min(J)){
          eval(parse(text=paste("x$hazLessThan",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',J[J<j],'Haz')]))",sep="")))
          eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz * (1-x$hazLessThan",j,")",sep="")))
        }else{
          eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz",sep="")))
        }
        x
      }, j=j)
    }
  }else if(is.null(glmFormula.event) & is.null(SLlibrary.event)){
    for(j in J){
      X <- model.matrix(~factor(dataList[[1]]$t)*dataList[[1]]$Z*factor(dataList[[1]]$strata))
      eval(parse(text=paste("beta <- solve(t(X)%*%X)%*%t(X)%*%dataList[[1]]$N",j,sep="")))
      dataList <- lapply(dataList, function(x,beta){
        X <- model.matrix(~factor(x$t)*x$Z*factor(x$strata))
        eval(parse(text=paste("x$Q",j,"Haz <- X%*%beta",sep="")))
        x
      },beta=beta)
    }
  }
  dataList
  
}



