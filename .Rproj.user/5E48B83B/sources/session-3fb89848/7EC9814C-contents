
# Run cross-validation on MSOMs


setwd("~/testRenvProject/Source_Files")  # Fill in as appropriate
# setwd("~/ArkCrossVal")

### Required packages:
# library(animation)  # for kfcv function
# library(boot)  # for logit functions
# library(R2jags)  # to run the models. needs to have Jags downloaded on cmputer. 
# library(mcmcplots)  # to check model convergence
# library(verification)  # to calc AUC statistic (used in functions)
# library(matrixStats)

# library(doParallel)  # to run folds in parallel
# library(foreach)  # to run folds in parallel
# library(parallel) # for mclapply
# cl <- makeCluster(n.folds)
# registerDoParallel(cl)

### Import the data and source functions -----------------------
#base::load("~/testRenvProject/Output_Files/2_OrganizeNewData/ArkData.Rdata")
#source("~/testRenvProject/Source_Files/ArkFunctions.R")



### Set up for all models, all cross-validation data sets ---------------
set.seed(2000)

# Number of coef/covariates of each model:
nPsiList <- c(17, 17, 15, 9, 7, 10, 8, 9, 6, 10, 10,
              7, 6, 9, 2, 2, 5, 5, 5, 5, 2, 2, 2,
              2, 2, 2, 5, 3)
nDetList <- c(28, rep(12, 3), 13, 11, rep(10, 4), 11,
              13, 12, 12, 12, 10, 12, 10, 13, 10, 11, 4, 11,
              4, 3, 4, 3, 4)
nModels <- length(nPsiList)

# Tell JAGS what variables to save for output
params <- c("muPsiCoef", "sigmaPsiCoef", 
            "muDetCoef", "sigmaDetCoef", 
            "detCoef", "psiCoef", "score01ij")

# Convert the SDs in the prior distributions to precisions for JAGs
muPsiCoef.tau <- 1 / 2.25 ^ 2
muDetCoef.tau <- 1 / 2.25 ^ 2
t.tau <- 1 / 2.25 ^ 2

# Tell JAGS what variables in the model represent data
jags.data <- list('Y', #'Ybinom',
                  'muPsiCoef.tau', 'muDetCoef.tau', 't.tau', # params for prior distributions
                  'nSites', 'nSurveys', 'nPossSpp', 
                  'nPsiCoef', 'nDetCoef',
                  'site.vars', 'X.array') 
# Y = 3d array of observations
# Ybinom= 2d array = sum of detections for each species and site
# Psi.used = Subset of matrix of site-level covariates-- changes for each model
# nPossSpp = 17 known native fish in the Ark
# X.array = survey-level covariates

# Provide initial values for each cross-val run
# (Because of the way jags.parallel works with do.call, we need the two functions)
jags.inits113 <- function(nPossSpp=17, nSites=113){
  list("muPsiCoef"=rnorm(nPsiCoef, 0, 1), 
       "muDetCoef"=rnorm(nDetCoef, 0, 1),
       "halfsigmaPsiCoef"=runif(nPsiCoef, 0, 2), 
       "halfsigmaDetCoef"=runif(nDetCoef, 0, 2),
       "Z"=matrix(1, nrow=17, ncol=113))
}
jags.inits112 <- function(nPossSpp=17, nSites=112){
  list("muPsiCoef"=rnorm(nPsiCoef, 0, 1), 
       "muDetCoef"=rnorm(nDetCoef, 0, 1),
       "halfsigmaPsiCoef"=runif(nPsiCoef, 0, 2), 
       "halfsigmaDetCoef"=runif(nDetCoef, 0, 2),
       "Z"=matrix(1, nrow=17, ncol=112))
}

# nChains; nIterations; nBurn-in; nThin
nc=3; 
ni = 15000 # 50 #
nb = 5000   # 
nt = 10  #

### Specific cross-validation set up ----------------------

n.folds <- 5
nSites.all <- nSites

## Randomize data order for CV tests
set.seed(2015)
kf = cumsum(c(1, kfcv(n.folds, nSites)))
randomized.sites <- sample(nSites)
cvY <- Y[, randomized.sites, ]
cvYbinom <- Ybinom[, randomized.sites]
cv.nSurveys <- nSurveys[randomized.sites]

cvSiteVars <- site.varsAll[, randomized.sites]
cvX.arrayAll <- X.arrayAll[, randomized.sites, ]

# Create list of model files:
models <- list()
models[[1]] <- "JAGS_models/KnownSppMSOM.Full.JAGS.R"
models[[2]] <- "JAGS_models/KnownSppMSOM.FullOccu.JAGS.R"
models[[3]] <- "JAGS_models/KnownSppMSOM.FullOccu2.JAGS.R"
models[[4]] <- "JAGS_models/KnownSppMSOM.NoCorrs.JAGS.R"
models[[5]] <- "JAGS_models/KnownSppMSOM.LandCoverSize.JAGS.R"
models[[6]] <- "JAGS_models/KnownSppMSOM.Combined.JAGS.R"
models[[7]] <- "JAGS_models/KnownSppMSOM.Reservoirs.JAGS.R"
models[[8]] <- "JAGS_models/KnownSppMSOM.Stream.JAGS.R"
models[[9]] <- "JAGS_models/KnownSppMSOM.Size.JAGS.R"
models[[10]] <- "JAGS_models/KnownSppMSOM.NoCorrs2.JAGS.R"
models[[11]] <- "JAGS_models/KnownSppMSOM.Combined2.JAGS.R"

models[[12]] <- "JAGS_models/KnownSppMSOM.LandCoverSize.JAGS.R"
models[[13]] <- "JAGS_models/KnownSppMSOM.AltSize.JAGS.R"
models[[14]] <- "JAGS_models/KnownSppMSOM.AltStream.JAGS.R"
models[[15]] <- "JAGS_models/KnownSppMSOM.SimpleSize.JAGS.R"
models[[16]] <- "JAGS_models/KnownSppMSOM.SimpleSize2.JAGS.R"
models[[17]] <- "JAGS_models/KnownSppMSOM.SimpleAltStream.JAGS.R"
models[[18]] <- "JAGS_models/KnownSppMSOM.SimpleAltStream2.JAGS.R"
models[[19]] <- "JAGS_models/KnownSppMSOM.SimpleLandCoverSize.JAGS.R"
models[[20]] <- "JAGS_models/KnownSppMSOM.SimpleLandCoverSize2.JAGS.R"

models[[21]] <- "JAGS_models/KnownSppMSOM.SimpleSizeCount.JAGS.R"
models[[22]] <- "JAGS_models/KnownSppMSOM.SimpleSizeCountNoYear.JAGS.R"
models[[23]] <- "JAGS_models/KnownSppMSOM.SimpleSize2.JAGS.R"

models[[24]] <- "JAGS_models/KnownSppMSOM.SimpleSizeCountNoYear.JAGS.R"
models[[25]] <- "JAGS_models/KnownSppMSOM.SimpleSize2noYr.JAGS.R"
models[[26]] <- "JAGS_models/KnownSppMSOM.SimpleSizeNoYr.JAGS.R"
models[[27]] <- "JAGS_models/KnownSppMSOM.SimpleAltStream2noYr.JAGS.R"
models[[28]] <- "JAGS_models/KnownSppMSOM.SimpleSizeNoYrConnect.JAGS.R"


### Run the cross-validation -----------------
Probs <- list()
cv.stats <- array(NA, dim=c(n.folds, nModels, 4))
out.list <- list()  # save all output.
## Create probability expressions
tmp <- createExpressionsArk()
expressPsi <- tmp$expressPsi
expressDet <- tmp$expressDet

for (fold in 1:n.folds){
# out.list <-  mclapply(1:n.folds, runCV, mc.cores=5)
# #out.list <-  foreach(fold = 1:n.folds, .packages="R2jags") %dopar% {
#   runCV <- function(){
  Probs[[fold]] <- list()
  test <- c(kf[fold]:(kf[fold + 1]-1) ) # test = "testing data"
  cv.nSurveys.test <- cv.nSurveys[test]
  testY <- cvY[, test, 1:max(cv.nSurveys.test)]
  
  ## save output in lists for organization:
  out <- list()
  
  # Use a different inits function depending on number of sites in training data
  if (nSites.all - length(test) == 113){  
    jags.inits <- jags.inits113
  }else{
  	  jags.inits <- jags.inits112
  }
  
  # Fit the models to the training data set --------
  # First re-create the data:
  Y            <- cvY[, -test, ]; 
  Ybinom       <- cvYbinom[, -test]; 
  nSites       <- nSites.all - length(test);
  nSurveys     <- cv.nSurveys[-test];
  site.vars    <- cvSiteVars[, -test];
  X.array      <- cvX.arrayAll[, -test, ]

  ### Fit Models in a loop 
  cat("Starting fold", fold, "\n")
  cat("Y dim = ", dim(Y), "\n")
  # print(jags.inits, "\n")
  for (i in 24:nModels){
    nPsiCoef <- nPsiList[[i]]  # elev + intecept
    nDetCoef <- nDetList[[i]] # year, elev
    out[[i]] <- do.call(jags.parallel, list(jags.data, jags.inits, params, 
                                            models[[i]],
                                            nc, ni, nb, nt));
#     out[[i]] <- jags(jags.data, jags.inits, params, 
#                      models[[i]], nc, ni, nb, nt);
     cat("Mod ", i, "done!\n")
  }
#  return(out)
  #save.image("~/testRenvProject/Output_Files/3_RunCrossVal/CrossValProbs2.RData") # save often in case it crashes
  
  # Process the output:
  for (i in 24:nModels) {  
    Probs[[fold]][[i]] <- derivePsiDetect(out[[i]], 
                                          expressPsi=expressPsi[[i]], 
                                          expressDet=expressDet[[i]], 
                                          nPossSpp=nPossSpp, 
                                          nSites=length(test), 
                                          nSurveys=cv.nSurveys.test )
    
    score01 <- mean( apply(out[[i]]$BUGSoutput$sims.list$score01ij, 
                           c(2, 3), sum, na.rm=T) )
    other.scores <- calcCrossVal(testY, Probs[[fold]][[i]], 
                                 nPossSpp=nPossSpp, nSites=length(test),              
                                 nSurveys=cv.nSurveys.test)
    cv.stats[fold, i, ] <- c(score01, other.scores )
    cat("0-1 score \t Deviance \t AUC \t Brier \n")
    print(round(cv.stats[fold, i, ], 2))   
    write.table(t(c(paste("fold", fold, sep=""), 
                    paste("mod", i, sep=""), 
                    cv.stats[fold, i, ])), 
                "~/testRenvProject/Output_Files/3_RunCrossVal/CrossValStats2.csv", append=T, sep=",", 
                row.names=FALSE, col.names=FALSE)
  }
  out.list[[i]] <- out
}

#save.image("~/testRenvProject/Output_Files/3_RunCrossVal/CrossValProbs2.RData") # save often in case it crashes


### ------------------
###  Calculate cross-validation statistics:
###            0-1 score, log score, Brier score, AUC
### ------------------

### And Summarize: ------------------------

# Cross-val stats are from saved .csv file.
crossVal.df <- read.csv("~/testRenvProject/Output_Files/3_RunCrossVal/CrossValStats2.csv", header=F)
### Table of CV values for all folds and all models: (Table C.1 in Appendix C)
# head(crossVal.df)

## Take the average over the folds to get the reported statistics (Table 2)
crossVal.means <- aggregate(crossVal.df[, -c(1, 2)], by=list(Model = crossVal.df$V2), mean)
names(crossVal.means)[2:5] <- c("score01", "deviance", "AUC", "Brier")
# crossVal.means[order(crossVal.means$deviance), ]
#

write.csv(crossVal.means, "~/testRenvProject/Output_Files/3_RunCrossVal/CrossValMeans2.csv", row.names = F)
#


### end of file. -------------------------------

