
###   This script imports data, runs all ten models, 
### and outputs in-sample model selection criterion.

## INPUTS:
## ArkData.RData


##OUTPUTS:
##ArkModels.RData

## Set directory:
setwd("~/testRenvProject/")  # Fill in as appropriate
# setwd("~/ArkModels")  # for server

## LIBRARIES:
# library(boot)  # for logit functions
# library(R2jags)  # to run the models
# library(mcmcplots)  # to check model convergence
# library(verification)  # to calc AUC statistic (used in functions)


### Import the data and source functions -----------------------

#base::load("~/testRenvProject/Output_Files/2_OrganizeNewData/ArkData.Rdata")
#source("~/testRenvProject/Source_Files/ArkFunctions.R")
# base::load("~/ArkCrossVal/ArkData.RData")  # for server
# source("~/ArkCrossVal/ArkFunctions.R")  # for server

### Covariate objects cheat sheet ---------------------

# site.vars[, 1] <- "X" 
# site.vars[, 2] <- "X2" 
# site.vars[, 3] <- "Y" 
# site.vars[, 4] <- "Y2" 
# site.vars[, 5] <- "ELEV" 
# site.vars[, 6] <- "ELEV2"
# site.vars[, 7] <- "CROPS" 
# site.vars[, 8] <- "DVLPD" 
# site.vars[, 9] <- "WTLNDS"
# site.vars[, 10] <- "UNCONNECT"
# site.vars[, 11] <- "SIZE"
# site.vars[, 12] <- "aboveP" 
# site.vars[, 13] <- "belowJM"  # between reservois=base group
# site.vars[, 14] <- "FTN" 
# site.vars[, 15] <- "MAIN" 
# site.vars[, 16] <- "PURG"
# 
# X.array[1, j, 1:K[j]] <- YDAY
# X.array[2, j, 1:K[j]] <- YDAY2
# X.array[3, j, 1:K[j]] <- ifelse(tmp$METHOD == "S", 1, 0)
# X.array[4, j, 1:K[j]] <- PASSNO
# X.array[5, j, 1:K[j]] <- ifelse(tmp$YEAR == 2008, 1, 0)
# X.array[6, j, 1:K[j]] <- ifelse(tmp$YEAR == 2009, 1, 0)
# X.array[7, j, 1:K[j]] <- ifelse(tmp$YEAR == 2010, 1, 0)
# X.array[8, j, 1:K[j]] <- ifelse(tmp$YEAR == 2011, 1, 0)
# X.array[9, j, 1:K[j]] <- ifelse(tmp$YEAR == 2012, 1, 0)
# X.array[10, j, 1:K[j]] <- ifelse(tmp$YEAR == 2013, 1, 0)
# X.array[11, j, 1:K[j]] <- ifelse(tmp$YEAR == 2014, 1, 0)
# X.array[12, j, 1:K[j]] <- ifelse(tmp$YEAR == 2015, 1, 0)
# X.array[13, j, 1:K[j]] <- log(tmp$TotalCt)

### Set up for all models ---------------
set.seed(3000)

# Number of coef/covariates of each model:
# Number of coef/covariates of each model:
nPsiList <- c(17, 17, 15, 9, 7, 10, 8, 9, 6, 10, 10,
              7, 6, 9, 2, 2, 5, 5, 5, 5, 2, 2, 2,
              2, 2, 2, 5, 3)
nDetList <- c(28, rep(12, 3), 13, 11, rep(10, 4), 11,
              13, 12, 12, 12, 10, 12, 10, 13, 10, 11, 4, 11,
              4, 3, 4, 3, 4)
# nPsiList <- 2 # c(2, 2, 2)
# nDetList <- 3 # c(11, 4, 10)
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
site.vars <- site.varsAll
X.array <- X.arrayAll

# Provide initial values for each model run
jags.inits <- function(nPossSpp=17, nSites=141){
  list("muPsiCoef"=rnorm(nPsiCoef, 0, 1), #2.25), 
       "muDetCoef"=rnorm(nDetCoef, 0, 1), #2.25),
       "halfsigmaPsiCoef"=runif(nPsiCoef, 0, 2), #5), 
       "halfsigmaDetCoef"=runif(nDetCoef, 0, 2), #5),
       "Z"=matrix(1, nrow=nPossSpp, ncol=nSites))
}

# nChains; nIterations; nBurn-in; nThin
nc=3; 
ni= 30000; # 30000 # 
nb= 10000;  # 10000  #  
nt= 10 # 10  # 100

# Create list of model files:
models <- list()
models[[1]] <- "Source_Files/JAGS_models/KnownSppMSOM.Full.JAGS.R"
models[[2]] <- "Source_Files/JAGS_models/KnownSppMSOM.FullOccu.JAGS.R"
models[[3]] <- "Source_Files/JAGS_models/KnownSppMSOM.FullOccu2.JAGS.R"
models[[4]] <- "Source_Files/JAGS_models/KnownSppMSOM.NoCorrs.JAGS.R"
models[[5]] <- "Source_Files/JAGS_models/KnownSppMSOM.LandCoverSize.JAGS.R"
models[[6]] <- "Source_Files/JAGS_models/KnownSppMSOM.MyGuess.JAGS.R"
models[[7]] <- "Source_Files/JAGS_models/KnownSppMSOM.Reservoirs.JAGS.R"
models[[8]] <- "Source_Files/JAGS_models/KnownSppMSOM.Stream.JAGS.R"
models[[9]] <- "Source_Files/JAGS_models/KnownSppMSOM.Size.JAGS.R"
models[[10]] <- "Source_Files/JAGS_models/KnownSppMSOM.NoCorrs2.JAGS.R"
models[[11]] <- "Source_Files/JAGS_models/KnownSppMSOM.Combined2.JAGS.R"

models[[12]] <- "Source_Files/JAGS_models/KnownSppMSOM.LandCoverSize.JAGS.R"
models[[13]] <- "Source_Files/JAGS_models/KnownSppMSOM.AltSize.JAGS.R"
models[[14]] <- "Source_Files/JAGS_models/KnownSppMSOM.AltStream.JAGS.R"
models[[15]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSize.JAGS.R"
models[[16]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSize2.JAGS.R"
models[[17]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleAltStream.JAGS.R"
models[[18]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleAltStream2.JAGS.R"
models[[19]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleLandCoverSize.JAGS.R"
models[[20]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleLandCoverSize2.JAGS.R"

models[[21]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSizeCount.JAGS.R"
models[[22]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSizeCountNoYear.JAGS.R"
models[[23]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSize2.JAGS.R"

models[[24]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSizeCountNoYear.JAGS.R"
models[[25]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSize2noYr.JAGS.R"
models[[26]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSizeNoYr.JAGS.R"
models[[27]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleAltStream2noYr.JAGS.R"
models[[28]] <- "Source_Files/JAGS_models/KnownSppMSOM.SimpleSizeNoYrConnect.JAGS.R"


### Fit Models!  ---------------------------------

## Create probability expressions
tmp <- createExpressionsArk()
expressPsi <- tmp$expressPsi
expressDet <- tmp$expressDet
## save output in lists for organization:
out <- list()
score01 <- vector(length=nModels)

for (i in 24:nModels){
  nPsiCoef <- nPsiList[[i]]  # elev + intecept
  nDetCoef <- nDetList[[i]] # year, elev
  set.seed(5000)
  out[[i]] <- do.call(jags.parallel, list(jags.data, jags.inits, params, 
                                          models[[i]],
                                          nc, ni, nb, nt));
  #     out[[i]] <- jags(jags.data, jags.inits, params, 
  #                      models[[i]], nc, ni, nb, nt);
  cat("Mod ", i, "done!\n")
  #save.image("~/testRenvProject/Output_Files/4_RunModels/ArkInSampleModels.RData") # save often in case it crashes
}



### ------------------
###  Calc in-sample model selection criterion:
###           Deviance, WAIC, CPO 
### ------------------

# First, gets MCMC samples:
#  --> Run all models
# base::load("~/Dropbox/Fish_PostDoc/ArkansasR/ArkModels/ArkInSampleModels.RData")
#base::load("~/testRenvProject/Output_Files/4_RunModels/ArkInSampleModels.RData")
# setwd("/Volumes/CPW_Work/Optimum Sampling/Ark_Optimal_bwa/Output_Files/")

## Calculate and output model selection for all five models
waic.out <- list()
Probs <- list()

for(i in 24:nModels){
  Probs[[i]] <- derivePsiDetect(out[[i]], 
                                expressPsi[[i]], expressDet[[i]],
                                nPossSpp, nSites, nSurveys)
  lik <- calcLik(Probs[[i]], Y, nMCMC=out[[i]]$BUGSoutput$n.sims, 
                 nSites=nSites, nSpp=nPossSpp)
  waic.out[[i]] <- waic.cpo(lik, nMCMC=out[[i]]$BUGSoutput$n.sims)
  score01[i] <- mean( apply(out[[i]]$BUGSoutput$sims.list$score01ij, 
                            c(2, 3), sum, na.rm=T) )
  print(i)
}

# and check model selection stats:
tmp <- stack(unlist(waic.out))#, FALSE))
# tmp
tmp$time <- rep(24:nModels, each=6)
modSelect.tab <- round(as.matrix(reshape(tmp, idvar="time", timevar="ind", 
                                         v.names="values", direction="wide")))
inSample.stats <- data.frame(Model=unlist(models)[24:28], modSelect.tab, score01[24:nModels])
inSample.stats$nDF <- nPsiList[24:nModels] + nDetList[24:nModels]
# inSample.stats
# inSample.stats[order(inSample.stats$nDF), ]
write.csv(inSample.stats, "~/testRenvProject/Output_Files/4_RunModels/InSampleStats.csv", row.names=F)
write.csv(mtcars, "~/testRenvProject/Output_Files/4_RunModels/InSampleStats.csv", row.names=F)
# end of file. ---------------------------------------------------------

