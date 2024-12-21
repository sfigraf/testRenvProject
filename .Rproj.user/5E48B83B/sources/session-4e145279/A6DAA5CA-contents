

###   This script imports data, the models, and the spatially-balanced sites
###  And outputs sites selected through OASD.
###  10 sites per year, for 3 years: 2016, 2017, 2018

## INPUTS:
## ArkData.RData
## ArkModels.RData
## ArkBalancedSites.RData


##OUTPUTS:
##ArkOASDsites.RData


## Set directory:
setwd("~/CPWOptimalSampling/Output_Files/8_SelectOASDsites")  # Fill in as appropriate
# setwd("~/ArkSelectSites")  # for server

## LIBRARIES:

# the data and source functions -----------------------
load("~/CPWOptimalSampling/Output_Files/2_OrganizeNewData/ArkData.Rdata")
load("~/CPWOptimalSampling/Output_Files/4_RunModels/ArkInSampleModels.RData")
load("~/CPWOptimalSampling/Output_Files/7_SelectBalancedSites/ArkBalancedSites.RData") # load("ArkBalancedSites.RData")
source("~/CPWOptimalSampling/Source_Files/ArkFunctions.R")
# load("ArkData.Rdata")
# load("ArkInSampleModels.RData")
# load("ArkBalancedSites.RData")
# source("ArkFunctions.R")

mod.no <- 25
mod <- out[[mod.no]]
rm(out)  # keep more workspace memory available
nPsiCoef <- nPsiList[[mod.no]]
nDetCoef <- nDetList[[mod.no]]

### Organize sites and covariates --------------------------
names(sites)  # sites that were already sampled
names(pred_sites)  # potential sites to be surveyed in future
max(sites$SiteID)
pred_sites$SiteID <- 1000 + 1:nrow(pred_sites) # all sites need an ID
## Add some columns to pred_sites to match orig sites data frame:
pred_sites$nSurveys <- 0  # potential sites have not been surveyed yet
pred_sites$aboveP <- ifelse(pred_sites$RESERVOIR == "aboveP", 1, 0)
pred_sites$belowJM <- ifelse(pred_sites$RESERVOIR == "belowJM", 1, 0)
pred_sites$X2 <- pred_sites$X ^ 2
pred_sites$Y2 <- pred_sites$Y ^2
AllSites <- rbind(sites, pred_sites[, c("SiteID", "UTMX", "UTMY", "X", "Y",
                                        "reachid", "pointid", "UNCONNECT",
                                        "GRAD", "ELEV", "ELEV2", "SIZE",
                                        "CROPS", "DVLPD", "WTLNDS", "RESERVOIR",
                                        "FTN", "MAIN", "PURG", "nSurveys",
                                        "aboveP", "belowJM", "X2", "Y2")])
nAllSites <- nrow(AllSites)
AllSites$OASDid <- 1:nAllSites

AllpsiX.mat <- cbind(1, AllSites[, c("SIZE")])
## Site-level covariate(s) needed for detection prob's:
SIZE <- AllSites[, "SIZE"] # c(sites[, c("SIZE")], pred_sites[, c("SIZE")])

# Check detection expression to build orig_detX.array:
expressDet[[mod.no]]

# Build array of covariates for original sites sampled:
orig_detX.array <- array(NA, dim=c(nDetCoef, nSites, maxSurveys))
orig_detX.array[1, , ] <- 1  # the intercept
orig_detX.array[2:(nDetCoef-1), , ] <- X.arrayAll[3, , ] #survey covariates
orig_detX.array[nDetCoef, , ] <- site.vars[11, ] # site covariate=StreamSize
# check that site-level cov filled in correctly:
# orig_detX.array[3, 3, ] 

# Covariate proxy for future sampling:
(meanCt <- mean(surveys$logTotalCt))

## Incorporate the balanced sites as already sampled:
## (but exclude the oversample)
extraSitesIndex <- which(AllSites$pointid %in% BalancedSites@data$pointid)[1:30]

### Exchange Algorithm Set Up -------------------------------------------------

# set up that might change:
nYears <- 3  # no. of years to pick OASD sites. 3 = 2016, 2017, and 2018
nRandom <- 1000  # how many random searches to begin algorithm with

nNbors <- 20 # how many n'bors to look at in the algorithm.
nAdd <- 10  # how many OASD sites do we want to sample each year?

#number of times to repeat the algorithm to get closer to optimal
( nStarts <- detectCores() / 2 )  # 8/2 = 4
nReps <- 2
# actual no. of times to repeat alg is nStarts * nReps
# nStarts = how many in parallel
# nReps = how many times sequntial


## Representative Species to include when calculating q
# repSpp <- which(natives %in% c("ARD", "LND", "SMM"))  # c(natives, "RCS") #
repSpp <- which(natives %in% c("BBH", "FHC", "SNF", "LND", "OSF",
                               "PKF", "RDS", "SAH", "SMM"))  # c(natives, "RCS") #
repSppIndices <- repSpp # 1:length(repSpp)

(nAllSites <- nrow(AllSites))
nOrigSites <- nrow(sites)  # == nSites
#
nMCMC <- mod$BUGSoutput$n.sims

# create matrix of nbors of each site:
allDist <- with(AllSites, iDist(UTMX, UTMY))
nbors <- matrix(NA, nrow=nAllSites, ncol=(nNbors + 1))
nbors[, 1] <- AllSites$OASDid
for(i in 1:nAllSites){
  nbors[i, -1] <- AllSites$OASDid[which(rank(allDist[i, ], 
                                             ties.method="random") <= nNbors )]
}

### Loop Set Up ------------------------------------

AllSurveysDetectProbs <- createDetectProbs(repSppIndices, mod, 
                                           nSurveys, nAllSites, nOrigSites, 
                                           nDetCoef, orig_detX.array, SIZE, meanCt)
AllSitesOccuProbs <- createOccuProbs(repSppIndices, mod, AllpsiX.mat, 
                                     nAllSites, nMCMC)

# # Add the OASD sites from 2016:
# sites2016Index <- qList2016[[3]]$allIndices[2, ]
# extraSitesIndex <- c(extraSitesIndex, sites2016Index)
# sites2017Index <- qList2017[[2]]$allIndices[2, ]
# extraSitesIndex <- c(extraSitesIndex, sites2017Index)

### Calc q for no/all/random sites sampled --------------------------

# no additional sites sampled:
high.q <- sum(apply(calc_q(repSppIndices, mod, 
                           AllSitesOccuProbs, AllSurveysDetectProbs,
                           nAllSites, extraSitesIndex,
                           nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# high.q  # 778.9192
# q per site, per species:  
# (high.q / length(repSppIndices) ) / nAllSites  # 0.13596

low.q <- sum(apply(calc_q(repSppIndices, mod, 
                          AllSitesOccuProbs, AllSurveysDetectProbs,
                          nAllSites, 1:nAllSites,
                          nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# low.q  # 237.8205
# (low.q / length(repSppIndices) ) / nAllSites  # 0.043
#

# and for a random design search
set.seed(1112)
RandomSites <- sample(1:nAllSites, nAdd, replace=T)
random.q <- sum(apply(calc_q(repSppIndices, mod, 
                             AllSitesOccuProbs, AllSurveysDetectProbs,
                             nAllSites, 
                             c(extraSitesIndex, RandomSites),
                             nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# random.q  # 769.866
# (random.q / length(repSppIndices) ) / nAllSites  # 0.138
#

### Exchange Algorithm Search -------------------------------------

OASDsites <- array(NA, dim=c(nYears, nAdd))
OASDq <- OASD.qInd <- vector(length=nYears)
qAll <- list()
SitesAll <- list()

for (yr in 1:nYears) {
  qAll[[yr]] <- vector(length=nReps*nStarts)
  SitesAll[[yr]] <- array(NA, dim=c(nReps*nStarts, nAdd))
  for(j in 1:nReps){
    
    #setup parallel backend to use 8 processors
    cl <- makeCluster(nStarts)
    registerDoParallel(cl)
    #qList <- list()
    #options(error=recover)
    set.seed(216 + j*10 + yr*100)
    qListTmp <- foreach(s=1:nStarts, .packages='boot') %dorng%  {
      
      # Pick random set of sites to start with:
      tmpIndex <- list()
      tmp.q <- vector(length=nRandom)
      for(i in 1:nRandom){
        # set.seed(s * 100 + i)
        tmpIndex[[i]] <- sample(1:nAllSites, nAdd, replace=F)
        
        ## Calculate q for starting set of sites:
        tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                                     AllSitesOccuProbs, AllSurveysDetectProbs,
                                     nAllSites, 
                                     c(extraSitesIndex, unlist(tmpIndex[[i]])),
                                     nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
      }
      new_q <- min(tmp.q)
      qIndex <- which(tmp.q == new_q)[1]
      allIndices <- tmpIndex[[qIndex]]
      all.q <- new_q
      
      prev_q <- 0  # to start the while loop
      previousIndex <- tmpIndex[[qIndex]]
      nIter <- 0  # keep track of how many it loops through the sites
      print(nIter)
      while (prev_q != new_q) {
        prev_q <- new_q
        nIter <- nIter + 1
        
        ## Recalculate q for each NN for each current design point
        for (k in 1:nAdd){
          optimalIndex <- previousIndex
          nbor.q <- vector(length=dim(nbors)[2])
          nbor.q[1] <- new_q
          for (j in 1:nNbors){
            # replace optimal site with one of its n'bors
            optimalIndex[k] <- nbors[previousIndex[k], j+1]  
            #         cat("Optimal Index", optimalIndex, "\n")
            ## Calculate q with sampling for new sites
            new.q <- calc_q(repSppIndices, mod, AllSitesOccuProbs, AllSurveysDetectProbs,
                            nAllSites, c(extraSitesIndex, optimalIndex),
                            nSurveys, nMCMC=nMCMC)
            nbor.q[j+1] <- sum( apply(new.q, 1, sum) )
            #         cat("Nbor.q", nbor.q, "\n")
            print(j)
          }
          new_q <- min(nbor.q)
          
            gc()
          
          previousIndex[k] <- nbors[previousIndex[k], which(nbor.q == new_q)[1]]
        }
        
        ## Save the indices
        allIndices <- rbind(allIndices, previousIndex)
        all.q <- c(all.q, new_q)
      }
      
      list(all.q=all.q, allIndices=allIndices, 
           min.q=new_q, bestSites = previousIndex)
    }
    stopCluster(cl)
    
    qTmp <- list()
    for (s in 1:nStarts){
      qTmp[[s]] <- qListTmp[[s]]$min.q
      SitesAll[[yr]][(j*nStarts - nStarts + 1):(j*nStarts), ][s, ] <- 
        unlist(qListTmp[[s]]$bestSites)
    }
    qAll[[yr]][(j*nStarts - nStarts + 1):(j*nStarts)] <- unlist(qTmp)
    cat("Rep ", j, "completed \n")
  }
  
  min.qInd <- which(qAll[[yr]] == min(qAll[[yr]]))[1]
  OASD.qInd[yr] <- min.qInd
  OASDq[yr] <- qAll[[yr]][min.qInd]
  OASDsites[yr, ] <- SitesAll[[yr]][min.qInd, ]
  
  extraSitesIndex <- c(extraSitesIndex, OASDsites[yr, ])
  cat("Year ", yr, "completed \n")  
  
  save.image("~/ArkAllSiteSelections.RData")
}

# save.image("~/ArkAllSiteSelections.RData")

## plot and export chosen sites --------------------------------
newsites2016 <- AllSites[SitesAll[[1]][5, ], ]
newsites2016b <- AllSites[SitesAll[[1]][6, ], ]
newsites2017 <- AllSites[OASDsites[2, ], ]
newsites2018 <- AllSites[OASDsites[3, ], ]
png("Ark_AllFutureSamplingLocations.png",  
    width=1024, height=768)
par(mfrow=c(1,1))
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n"))
# add locations already sampled:
points(sites$UTMX, sites$UTMY, pch=16, col=grey(0.4))
## And add in new OASD sites:
points(newsites2016$UTMX, newsites2016$UTMY, pch=8, cex=1.5, col=1)
points(newsites2017$UTMX, newsites2017$UTMY, pch=8, cex=1.5, col="red")
points(newsites2018$UTMX, newsites2018$UTMY, pch=8, cex=1.5, col="cyan")
## And add GRTS sites:
with(subset(BalancedSites@data, panel=="Year1"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
with(subset(BalancedSites@data, panel=="Year2"), 
     points(xcoord, ycoord, col="red", pch=16, cex=1.5))
with(subset(BalancedSites@data, panel=="Year3"), 
     points(xcoord, ycoord,  col="cyan", pch=16, cex=1.5))
dev.off()

#(more plotting and output exporting in "ArkPlotFutureSites.R" file)


# tmp1 <- subset(AllSites, SIZE==1)
# site1 <- sample(tmp1$OASDid, 1)
# tmp2 <- subset(AllSites, SIZE==2)
# sites2 <- sample(tmp2$OASDid, 9)
# tmp.index <- which(AllSites %in% c(site1, site2))
# test.q <- sum(apply(calc_q(repSppIndices, mod, 
#                              AllSitesOccuProbs, AllSurveysDetectProbs,
#                              nAllSites, 
#                              c(extraSitesIndex, site1, sites2),
#                              nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# test.q  # 769.866
# (test.q / length(repSppIndices) ) / nAllSites  # 0.138
# 

# end of file. ----------------------------------