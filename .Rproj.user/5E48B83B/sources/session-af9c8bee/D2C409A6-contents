
### OASD Functions -------------------------------
# Returns a 4-D array of detection probabilities:
#  for each "representative" species, each survey of each site, and each MCMC iteration
# Adds 9 new surveys to each site. (3 new surveys * max 3 samp occasions)
# Detect probs are then included (or not)
#  via the design criterion calculation function.
createDetectProbs <- function(repSppIndices, mod, 
                              nSurveys, nAllSites, nOrigSites,
                              nDetCoef,
                              orig_detX.array, 
                              SIZE, meanCt){
  # nSurveys = original vector of nSurveys.
  # nAllSites = 141 OrigSites + 479 potential new Pred sites
  # nOrigSites = nSites = 141 Orig Sites already sampled
  # nDetCoef = how many coeff associated with the detect probs?
  # orig_detX.array = original array of survey/detection covariates
  # SIZE = site-level covariate.  == AllSites$SIZE
  # meanCt = mean of log of TotalCt per occasion.  Use mean for future surveys.
  
  # Set up the nSurveys (changes if old sites are re-sampled):
  origMax_nSurveys <- max(nSurveys)
  nPredSites <- nAllSites - nOrigSites # new sites that will be sampled.
  new.nSurveys <- c(nSurveys, rep(0, nPredSites))
  potential.nSurveys <- new.nSurveys + 9 # assume 3 new surveys per site * max 3 visits
  (newMax_nSurveys <- max(potential.nSurveys))
  
  ## Create X_detect array for new surveys:
  detectX.array <- array(NA, dim=c(nDetCoef, nAllSites, newMax_nSurveys)) 
  detectX.array[, 1:nOrigSites, 1:origMax_nSurveys] <- orig_detX.array  # max of 24 surveys originally
  for (i in 1:nAllSites) {
    detectX.array[1, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- 1  # intercept
    #    detectX.array[2:3, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- 0  # YDAY, YDAY2
    
    detectX.array[2, i, (new.nSurveys[i] +1):(new.nSurveys[i] +2)] <- 1  # 2 electro passes
    detectX.array[2, i, (new.nSurveys[i] +3)] <- 0  # then 1 seine pass
    detectX.array[2, i, (new.nSurveys[i] +4):(new.nSurveys[i] +5)] <- 1  # 2 electro passes
    detectX.array[2, i, (new.nSurveys[i] +6)] <- 0  # 1 seine pass
    detectX.array[2, i, (new.nSurveys[i] +7):(new.nSurveys[i] +8)] <- 1  # 2 electro passes
    detectX.array[2, i, (new.nSurveys[i] +9)] <- 0  # 1 seine pass
    #    detectX.array[5, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- c(1, 2, 1) # Add PassNo
    
    ## assume detections all from year=2012 b/c that had lowest variability:
    #     detectX.array[6:8, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- 0  # 2009, 2010, 2011
    #     detectX.array[9, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- 1  # 2012
    #     detectX.array[10:12, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- 0  # 2013, 2014, 2015
    detectX.array[3, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- meanCt # assume mean count
#    detectX.array[4, i, (new.nSurveys[i] +1):potential.nSurveys[i]] <- SIZE[i]
  }
  
  ## Now create detect prob's using model output:
  AllSurveysDetectProbs <- array(NA, dim=c(length(repSppIndices), 
                                           nAllSites, newMax_nSurveys, nMCMC))
  
  # Loop over the species:
  for (i in seq_along(repSppIndices)) {  
    # Loop over the sites:
    for (j in 1:nAllSites) {  
      # Loop over the surveys:
      for (k in 1:potential.nSurveys[j]) {  
        AllSurveysDetectProbs[ i, j, k, ] <- 
          inv.logit( t(as.matrix(detectX.array[, j, k])) %*% 
                       t(mod$BUGSoutput$sims.list$detCoef[, repSppIndices[i], ]))
      }
    }
    print(i)
  }
  return(AllSurveysDetectProbs)
}

createOccuProbs <-function(repSppIndices, mod, psiX.mat, nAllSites, nMCMC){
  # Returns a 3-D array of occupancy probabilities:
  #  for each "representative" species, each site, and each MCMC iteration
  
  # repSppIndices = indices of species to follow
  # mod = JAGS output
  # psiX.mat= Occupancy covariate values at all possible sites
  # nAllSites = 141 previously sampled sites + 479 new sites
  # nMCMC = 3000 samples
  
  occuMat <- array(NA, dim=c(length(repSppIndices), nAllSites, nMCMC))
  for (i in seq_along(repSppIndices)) {  
    occuMat[i, , ] <- inv.logit( as.matrix(psiX.mat) %*% 
                                   t(mod$BUGSoutput$sims.list$psiCoef[, repSppIndices[i], ]) )
  }
  return(occuMat)
}


calc_q <- function (repSppIndices, mod, AllSitesOccuProbs, AllSurveysDetectProbs, 
                    nAllSites, newSiteIndex,
                    nSurveys, nMCMC) {
  # Function calculates the design criterion for all representative species.
  # the returned q has dim = nRepSpecies X nAllSites
  
  # repSppIndices = indices of species to follow
  # mod = JAGS output
  # AllSitesOccuProbs = Array: rep. species X nAllSites X nMCMC matrix of occupancy prob's
  # AllSurveysDetectProbs = Detect prob. values at all possible sites,all possible surveys
  # nAllSites = 113 previously sampled sites + ~585 new sites
  # newSiteIndex = indices of sites associated with current design, to be sampled in future
  # nSurveys = vector of original no. of surveys per site, for resampled sites
  
  # Set up the new nSurveys (changes if old sites are re-sampled):
  # nPredSites = nAllSites - previously sampled sites
  nOrigSites <- length(nSurveys)
  nPredSites <- nAllSites - nOrigSites
  design.nSurveys <- c(nSurveys, rep(0, nPredSites))
  for(m in newSiteIndex){
    design.nSurveys[m] <- design.nSurveys[m] + 3 # assume 3 new surveys per site
  }
  #  (newMax_nSurveys <- max(design.nSurveys))
  
  # 9 species of interest to follow:
  cum.det <- array(1, dim=c(length(repSppIndices), nAllSites, nMCMC))
  
  # Loop over the species:
  for (i in seq_along(repSppIndices)) {  
    # Loop over the sites:
    for (m in 1:nAllSites) {  
      # Loop over the surveys, if the site was surveyed:
      if (design.nSurveys[m] > 0){
        for (j in seq_along(design.nSurveys[m])) {  
          cum.det[i, m, ] <- (1 - AllSurveysDetectProbs[i, m, j, ]) * cum.det[i, m, ]
        }
      }
    }
  }
  
  psi.tilde <- AllSitesOccuProbs * cum.det / 
    (AllSitesOccuProbs * cum.det + 1 - AllSitesOccuProbs) 
  
  # take mean over the MCMC samples:
  q <- apply(psi.tilde * (1 - psi.tilde), c(1, 2), mean, na.rm=T)  
  
  # And set Var(q) = 0 for any site where species was detected.
  q[, 1:nOrigSites] <- ifelse(Ybinom[repSppIndices, ] > 0, 0, q[, 1:nOrigSites])
  
  return(q)
}


### Function to run MSOMs ---------------

# Function to calc. mode of posterior from mcmc sample
mode.est <- function(z){
  if(var(z)==0) return(z[1])  # code-thing from Devin.  To make sure no error
  bp <- boxplot(z, plot=FALSE) # calc to remove extreme values from density fxn
  dens <- density(z, from=bp$stats[1,1], to=bp$stats[5,1])
  return(dens$x[dens$y==max(dens$y)])
}

### Probability expressions  ------------------

## Input the logit-function expressions used in JAGS 
## (because data is 2+ dimensions, cannot use matrix multiplication easily)

# The expressions are required to derive quantities 
#  that were not saved with the MCMC
# mod.type must equal "AllData" or "CrossVal"
createExpressionsArk <- function(mod.type="AllData"){
  expressPsi <- list()
  expressDet <- list()
    expressPsi[[1]] <- expression(inv.logit(psiCoef[, i, 1]  + psiCoef[, i, 2] * site.vars[1, j] +
                                              psiCoef[, i, 3] * site.vars[2, j] + psiCoef[, i, 4]*site.vars[3, j] +
                                              psiCoef[, i, 5]*site.vars[4, j] + psiCoef[, i, 6]*site.vars[5, j] +
                                              psiCoef[, i, 7]*site.vars[6, j] + psiCoef[, i, 8]*site.vars[7, j] +
                                              psiCoef[, i, 9]*site.vars[8, j] + psiCoef[, i, 10]*site.vars[9, j] +
                                              psiCoef[, i, 11]*site.vars[10, j] + psiCoef[, i, 12]*site.vars[11, j] +
                                              psiCoef[, i, 13]*site.vars[12, j] + psiCoef[, i, 14]*site.vars[13, j] +
                                              psiCoef[, i, 15]*site.vars[14, j] + psiCoef[, i, 16]*site.vars[15, j] +
                                              psiCoef[, i, 17]*site.vars[16, j]))
    expressPsi[[2]] <- expression(inv.logit( psiCoef[, i, 1]  + psiCoef[, i, 2] * site.vars[1, j] +
                                               psiCoef[, i, 3] * site.vars[2, j] + psiCoef[, i, 4]*site.vars[3, j] +
                                               psiCoef[, i, 5]*site.vars[4, j] + psiCoef[, i, 6]*site.vars[5, j] +
                                               psiCoef[, i, 7]*site.vars[6, j] + psiCoef[, i, 8]*site.vars[7, j] +
                                               psiCoef[, i, 9]*site.vars[8, j] + psiCoef[, i, 10]*site.vars[9, j] +
                                               psiCoef[, i, 11]*site.vars[10, j] + psiCoef[, i, 12]*site.vars[11, j] +
                                               psiCoef[, i, 13]*site.vars[12, j] + psiCoef[, i, 14]*site.vars[13, j] +
                                               psiCoef[, i, 15]*site.vars[14, j] + psiCoef[, i, 16]*site.vars[15, j] +
                                               psiCoef[, i, 17]*site.vars[16, j]))
    expressPsi[[3]] <- expression(inv.logit(
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
        psiCoef[, i, 5]*site.vars[6, j] + psiCoef[, i, 6]*site.vars[7, j] +
        psiCoef[, i, 7]*site.vars[8, j] + psiCoef[, i, 8]*site.vars[9, j] +
        psiCoef[, i, 9]*site.vars[10, j] + psiCoef[, i, 10]*site.vars[11, j] +
        psiCoef[, i, 11]*site.vars[12, j] + psiCoef[, i, 12]*site.vars[13, j] +
        psiCoef[, i, 13]*site.vars[14, j] + psiCoef[, i, 14]*site.vars[15, j] +
        psiCoef[, i, 15]*site.vars[16, j] ))
    expressPsi[[4]] <- expression(inv.logit(
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
        psiCoef[, i, 5]*site.vars[10, j] +
        psiCoef[, i, 6]*site.vars[12, j] + psiCoef[, i, 7]*site.vars[13, j] +
        psiCoef[, i, 8]*site.vars[14, j] + psiCoef[, i, 9]*site.vars[15, j] ))
    expressPsi[[5]] <- expression(inv.logit(
        psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + 
        psiCoef[, i, 4]*site.vars[7, j] +
        psiCoef[, i, 5]*site.vars[8, j] + psiCoef[, i, 6]*site.vars[9, j] +
        psiCoef[, i, 7]*site.vars[11, j] ))
    expressPsi[[6]] <- expression(inv.logit(
        psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
        psiCoef[, i, 5]*site.vars[8, j] + psiCoef[, i, 6]*site.vars[9, j] +
        psiCoef[, i, 7]*site.vars[10, j] + psiCoef[, i, 8]*site.vars[12, j] +
        psiCoef[, i, 9]*site.vars[13, j]+ psiCoef[, i, 10]*site.vars[14, j] ))
    expressPsi[[7]] <- expression(inv.logit(
        psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
        psiCoef[, i, 5]*site.vars[6, j] + psiCoef[, i, 6]*site.vars[10, j] + 
        psiCoef[, i, 7]*site.vars[12, j] + psiCoef[, i, 8]*site.vars[13, j]  )) 
    expressPsi[[8]] <- expression(inv.logit(
        psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
        psiCoef[, i, 5]*site.vars[6, j] + psiCoef[, i, 6]*site.vars[10, j] + 
        psiCoef[, i, 7]*site.vars[14, j] + psiCoef[, i, 8]*site.vars[15, j] +
        psiCoef[, i, 9]*site.vars[16, j] ))     
    expressPsi[[9]] <- expression(inv.logit(    
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
      psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
      psiCoef[, i, 5]*site.vars[6, j] + psiCoef[, i, 6]*site.vars[11, j] ))
    expressPsi[[10]] <- expression(inv.logit(
        psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[4, j] + psiCoef[, i, 4]*site.vars[5, j] +
        psiCoef[, i, 5]*site.vars[10, j] + 
        psiCoef[, i, 6]*site.vars[12, j] + psiCoef[, i, 7]*site.vars[13, j] +
        psiCoef[, i, 8]*site.vars[14, j] + psiCoef[, i, 9]*site.vars[16, j] +
        psiCoef[, i, 10]*site.vars[11, j] )) 
    expressPsi[[11]] <- expression(inv.logit(
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
        psiCoef[, i, 3]*site.vars[5, j] +
        psiCoef[, i, 4]*site.vars[9, j] + psiCoef[, i, 5]*site.vars[10, j] +
        psiCoef[, i, 6]*site.vars[12, j] + psiCoef[, i, 7]*site.vars[13, j] +
        psiCoef[, i, 8]*site.vars[14, j]+ psiCoef[, i, 9]*site.vars[15, j] +
        psiCoef[, i, 10]*site.vars[16, j]  )) 

    expressPsi[[12]] <- expression(inv.logit(    
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[3, j] +
      psiCoef[, i, 3]*site.vars[4, j] + 
      psiCoef[, i, 4]*site.vars[7, j] +
      psiCoef[, i, 5]*site.vars[8, j] + psiCoef[, i, 6]*site.vars[9, j] +
      psiCoef[, i, 7]*site.vars[11, j]  ))
    expressPsi[[13]] <- expression(inv.logit(  
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[1, j] +
        psiCoef[, i, 3]*site.vars[2, j] + psiCoef[, i, 4]*site.vars[3, j] +
        psiCoef[, i, 5]*site.vars[4, j] +
        psiCoef[, i, 6]*site.vars[11, j] ))
    expressPsi[[14]] <- expression(inv.logit(   
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[1, j] +
        psiCoef[, i, 3]*site.vars[2, j] + psiCoef[, i, 4]*site.vars[3, j] +
        psiCoef[, i, 5]*site.vars[4, j] + psiCoef[, i, 6]*site.vars[10, j] +
        psiCoef[, i, 7]*site.vars[14, j] + psiCoef[, i, 8]*site.vars[15, j] +
        psiCoef[, i, 9]*site.vars[16, j]))
    expressPsi[[15]] <- expression(inv.logit(    
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[11, j] ))
    expressPsi[[16]] <- expression(inv.logit( 
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[11, j] ))
    expressPsi[[17]] <- expression(inv.logit(  
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[10, j] +
        psiCoef[, i, 3]*site.vars[14, j] + psiCoef[, i, 4]*site.vars[15, j] +
        psiCoef[, i, 5]*site.vars[16, j]))
    expressPsi[[18]] <- expression(inv.logit( 
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[10, j] +
        psiCoef[, i, 3]*site.vars[14, j] + psiCoef[, i, 4]*site.vars[15, j] +
        psiCoef[, i, 5]*site.vars[16, j]))
    expressPsi[[19]] <- expression(inv.logit(   
      psiCoef[, i, 1]  + 
        psiCoef[, i, 2]*site.vars[7, j] +
        psiCoef[, i, 3]*site.vars[8, j] + psiCoef[, i, 4]*site.vars[9, j] +
        psiCoef[, i, 5]*site.vars[11, j]  ))
    expressPsi[[20]] <- expression(inv.logit(    
      psiCoef[, i, 1]  + 
        psiCoef[, i, 2]*site.vars[7, j] +
        psiCoef[, i, 3]*site.vars[8, j] + psiCoef[, i, 4]*site.vars[9, j] +
        psiCoef[, i, 5]*site.vars[11, j]  )) 
    
    expressPsi[[21]] <- expressPsi[[16]]
    expressPsi[[22]] <- expressPsi[[16]]
    expressPsi[[23]] <- expressPsi[[16]]
                      
    expressPsi[[24]] <- expressPsi[[22]]  
    expressPsi[[25]] <- expressPsi[[22]]
    expressPsi[[26]] <- expressPsi[[22]]
    expressPsi[[27]] <- expression(inv.logit(    
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[10, j] +
        psiCoef[, i, 3]*site.vars[14, j] + psiCoef[, i, 4]*site.vars[15, j] +
        psiCoef[, i, 5]*site.vars[16, j] ))
    expressPsi[[28]] <- expression(inv.logit(  
      psiCoef[, i, 1]  + psiCoef[, i, 2]*site.vars[11, j] +
        psiCoef[, i, 3]*site.vars[10, j] ))

    
    expressDet[[1]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[3, j, k] +
        detCoef[, i, 5]*X.array[4, j, k] + 
        detCoef[, i, 6]*X.array[6, j, k] + detCoef[, i, 7]*X.array[7, j, k] +
        detCoef[, i, 8]*X.array[8, j, k]  + detCoef[, i, 9]*X.array[9, j, k]  +
        detCoef[, i, 10]*X.array[10, j, k]  + detCoef[, i, 11]*X.array[11, j, k]  +
        detCoef[, i, 12]*X.array[12, j, k] +
        detCoef[, i, 13]*site.vars[1, j] + detCoef[, i, 14]*site.vars[2, j] +
        detCoef[, i, 15]*site.vars[3, j] + detCoef[, i, 16]*site.vars[4, j]  +
        detCoef[, i, 17]*site.vars[5, j]  + detCoef[, i, 18]*site.vars[6, j] +
        detCoef[, i, 19]*site.vars[7, j]  + detCoef[, i, 20]*site.vars[8, j] +
        detCoef[, i, 21]*site.vars[9, j]  + detCoef[, i, 22]*site.vars[10, j] +
        detCoef[, i, 23]*site.vars[11, j]  + detCoef[, i, 24]*site.vars[12, j] +
        detCoef[, i, 25]*site.vars[13, j] + detCoef[, i, 26]*site.vars[14, j] +
        detCoef[, i, 27]*site.vars[15, j] + detCoef[, i, 28]*site.vars[16, j]))
    expressDet[[2]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[3, j, k] +
        detCoef[, i, 5]*X.array[4, j, k] + 
        detCoef[, i, 6]*X.array[6, j, k] + detCoef[, i, 7]*X.array[7, j, k] +
        detCoef[, i, 8]*X.array[8, j, k]  + detCoef[, i, 9]*X.array[9, j, k]  +
        detCoef[, i, 10]*X.array[10, j, k]  + detCoef[, i, 11]*X.array[11, j, k]  +
        detCoef[, i, 12]*X.array[12, j, k]  ))
    
    expressDet[[3]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[3, j, k] +
        detCoef[, i, 5]*X.array[6, j, k] + detCoef[, i, 6]*X.array[7, j, k] +
        detCoef[, i, 7]*X.array[8, j, k]  + detCoef[, i, 8]*X.array[9, j, k]  +
        detCoef[, i, 9]*X.array[10, j, k]  + detCoef[, i, 10]*X.array[11, j, k]  +
        detCoef[, i, 11]*X.array[12, j, k]  + detCoef[, i, 12]*site.vars[11, j]))
    expressDet[[4]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[3, j, k] +
        detCoef[, i, 5]*X.array[4, j, k] +
        detCoef[, i, 6]*X.array[6, j, k] + detCoef[, i, 7]*X.array[7, j, k] +
        detCoef[, i, 8]*X.array[8, j, k]  + detCoef[, i, 9]*X.array[9, j, k]  +
        detCoef[, i, 10]*X.array[10, j, k]  + detCoef[, i, 11]*X.array[11, j, k]  +
        detCoef[, i, 12]*X.array[12, j, k] ))
    expressDet[[5]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[3, j, k] +
        detCoef[, i, 5]*X.array[4, j, k] +
        detCoef[, i, 6]*X.array[6, j, k] + detCoef[, i, 7]*X.array[7, j, k] +
        detCoef[, i, 8]*X.array[8, j, k]  + detCoef[, i, 9]*X.array[9, j, k]  +
        detCoef[, i, 10]*X.array[10, j, k]  + detCoef[, i, 11]*X.array[11, j, k]  +
        detCoef[, i, 12]*X.array[12, j, k]  + detCoef[, i, 13]*site.vars[11, j]))
    expressDet[[6]] <- expression(inv.logit(
        detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[4, j, k] +
        detCoef[, i, 5]*X.array[6, j, k] + detCoef[, i, 6]*X.array[7, j, k] +
        detCoef[, i, 7]*X.array[8, j, k]  + detCoef[, i, 8]*X.array[9, j, k]  +
        detCoef[, i, 9]*X.array[10, j, k]  + detCoef[, i, 10]*X.array[11, j, k]  +
        detCoef[, i, 11]*X.array[12, j, k]))
    
    expressDet[[7]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[3, j, k] +
        detCoef[, i, 4]*X.array[6, j, k] + detCoef[, i, 5]*X.array[7, j, k] +
        detCoef[, i, 6]*X.array[8, j, k]  + detCoef[, i, 7]*X.array[9, j, k]  +
        detCoef[, i, 8]*X.array[10, j, k]  + detCoef[, i, 9]*X.array[11, j, k]  +
        detCoef[, i, 10]*X.array[12, j, k]))
    expressDet[[8]] <- expressDet[[7]]
    expressDet[[9]] <- expressDet[[7]]
    expressDet[[10]] <- expressDet[[7]]
    expressDet[[11]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
      detCoef[, i, 3]*X.array[4, j, k] +
      detCoef[, i, 4]*X.array[6, j, k] + detCoef[, i, 5]*X.array[7, j, k] +
      detCoef[, i, 6]*X.array[8, j, k]  + detCoef[, i, 7]*X.array[9, j, k]  +
      detCoef[, i, 8]*X.array[10, j, k]  + detCoef[, i, 9]*X.array[11, j, k]  +
      detCoef[, i, 10]*X.array[12, j, k]  + detCoef[, i, 11]*site.vars[11, j]))
    
    expressDet[[12]] <- expressDet[[5]]
    expressDet[[13]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[2, j, k] + detCoef[, i, 4]*X.array[3, j, k] +
        detCoef[, i, 5]*X.array[6, j, k] + detCoef[, i, 6]*X.array[7, j, k] +
        detCoef[, i, 7]*X.array[8, j, k]  + detCoef[, i, 8]*X.array[9, j, k]  +
        detCoef[, i, 9]*X.array[10, j, k]  + detCoef[, i, 10]*X.array[11, j, k]  +
        detCoef[, i, 11]*X.array[12, j, k]  + detCoef[, i, 12]*site.vars[11, j] ))
    expressDet[[14]] <- expressDet[[13]]
    expressDet[[15]] <- expressDet[[13]]
    expressDet[[16]] <- expression(inv.logit(
          detCoef[, i, 1]  + detCoef[, i, 2]*X.array[3, j, k] +
            detCoef[, i, 3]*X.array[6, j, k] + detCoef[, i, 4]*X.array[7, j, k] +
            detCoef[, i, 5]*X.array[8, j, k]  + detCoef[, i, 6]*X.array[9, j, k]  +
            detCoef[, i, 7]*X.array[10, j, k]  + detCoef[, i, 8]*X.array[11, j, k]  +
            detCoef[, i, 9]*X.array[12, j, k]  + detCoef[, i, 10]*site.vars[11, j]))
    expressDet[[17]] <- expressDet[[15]]
    expressDet[[18]] <- expressDet[[16]]
    expressDet[[19]] <- expressDet[[12]]
    expressDet[[20]] <- expressDet[[18]]
    
    expressDet[[21]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[3, j, k] +
        detCoef[, i, 3]*X.array[6, j, k] + detCoef[, i, 4]*X.array[7, j, k] +
        detCoef[, i, 5]*X.array[8, j, k]  + detCoef[, i, 6]*X.array[9, j, k]  +
        detCoef[, i, 7]*X.array[10, j, k]  + detCoef[, i, 8]*X.array[11, j, k]  +
        detCoef[, i, 9]*X.array[12, j, k]  + detCoef[, i, 10]*X.array[13, j, k]  + 
        detCoef[, i, 11]*site.vars[11, j] ))
    expressDet[[22]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[3, j, k] +
        detCoef[, i, 3]*X.array[13, j, k]  + detCoef[, i, 4]*site.vars[11, j] ))
    expressDet[[23]] <- expressDet[[16]]
    
    expressDet[[24]] <- expressDet[[22]]
    expressDet[[25]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[3, j, k] +
        detCoef[, i, 3]*site.vars[11, j] ))
    expressDet[[26]] <- expression(inv.logit(
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[1, j, k] +
        detCoef[, i, 3]*X.array[3, j, k] +
        detCoef[, i, 4]*site.vars[11, j] ))
    expressDet[[27]] <- expression(inv.logit(    
      detCoef[, i, 1]  + detCoef[, i, 2]*X.array[3, j, k] +
        detCoef[, i, 3]*site.vars[11, j] ))
    expressDet[[28]] <- expressDet[[22]]

  return( list(expressPsi=expressPsi, expressDet=expressDet))
}


### Derived Psi's, Detect's from Jags output  ---------------
# Use the distributions associated with each coefficient to recreate
#  the posteriors associated with the coefficients for each species
derivePsiDetect <- function(jags.out, expressPsi, expressDet, 
                            nPossSpp, nSites, nSurveys){
  # jags.out = MCMC output from JAGS
  # expressPsi = logit function formula to calculate occupancy prob's.
  #              (equivalent to X *%* beta but must be written as vectors for JAGS)
  # expressDet = logit function formula to calculate detection prob's.
  #              (equivalent to X *%* beta but must be written as vectors for JAGS)
  # nPossSpp = 16 native species in the system (24 were detected)

  require(boot)  # for inv.logit
  (maxSurveys <- max(nSurveys))
  n.saved <- jags.out$BUGSoutput$n.sims
  psiCoef <- jags.out$BUGSoutput$sims.list$psiCoef  # n.saved X nPossSpp X nPsiCoef  
  psiS <- array(NA, dim=c(n.saved, nPossSpp, nSites) )
  detCoef <- jags.out$BUGSoutput$sims.list$detCoef  # n.saved X nPossSpp X nDetCoef
  detectS <- array(NA, dim=c(n.saved, nPossSpp, nSites, maxSurveys) )
  
  integrated.probs <- array(NA, dim=c(n.saved, nPossSpp, nSites, maxSurveys) )
  for (i in 1:nPossSpp) {
    for (j in 1:nSites) {
      psiS[ , i, j] <- eval(expressPsi)  
      for (k in 1:nSurveys[j]) {
        detectS[ , i, j, k] <-  eval(expressDet)
        integrated.probs[ , i, j, k] <- psiS[ , i, j] * detectS[ , i, j, k] 
      }
    }
  }
  return(list(psiS=psiS, detectS=detectS, integrated.probs=integrated.probs))
}

### Calculate likelihoods from JAGS output  ---------------
## Likelihood
dmix <- function(y, p, psi) {
  require(matrixStats)  # for rowProds() functions
  # out <- rep(0, length(y))
  #zero.idx <- (y == 0)
  z.tmp <- ifelse(rowSums(y, na.rm=T) > 0, 1, 0)
  out <- rep(0, dim(y)[1])  # dim(y)[1] == nSites
  zero.idx <- (z.tmp == 0)
  out[zero.idx] <- 1 - psi[zero.idx] + 
    psi[zero.idx] * rowProds(1-p[zero.idx, ], na.rm=T) #^ J[zero.idx]
  if( sum(!zero.idx) > 1) {
    out[!zero.idx] <- psi[!zero.idx] * 
      rowProds(dbinom(y[!zero.idx, ], 1, p[!zero.idx, ]), na.rm=T)
  } else{
    out[!zero.idx] <- psi[!zero.idx] * 
      prod(dbinom(y[!zero.idx, ], 1, p[!zero.idx, ]), na.rm=T)
  }
  out
}

# Likelihood value for all spp, sites, MCMC
calcLik <- function(probs.out, Y, nMCMC=3000, nSites=113, nSpp=16){
  lik <- array(NA, dim=c(nSpp, nMCMC, nSites))
  for(i in 1:nSpp){
    for(s in 1:nMCMC){
      lik[i, s, ] <- dmix(y=Y[i, , ], 
                          p=probs.out$detectS[s, i, , ], 
                          psi=probs.out$psiS[s, i, ])
    }
  }
  lik
}

### Calculate WAIC, CPO from JAGS output  ---------------
waic.cpo <- function(lik, nMCMC=3000){
  
  ## First calculate deviance:
  meanDev <- mean(-2 * apply( log(lik), 2, sum ))
  
  ## Then calc WAIC
  #  = the log of the mean, mean taken over the S MCMC iterations.
  ## and then sum over the nSpp, nSites
  (ellpd <- sum(log( apply(lik, c(1, 3), mean) ) ))
  
  #  (lppd <- -2 * sum( log(apply(lik, c(2, 3), mean)), na.rm=T ) )
  pd <- sum( apply( log(lik), c(1, 3), var, na.rm=T ) )
  pd2 <- 2 * sum( (log( apply(lik, c(1, 3), mean) ) -
                     apply( log(lik), c(1, 3), mean) ) )
  waic1 <- -2*ellpd + 2 * pd
  waic2 <- -2*ellpd + 2 * pd2
  
  
  # And calculate CPO
  cpo.i <- nMCMC / apply(1 / lik, c(1, 3), sum)
  cpo <- -1 * sum( log(cpo.i), na.rm=T)
  
  return(list(pD=pd, pD2=pd2, waic1=waic1, waic2=waic2,
              cpo=cpo, meanDev=meanDev))
}

### Calculate cross-validation stats  ---------------
calcCrossVal <- function(dataY, derivedProbs, nPossSpp, nSites, nSurveys){
  require(verification)  # to calc AUC statistic
  nMCMC <- dim(derivedProbs[[1]])[1]
  print(nMCMC)
  ## Calc the log score:
  devS <- calcLik(derivedProbs, dataY, nMCMC=nMCMC, nSites=nSites)
  LogScore <- mean(-2 * apply( log(devS + 0.00001), 2, sum)) 
  
  
  ## Then calculate AUC stat and Brier score:
  Brier <- list()
  AUC <- vector(length=nMCMC)
  for (s in 1:nMCMC){
    AUC[s] <- roc.area(as.vector(dataY), 
                       as.vector(derivedProbs$integrated.probs[s, , , ]))$A
    Brier[[s]] <- dataY * (1 - derivedProbs$integrated.probs[s, , , ])^2 + 
      (1 - dataY) * (derivedProbs$integrated.probs[s, , , ]) ^ 2
  }
  #  summary(AUC) 
  AUCmean <- mean(AUC)
  
  # Brier score:
  Briermean <- mean( unlist( lapply(Brier, sum, na.rm=T) ) )
  
  return(c(LogScore=LogScore, AUC=AUCmean, Brier=Briermean))
}


# end of file. ----------------------------