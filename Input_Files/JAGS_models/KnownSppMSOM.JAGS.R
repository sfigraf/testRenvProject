model {   
  ### PRIORS 
  for(m in 1:nPsiCoef){
    muPsiCoef[m] ~ dnorm(0, muPsiCoef.tau)    
    halfsigmaPsiCoef[m] ~ dt(0, t.tau, 1)  # dhalfcauchy(2.25)  # half-Cauchy distribution
    sigmaPsiCoef[m] <- abs(halfsigmaPsiCoef[m])
    tauPsiCoef[m] <- pow(sigmaPsiCoef[m], -2)
  }
  for(m in 1:nDetCoef){
    muDetCoef[m] ~ dnorm(0, muDetCoef.tau)
    halfsigmaDetCoef[m] ~ dt(0, t.tau, 1)  # dhalfcauchy(2.25)   # half-Cauchy distribution
    sigmaDetCoef[m] <- abs(halfsigmaDetCoef[m])
    tauDetCoef[m] <- 1 / (sigmaDetCoef[m] * sigmaDetCoef[m])
  }
  
  ### LIKELIHOOD
  for (i in 1:nPossSpp) {
    for (m in 1:nPsiCoef){
      psiCoef[i, m] ~ dnorm(muPsiCoef[m], tauPsiCoef[m])
    }
    for (m in 1:nDetCoef){
      detCoef[i, m] ~ dnorm(muDetCoef[m], tauDetCoef[m])
    }
 
    for (j in 1:nSites) {
      logit(psi[i, j]) <- psiCoef[i, 1] + psiCoef[i, 2] * Psi.elev[j] 
        #         psiCoef[i, 3] * Psi.I25[j] + psiCoef[i, 4]*Psi.elevI25[j] +
        #         psiCoef[i, 5]*Psi.dvlpd[j] + psiCoef[i, 6]*Psi.PondY[j] +
        #         psiCoef[i, 7]*Psi.PondInt[j]
      Z[i, j] ~ dbern(psi[i, j])
        
      for (k in 1:nSurveys[j]) {
        logit(detect[i, j, k]) <- detCoef[i, 1] + detCoef[i, 2]*X.2009[j, k] +
          detCoef[i, 3]*X.2010[j, k] + detCoef[i, 4]*X.2011[j, k] +
          detCoef[i, 5]*X.2012[j, k] + detCoef[i, 6]*X.2013[j, k] +
          detCoef[i, 7]*X.elev[j, k] # + detCoef[i, 8]*X.rca[j, k]  #X.PondY[j, k]  + 
          #           detCoef[i, 5]*X.PondInt[j, k] + detCoef[i, 6]*X.elev[j, k]  + 
          #           detCoef[i, 7]*X.rca[j, k] + detCoef[i, 8]*X.crops[j, k]  + 
          #           detCoef[i, 9]*X.maxct[j, k] + detCoef[i, 10]*X.yday[j, k] 
        true.detect[i, j, k] <- detect[i, j, k] * Z[i, j]
        Y[i, j, k] ~ dbern(true.detect[i, j, k])   
        Ypred[i, j, k] ~ dbern(true.detect[i, j, k])
        score01.all[i, j, k] <- abs(Y[i, j, k] - Ypred[i, j, k])
      }
      score01ij[i, j] <- sum(score01.all[i, j, 1:nSurveys[j]])  ## 0-1 Loss function
    }          
  }

}

