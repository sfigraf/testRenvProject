
library(ggplot2)  # for plots
library(vcd)  # for better colors in plots
library(reshape)  # for melt() function

# setwd("~/Dropbox/Fish_PostDoc/ArkansasR/ArkModels")  # Fill in as appropriate
#SG: not loading .rdata bc it doesn't work as well in Rmarkdown. if all scripts ran in succession, shouldn't be a need to call .Rdata files
# base::load("~/testRenvProject/Output_Files/4_RunModels/ArkInSampleModels.RData")
#source("~/testRenvProject/Source_Files/ArkFunctions.R")

## Best Model = MIDDLE-2
mod.no <- 25
mod <- out[[mod.no]]  # plots of medium model for now.
nPsiCoef <- nPsiList[[mod.no]]
nDetCoef <- nDetList[[mod.no]]


### Double-check convergence plots: -------------------

library(mcmcplots)  # for convergence plots

# denplot(mod, parms = c("muPsiCoef", "muDetCoef","sigmaPsiCoef", "sigmaDetCoef"))
# traplot(mod, parms = c("muPsiCoef", "muDetCoef", "sigmaPsiCoef", "sigmaDetCoef"))
# mcmcplot(mod, parms = c("muPsiCoef", "muDetCoef","sigmaPsiCoef", "sigmaDetCoef")) #,
         # dir="~/Dropbox/Fish_PostDoc/MultiSpp/write_up/MCMCplots/")
# caterplot(mod, parms = c("muPsiCoef", "muDetCoef","sigmaPsiCoef", "sigmaDetCoef"))


### Calculate deviance residuals -----------------------

natives[17] <- "RCS"
mod.probs <- derivePsiDetect(mod, expressPsi[[mod.no]], expressDet[[mod.no]], 
                             nPossSpp, nSites, nSurveys)
# mod.resids <- calcLik(mod, Y, mod.probs$integrated.probs, maxSurveys=max(nSurveys))
mod.resids <- calcLik(mod.probs, Y, nMCMC=mod$BUGSoutput$n.sims,
                      nSites=nSites, nSpp=nPossSpp)

d.resids <- -1 * log(mod.resids)
#str(d.resids)  # n.iter X n.spp X n.sites X nsurveys
devIJ <- aperm(d.resids, c(2, 1, 3))# apply(d.resids, c(1, 2, 3), sum, na.rm=T)
dimnames(devIJ)[[2]] <- natives # list(2:3, natives, sites$SiteID)
dimnames(devIJ)[[3]] <- sites$SiteID # 
d.residsWide <- as.data.frame(devIJ)
d.residsLong <- stack(d.residsWide)
d.residsLong$spp <- substr(d.residsLong$ind, 1, 3)
d.residsLong$site <- substring(d.residsLong$ind, 5)
# head(d.residsLong)
d.residsMedian <- with(d.residsLong, tapply(values, list(site, spp), mode.est))
d.residsMedian <- melt(d.residsMedian, varnames=c("site", "spp") )
d.residsMedian$values2 <- cut(d.residsMedian$value, 
                              breaks=c(-Inf, 2, 5, 8, 10, Inf), right=FALSE)
                              #breaks=c(-Inf, -1:10, Inf), right=FALSE)
names(d.residsMedian)[3:4] <- c("deviance", "deviance.cut")



### Plots of residuals  -------------------------------

setwd(file.path(here(), "Output_Files/6_BestModelResidPlots"))  # Fill in as appropriate

# Summarize deviances for plotting.

# Plot of all residuals
g <- ggplot(d.residsMedian, aes(spp, as.factor(site)) )
g.dev <- g + geom_tile(aes(fill=deviance.cut), color="white") + 
  scale_fill_brewer(palette = "YlGnBu")
g.dev
ggsave("Deviances_all.jpg", width=7, height=10, units="in")
#
# check out the sites with high deviances:
# subset(sites, sites$SiteID %in% c(44, 73, 100, 103, 116))
#  All have > 4 Surveys per sites.  Site 100= 17 surveys.

#  Without the cuts 
gd2 <- g + geom_tile(aes(fill=deviance), color="white") + 
  scale_fill_gradient2(midpoint=3, low="white")
gd2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Residuals by species
g.spp <- ggplot(d.residsMedian, aes(spp, y=deviance)) + geom_boxplot()
g.spp
# g.spp <- qplot(spp, deviance, data=d.residsMedian, geom="boxplot", alpha = I(1/5)) + theme_bw(base_family= 'Helvetica') # qplot was removed from package, next line uses the ggplot
g.spp <- ggplot(data=d.residsMedian, aes(x=spp, y=deviance)) + geom_boxplot(alpha = I(1/5)) + theme_bw(base_family= 'Helvetica')
g.spp  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Deviances_spp.jpg", width=9, height=7, units="in")
# tmp <- subset(d.residsMedian, deviance < 10)  #30)
# g.spp %+% tmp
#

# Residuals by site
d.residsMedian$SiteBreaks <- ifelse(d.residsMedian$site <= 16, 1, 0)
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 16 & d.residsMedian$site <= 32)] <- 2
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 32 & d.residsMedian$site <= 48)] <- 3
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 48 & d.residsMedian$site <= 64)] <-  4
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 64 & d.residsMedian$site <= 80)] <-5
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 80 & d.residsMedian$site <= 96)] <-6
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 96 & d.residsMedian$site <= 112)] <-7
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 112 & d.residsMedian$site <= 128)] <-8
d.residsMedian$SiteBreaks[which(d.residsMedian$site > 128 & d.residsMedian$site <= 180)] <-9
# g.sites <- qplot(as.factor(site), deviance, data=d.residsMedian, geom="boxplot", alpha=I(1/5))
g.sites <- ggplot(data=d.residsMedian, aes(x=as.factor(site), y=deviance)) + geom_boxplot(alpha=I(1/5))
tmp <- subset(d.residsMedian, deviance < 30)
g.sites +   facet_wrap( ~ SiteBreaks, scales="free_x") + theme_bw(base_family= 'Helvetica')  #
# g.sites %+% tmp +   facet_wrap( ~ SiteBreaks, scales="free_x") + theme_bw(base_family= 'Helvetica')
ggsave("Deviances_sites.jpg", width=8, height=6, units="in")
#

### ----------------
### Plots of residuals versus covariates
### ----------------


# Residuals versus covariates:
# Look for patterns
resid.cov <- merge(d.residsMedian, sites, by.x="site", by.y="SiteID")
# dim(resid.cov)
# names(resid.cov)
wCov <- resid.cov
with(wCov, plot(ELEV, deviance))  #
with(wCov, plot(X, deviance))  #
with(wCov, plot(Y, deviance))  # slightly higher deviance at higher Y
with(wCov, plot(jitter(CROPS), deviance))# slightly higher dev with less crops
with(wCov, plot(jitter(DVLPD), deviance))# high dev at low DVLPD
with(wCov, plot(jitter(WTLNDS), deviance))# dev is spread out
par(mar=c(2,2,2,2), mfrow=c(2, 1))
with(wCov, plot(jitter(aboveP), deviance))#
with(wCov, plot(jitter(belowJM), deviance))#
with(wCov, plot(jitter(UNCONNECT), deviance))  #
with(wCov, plot(jitter(SIZE), deviance))  #higher dev at size==4?
par(mfrow=c(1,1))
#



# end of file.
