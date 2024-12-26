
## This file takes the SurveyDesignStreams,
##  SurveyDesignsStreamsAsPoints, and the snappedSites
##   files from arcGIS (i.e., the stream network and sites)
##  and organizes them for EDA and R model-fitting.

## INPUTS:
## SurveyDesignStreams.dbf (shapefile) (From ArcGIS, this file doesn't change unless sampling frame and covariates change)
## SurveyDesignStreamsAsPoints.dbf (shapefile) (From ArcGIS, this file doesn't change unless sampling frame and covariates change)
## snappedSurveys.dbf (This is from ArcGIS.  This file is a new file each year.  Has the updated sample locations and data each year.)
## SampOccasions.csv file

## OUTPUTS
## ArkData.RData R workspace
# (all_sites, pred_sites, sites, SampOccasions,
#  SurveyArray, MethodArray, PassNoArray,
#  natives, nativeSppIndex, 
#  site.means, site.sds, occ.means, occ.sds,
#  samp.dat = long format of all survey data)

# Other output: EDA tables and figures

## LIBRARIES:
library(foreign) # to read and write .dbf files
library(xtable)  # to easily create LaTex tables
library(RColorBrewer)  # for color choices
library(ggplot2) # plotting functions
library(sf) # new (2024) spatial data packages
library(sp) # older spatial data packages
library(SDMTools) # to add scale bar to maps
#SG: needs rtools to install
#install.packages('https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz')
library(OpenStreetMap) # to get an underlying satellite image?

# Set the directory to input files:
setwd("~/testRenvProject/")

### STEP 1: Organize site-level info and covariates --------

## Read in ArcGIS Stream Network info 
streams <- read.dbf("Input_Files/SurveyDesignStreams.dbf")
# str(streams)
## Change INT/PER to more proper Connect/Unconnect
streams$StreamType <- ifelse(streams$StreamType=="INT", "UNCONNECT", "CONNECT")

stream_pts <- read.dbf("Input_Files/SurveyDesignStreamsAsPoints.dbf")
# str(stream_pts)
# Avoid name overlap with streams object
stream_pts <- subset(stream_pts, select=-OBJECTID)

# Add stream info to the stream_pts file
all_sites <- merge(stream_pts, streams)
# dim(all_sites)
# make sure dim = 47848 = dim of stream_pts object 
# summary(all_sites)
# rm(streams, stream_pts)
#
## Turn Land Cover Counts into Proportions:
all_sites$CROPS <- all_sites$cropsCount / 1256
all_sites$DVLPD <- all_sites$dvlpdCount / 314
all_sites$WTLNDS <- all_sites$wtlndsCoun / 314

## Create the rest of the covariates in proper form:
all_sites$UNCONNECT <- ifelse(all_sites$StreamType=="UNCONNECT", 1, 0)
all_sites$log_grad <- log(all_sites$gradient + 0.05)

(site.means <- apply(all_sites[, c("UTMX", "UTMY", "gradient", "log_grad",
                                   "elevation", "cropsCount", "dvlpdCount",
                                   "wtlndsCoun")], 2, mean, na.rm=T))
(site.sds <- apply(all_sites[, c("UTMX", "UTMY", "gradient", "log_grad",
                                 "elevation", "cropsCount", "dvlpdCount",
                                 "wtlndsCoun")], 2, sd, na.rm=T))
all_sites$Y <- (all_sites$UTMY - site.means["UTMY"]) / site.sds["UTMY"]
all_sites$X <- (all_sites$UTMX - site.means["UTMX"]) / site.sds["UTMX"]
all_sites$ELEV <- (all_sites$elevation - site.means["elevation"]) / site.sds["elevation"]
all_sites$GRAD <- all_sites$log_grad # (all_sites$log_grad - site.means["log_grad"]) / site.sds["log_grad"]

all_sites$SIZE <- all_sites$StreamSize
all_sites$RESERVOIR <- as.factor(all_sites$ReservoirL)
all_sites$FTN <- with(all_sites, as.numeric(FountainCr) - 1)
# all_sites$FTNLAT <- with(all_sites, (as.numeric(FountainCr) - 1) *LAT)
all_sites$ELEV2 <- all_sites$ELEV ^ 2
all_sites$MAIN <- with(all_sites, as.numeric(MAIN) - 1)
all_sites$PURG <- with(all_sites, as.numeric(PURG) - 1)
#

# Look at covariate correlations
# round(cor(all_sites[, c("X", "Y", "elevation", "gradient",
#                         "CROPS", "DVLPD", "WTLNDS")]), 2)
tmp <- all_sites
tmp$StreamType <- as.numeric(as.factor(all_sites$StreamType))
tmp$RESERVOIR <- as.numeric(all_sites$RESERVOIR)
tmp$FTN <- as.numeric(all_sites$FTN)
tmp$PURG <- as.numeric(all_sites$PURG)
tmp$MAIN <- as.numeric(all_sites$MAIN)
# round(cor(tmp[, c("UTMX", "UTMY", "elevation", "gradient",
#                   "CROPS", "DVLPD", "WTLNDS",
#                   "StreamType", "StreamSize", 
#                   "RESERVOIR", "FTN", "PURG", "MAIN")]), 2)
# UTMX and elevation have cor=-0.81
# FTN and DVLPD have a cor=0.58
# StreamSize and WTLNDS have cor=0.57
# StreamSize and StreamType have cor=0.42
# StreamSize and MAIN have cor=0.75
# UTMY and PURG have cor=-0.59
# WTLNDS and MAIN have cor=0.56
# CROPS and UTMX have cor=0.45
# elevation and CROPS have cor=-0.47
rm(tmp)

## Stream info for predictions 
# dim(all_sites)
# Choose a subset of all sites for predictions and 
#   design criterion:
pred_sites <- all_sites[seq(1, nrow(all_sites), 100), ]
# dim(pred_sites)
# 479 sites to choose from.


### STEP 2: Add SampOccasion-level info and covariates --------

# Import the snapped locations to get correct covariate info:
#survey data snapped to correct location in the river (
#one that works is snapped data from july 6 2016
#new snaped surveys that brian added doesn't work
#original 2016 snapped surveys
#columns 20
#
snapped_sites <- read.dbf("Input_Files/snappedSurveys.dbf")

## After comparing data to stream network, need to edit one site's location:
snapped_sites[snapped_sites$SurveyID==43945, "UTMX"] <- 620287
snapped_sites[snapped_sites$SurveyID==43945, "UTMY"] <- 4167592
snapped_sites[snapped_sites$SurveyID==45532, "UTMX"] <- 620287
snapped_sites[snapped_sites$SurveyID==45532, "UTMY"] <- 4167592
#
#names(snapped_sites)
# survey ID changed to Species Code
#SG: 
#
names(snapped_sites)[20] <- "SpeciesCode"

## Remove some sites that are still hanging around and shouldn't be:
## Sites are repeated counts of other surveys, but slightly diff coords:
#SG: changedSurveyID here to Species Code, but also code 25784 doesn't seem to even be in here
snapped_sites <- subset(snapped_sites, !(SurveyID %in% c(25784, 23897)))
# dim(snapped_sites)  # 3197 obs.
#

## IF YOU GET TO END AND NEED TO REMOVE MORE SAMP OCCASIONS -------------
# Run the "Ark_OrgSampOccForExport.R" on the new/updated data spreadsheet.
# tmp1 <- include_dat[, c("mySurveyID", "Catch", 
#                         "SpeciesCodee", "sameStatus")]
# tmp2 <- snapped_sites[, c("WaterName", "SampleDate", "UTMX", "UTMY",
#                           "Source", "mySurveyID")]
# tmp2$mySurveyID <- as.character(tmp2$mySurveyID)
# tmp <- merge(tmp1, tmp2, all=T)
# tmp <- unique(tmp)
# tmp$mySurveyID <- as.factor(tmp$mySurveyID)
# summary(tmp)
# dim(tmp)
# dim(include_dat)  # should match, then:
# snapped_sites <- tmp
# -----------------------------------

# Add covariate info associated with each sampled site: 
# names(snapped_sites)
# names(all_sites)
#SG: should samp_dat have data here?
#left_join works, merge() doesnt...
samp_dat <- base::merge(all_sites, snapped_sites)
samp_dat <- merge(snapped_sites, all_sites)
#samp_dat <- left_join(all_sites, snapped_sites)

# dim(snapped_sites) #3188  3197...
# dim(samp_dat) #3188  3197...
# Check any sites that don't get merged:
#snapped_sites[!(snapped_sites$UTMX %in% samp_dat$UTMX), ]
# eh, OK.
# rm(snapped_sites)
#
# Only keep relevant columns:
# sort(names(samp_dat))
samp_dat <- unique(samp_dat[, c("UTMX", "UTMY", "mySurveyID",
                                "SampleDate",
                                #"Location", "Elevation", #"StationLength", "AvgWidth",
                                #"StationCode", "SiteType", "StreamType",
                                "Catch",  "sameStatus", "SpeciesCode", 
                                "Source", # whether site is USGS
                                "reachid", "pointid", "X", "Y", 
                                "UNCONNECT", "GRAD", "ELEV","ELEV2", 
                                "SIZE", # "FTNLAT", 
                                "CROPS", "DVLPD", "WTLNDS", "RESERVOIR",
                                "FTN", "MAIN", "PURG")] )
# dim(samp_dat)  # 3188.  
## Add survey-specific covariates ----------------------
# summary(samp_dat)
# Add YEAR as a factor covariate
# head(strptime(as.character(samp_dat$SampleDate), 
#            format="%d-%b-%Y")$yday)
# head(samp_dat$SampleDate)
# head(strptime(as.character(samp_dat$SampleDate), 
#               format="%d-%b-%Y")$year)
#SG: All my sample dates are NA
samp_dat$year <- strptime(as.character(samp_dat$SampleDate), 
                          format="%d-%b-%Y")$year - 100 + 2000
samp_dat$YEAR <- as.factor(samp_dat$year)
#
samp_dat$month <- strptime(as.character(samp_dat$SampleDate), 
                           format="%d-%b-%Y")$mon + 1
samp_dat$yday <- strptime(as.character(samp_dat$SampleDate), 
                          format="%d-%b-%Y")$yday
# subset(samp_dat, is.na(year))
# survey 39699 = No Data.  OK, Comments say spot sample for ARD.
samp_dat <- subset(samp_dat, !is.na(year))

samp_dat$YDAY <- (samp_dat$yday - mean(samp_dat$yday)) / sd(samp_dat$yday)
samp_dat$YDAY2 <- samp_dat$YDAY ^ 2

samp_dat$METHOD <- NA
samp_dat$METHOD[grep("E", samp_dat$sameStatus)] <- "E"
samp_dat$METHOD[grep("S", samp_dat$sameStatus)] <- "S"
samp_dat$METHOD[grep("D", samp_dat$sameStatus)] <- "D"
# unique(subset(samp_dat[, c("mySurveyID", "sameStatus", "METHOD")], is.na(METHOD))  )
# should be 0 rows.

samp_dat$PASSNO <- NA
samp_dat$PASSNO[grep("1", samp_dat$sameStatus)] <- 1
samp_dat$PASSNO[grep("2", samp_dat$sameStatus)] <- 2
samp_dat$PASSNO[grep("3", samp_dat$sameStatus)] <- 3
samp_dat$PASSNO[grep("4", samp_dat$sameStatus)] <- 4
samp_dat$PASSNO[grep("5", samp_dat$sameStatus)] <- 5
samp_dat$PASSNO[grep("6", samp_dat$sameStatus)] <- 6
# summary(samp_dat$PASSNO)
# subset(samp_dat[, c("mySurveyID", "SpeciesCode", "sameStatus", "METHOD", "PASSNO")], is.na(PASSNO))
# should be 0 rows

# samp_dat[samp_dat$mySurveyID %in% samp_dat$mySurveyID[grep("DIP", samp_dat$sameStatus)], ]
samp_dat$PASSNO[grep("DIP", samp_dat$sameStatus)] <- 1
#
# names(samp_dat)
# levels(as.factor(samp_dat$METHOD))
#
# summary(samp_dat)
#
# Double-check surveys that happened way early in the year:
# summary(subset(samp_dat, yday < 103))
# All from 2006-- they will eventually be removed.  Good.

## Add in total count per occasion
# length(unique(samp_dat$mySurveyID)) # 378
#SG: brin said warning is OK here
tmp <- with(samp_dat, tapply(as.numeric(as.character(Catch)), 
                             as.character(mySurveyID), sum, na.rm=T))
# head(tmp)
# length(tmp)  #378
#
tmp2 <- as.data.frame(tmp)
names(tmp2)[1] <- "TotalCt"
tmp2$mySurveyID <- row.names(tmp2)
# head(tmp2)
tmp23 <- merge(samp_dat, tmp2)
# dim(tmp23)  # 3185
# dim(samp_dat)  #3185
samp_dat <- tmp23
samp_dat$logTotalCt <- log(samp_dat$TotalCt + 1)
#
rm(tmp, tmp2, tmp23)
#




# Clean up surveys based on catch data -------------

## Clean-up the catch info:
# levels(samp_dat$SpeciesCode)
# CPK does not exist.  Should be common carp, CPP
samp_dat$SpeciesCode[samp_dat$SpeciesCode == "CPK"] <- "CPP"
# OTC does not exist.  Should be other warmwater species, OTS
samp_dat$SpeciesCode[samp_dat$SpeciesCode == "OTC"] <- "OTS"
# RSD does not exist.  Should be Red shiner, RDS
samp_dat$SpeciesCode[samp_dat$SpeciesCode == "RSD"] <- "RDS"

# Remove catches with non-specific codes:
# samp_dat <- subset(samp_dat, !(SpeciesCode %in% c("YOY", "UNK", "OTS")))

# Double-check NULL Catches:
tmp <-   subset(samp_dat, Catch == "---")
sampNULL <- subset(samp_dat, Catch == "---")$mySurveyID
# samp_dat[samp_dat$mySurveyID %in% sampNULL, ]
# ## Remove NULL Catches that don't make sense
# can't tell which to keep--throwout both
samp_dat <- samp_dat[samp_dat$mySurveyID != "718016_12-Aug-2008", ]
rm(tmp, sampNULL)
#
# levels(as.factor(as.character(samp_dat$SpeciesCode)))
# Some fish codes do not exist. Check them out.
tmp <- subset(samp_dat, SpeciesCode == "FHM")$mySurveyID
# samp_dat[samp_dat$mySurveyID %in% tmp, ]
tmp <- subset(samp_dat, SpeciesCode == "MNW")$mySurveyID
# samp_dat[samp_dat$mySurveyID %in% tmp, ]
tmp <- subset(samp_dat, SpeciesCode == "SUF")$mySurveyID
# samp_dat[samp_dat$mySurveyID %in% tmp, ]

## Check out surveys with trout and see if they should be excluded:
tmp <- samp_dat[samp_dat$SpeciesCode %in% c("HXC", "LOC", "RBT", "RXN", "TRT"), ]$mySurveyID
tmp2 <- samp_dat[samp_dat$mySurveyID %in% tmp, c("mySurveyID", "Catch", "SpeciesCode", 
                                                 "sameStatus", "ELEV")]
cbind(tmp2, elev=(site.sds["elevation"] * tmp2$ELEV) + site.means["elevation"])

# Sites from Grape Creek (471536..., 476218) are below 6000ft--keep. 
# Beaver creek(497594, 497983, )-- less than 6000 ft--keep
# Beaver creek(498477, 498743)-- greater than 6000ft-- remove.
# One site (507861...) from higher elev St. Charles-- Foutz says keep St. Charles (5/2016 mtg)
# Keep ABV Diversion, (510582...)
# Ftn Creek #4 sampling (513954...)-high elev and only trout-- Foutz says keep all Ftn Crk (5/2016 mtg)
# Keep Ftn Creek #3 (514396..., 515497...)
# Keep 517934..., caught many plains fish, few trout
# Keep 517862..., caught many plains fish, few trout
# Keep 518221..., caught many plains fish, few trout
# Keep 524655, 525007, 529933, 531810, 533299-- all on Ark at lower elev
# Keep 535464-- on Ftn Crk at lower elev
# Keep 543498-- on Purg at lower elev
tmp2 <- tmp[unlist(lapply(c("498743", "498477"),  function(e) grep(e, tmp)))]
tmp2 <- unique(tmp2)
samp_dat <- subset(samp_dat, !(mySurveyID %in% tmp2))
# dim(samp_dat) # now 3161  
rm(tmp, tmp2)
#





# Exclude surveys ------------------------------------
# Exclude USGS surveys: 
# levels(samp_dat$Source)
samp_dat <- subset(samp_dat, Source != "USGS")
# dim(samp_dat) # 2645 

## Exclude sampling before 2008 to match South Platte:
samp_dat <- subset(samp_dat, year >= 2008)
# dim(samp_dat) #2098  

## How many dipnet passes?
tmp <- subset(samp_dat, METHOD == "D")
# dim(unique(subset(tmp, select=-c(SpeciesCode, Catch)))) #14 passes
# How many dip net passes?
# dim(unique(subset(samp_dat, METHOD == "D")[, c("PASSNO", "mySurveyID")]))
# 13
## Exclude dipnet passes:
samp_dat <- subset(samp_dat, METHOD != "D")

# Object of samp occasion info only 
SampOccasions <- unique(samp_dat[, c("UTMX", "UTMY", "X", "Y", 
                                     "mySurveyID", 
                                     "year", "YEAR", 
                                     "yday", "YDAY", "YDAY2",
                                     "reachid", "pointid", #"FTNLAT", 
                                     "UNCONNECT", "GRAD", "ELEV", "ELEV2", 
                                     "SIZE",
                                     "CROPS", "DVLPD", "WTLNDS", "RESERVOIR",
                                     "FTN", "MAIN", "PURG")])
# dim(SampOccasions)  # 202 w/o USGS  
#
# mini-EDA of Samp Occasions
## How many surveys are left?
# length(unique(SampOccasions$mySurveyID))  # 201  
# above length SHOULD MATCH nrow(SampOccasions) 
# If not, check which ones are repeat:
# sort(SampOccasions$mySurveyID)
# samp_dat[samp_dat$mySurveyID == "626725_9-Oct-2013", ]
# First set has 0 elevation, otherwise counts are repeats.
# remove repeat survey
samp_dat <- samp_dat[!(samp_dat$mySurveyID == "626725_9-Oct-2013" & 
                         samp_dat$pointid == 6626), ]
# 2004 obs.  2060...
SampOccasions <- unique(samp_dat[, c("UTMX", "UTMY", "X", "Y", 
                                     "mySurveyID", 
                                     "year", "YEAR", 
                                     "yday", "YDAY", "YDAY2",
                                     "reachid", "pointid", #"FTNLAT", 
                                     "UNCONNECT", "GRAD", "ELEV", "ELEV2", 
                                     "SIZE",
                                     "CROPS", "DVLPD", "WTLNDS", "RESERVOIR",
                                     "FTN", "MAIN", "PURG")])
# dim(SampOccasions)  # now 201 w/o USGS 

## Exclude single-pass surveys?
## (Commented out-- they are included)
# SampOccasions <- multiplePass_occasions
# samp_dat <- multiplePass_dat

# Create DF of site-only info:
sites <- unique(subset(SampOccasions, 
                       select=-c(mySurveyID, year, YEAR, 
                                 yday, YDAY, YDAY2)))
# dim(sites)  #  141 
# summary(sites)
sites$SiteID <- 1:nrow(sites)

# Add SiteID to SampOccasions
SampOccasions <- merge(sites, SampOccasions) #still 201 rows

# Add SiteID to samp_dat
samp_dat2 <- merge(samp_dat, sites, all.x=T)
# summary(samp_dat2)  # if all SiteIDs are filled in, then
samp_dat <- samp_dat2
rm(samp_dat2)

# and add nSurveys column to sites:
surveys <- unique(subset(samp_dat, select=-c(SpeciesCode, Catch)))
nSurveys <- with(surveys, aggregate(sameStatus, by=list(SiteID), length))
names(nSurveys) <- c("SiteID", "nSurveys")
# How many passes/surveys in total?
# sum(nSurveys$nSurveys) # 405
SampOccasions <- merge(SampOccasions, nSurveys) # 201
sites <- merge(sites, nSurveys) # 141

# And add nSurveys per SampOccasions to SampOcc object:
SampnSurveys <- with(surveys, aggregate(sameStatus, by=list(mySurveyID), length))
names(SampnSurveys) <- c("mySurveyID", "nSurveys")
# How many passes/surveys in total?
# sum(SampnSurveys$nSurveys) # 405
SampOccasions <- merge(SampOccasions, SampnSurveys) #110...
rm(SampnSurveys)

# DOUBLE-CHECK larger nSurveys against Excel spreadsheets:
# summary(nSurveys)
subset(SampOccasions[, c("mySurveyID", "SiteID", "nSurveys")], nSurveys > 4)
rm(nSurveys)

### STEP 3:  Survey-level objects ------------------
spp <- sort(as.character(unique(samp_dat$SpeciesCode)))
# spp

natives <- c("ARD", "BBH", "STR", "CCF", "FMW", "FHC",
             "SNF", "LND", "OSF", "PKF", "PMW", "RDS",
             "SAH", "SRD", "SMM", "WHS")
(nDetected <- length(natives)) #16 native species *detected* (17 in basin)
nPossSpp <- 17

# MethodArray and PassNo Array:
(nSites <- nrow(sites))
(maxSurveys <- max(sites$nSurveys)) # max(SampOccasions$nSurveys))
methodArray <- array(NA, dim=c(nSites, maxSurveys))
passNoArray <- methodArray
surveyID.mat <- methodArray
# nrow(surveys)
# nrow should be 405 = Total number of passes
surveys$passID <- 1:nrow(surveys)
# 13 relevant survey-level covariates
X.arrayAll <- array(NA, dim=c(13, nSites, maxSurveys))
K <- vector(length=nSites)
for (j in 1:nSites){
  tmp <- subset(surveys, SiteID==sites$SiteID[j])
  K[j] <- length(tmp$mySurveyID)
  surveyID.mat[j, 1:K[j]] <- as.character(tmp$mySurveyID)
  methodArray[j, 1:K[j]] <- tmp$METHOD
  passNoArray[j, 1:K[j]] <- tmp$PASSNO
  X.arrayAll[1, j, 1:K[j]] <- tmp$YDAY
  X.arrayAll[2, j, 1:K[j]] <- tmp$YDAY2
  X.arrayAll[3, j, 1:K[j]] <- ifelse(tmp$METHOD == "S", 1, 0)
  X.arrayAll[4, j, 1:K[j]] <- tmp$PASSNO
  X.arrayAll[5, j, 1:K[j]] <- ifelse(tmp$YEAR == 2008, 1, 0)
  X.arrayAll[6, j, 1:K[j]] <- ifelse(tmp$YEAR == 2009, 1, 0)
  X.arrayAll[7, j, 1:K[j]] <- ifelse(tmp$YEAR == 2010, 1, 0)
  X.arrayAll[8, j, 1:K[j]] <- ifelse(tmp$YEAR == 2011, 1, 0)
  X.arrayAll[9, j, 1:K[j]] <- ifelse(tmp$YEAR == 2012, 1, 0)
  X.arrayAll[10, j, 1:K[j]] <- ifelse(tmp$YEAR == 2013, 1, 0)
  X.arrayAll[11, j, 1:K[j]] <- ifelse(tmp$YEAR == 2014, 1, 0)
  X.arrayAll[12, j, 1:K[j]] <- ifelse(tmp$YEAR == 2015, 1, 0)
  X.arrayAll[13, j, 1:K[j]] <- tmp$logTotalCt
}  
rm(tmp, K)

# SurveyArray:
Y <- array(NA, dim=c(nPossSpp, nSites, maxSurveys))
for (i in 1:nPossSpp){
  # fill in array with 0s for every pass as default value:
  Y[i, , ] <- 0 * passNoArray 
  tmp <- subset(samp_dat, SpeciesCode == natives[i])
  for(j in 1:nSites){
    for (k in 1:sites$nSurveys[j]) {
      if (surveyID.mat[j, k] %in% tmp$mySurveyID){
        Y[i, j, k] <- 1
      }
    }
  }
}
# str(Y)
#i = spp, j= sites, k=survey
rm(tmp, i, j, k)
## Create Z matrix of known occurrences
Ybinom <- apply(Y, c(1, 2), sum, na.rm=T)
Z <- ifelse(Ybinom > 0, 1, NA)
# t(Z)
##
### STEP 4: Add'l arrays and data frames for JAGS models -----------
sites$aboveP <- ifelse(sites$RESERVOIR=="abovePueblo", 1, 0)
sites$belowJM <- ifelse(sites$RESERVOIR=="belowJM", 1, 0)
sites$X2 <- sites$X ^ 2
sites$Y2 <- sites$Y ^ 2
site.varsAll <- t(as.matrix(sites[, c("X", "X2", "Y", "Y2",
                                      "ELEV", "ELEV2",
                                      "CROPS", "DVLPD", "WTLNDS",
                                      "UNCONNECT", "SIZE",
                                      "aboveP", "belowJM",
                                      "FTN", "MAIN", "PURG")]))
# summary(t(site.varsAll))
nSurveys <- sites$nSurveys





#
### Save created objects ----------------------------
# First remove extraneous objects:
rm(streams, stream_pts)

save.image("Output_Files/2_OrganizeNewData/ArkData.Rdata")
#














### ---------------------------------------------------------------
### Exploratory Data Analysis (EDA) tables and figures of site and survey info ----------------
setwd("Output_Files/2_OrganizeNewData/EDA")

# nrow(surveys)  # 412
# nrow(sites)  # 143
# nrow(SampOccasions) # 203
# length(unique(surveys$mySurveyID))  # 203
#
### TABLES --------------------------------------------------------

### Table of stream types sampled each year:
# as.matrix(with(sites, table(UNCONNECT)))
tmp <- as.matrix(with(SampOccasions, table(year, UNCONNECT)))
tmp2 <- rbind(tmp, t(as.matrix(with(SampOccasions, table(UNCONNECT)))))
xtmp <- xtable(tmp2, caption="The stream types sampled each year.  
               CONNECT means the stream segment was less than 2000 m from 
               the perennial stream network, and UNCONNECT means the 
              stream segment was more than 2000 m from the perennial
               stream network.")
# xtmp
sink("StreamTypeByYear.tex")
print(xtmp, caption.placement="top")
sink()


### Table of Stream SIZES sampled each year:
# as.matrix(with(sites, table(SIZE)))
tmp <- as.matrix( with(SampOccasions, table(year, SIZE)) )
tmp2 <- rbind(tmp, t(with(SampOccasions, table(SIZE))))
xtmp <- xtable(tmp2, caption="Stream sizes that were sampled each year.  
               A stream that is size 1 is an intermittent stream while 
               a stream of size 4 is the main stem of the Arkanas River, 
               below the Pueblo reservoir. The stream sizes follow the 
               Strahler stream order rules, but were based on the 
               sampling frame used in this design.")
# xtmp
# Check out 2012, size=3. lots of samples!
# SampOccasions[SampOccasions$year== "2012" & SampOccasions$SIZE == 3, 
#               c("SiteID", "UTMX", "UTMY", "mySurveyID", "nSurveys")]
sink("StreamSizeByYear.tex")
print(xtmp, caption.placement="top")
sink()

### Table of Sites location between Reservoirs, by year:
# as.matrix(with(sites, table(RESERVOIR)))
tmp <- rbind(as.matrix(with(SampOccasions, table(year, RESERVOIR)))[, c(1, 3, 2)],
             as.matrix(t(with(SampOccasions, table(RESERVOIR))))[, c(1, 3, 2)])
xtmp <- xtable(tmp, caption="Locations of sampling occasions each year, 
               broken down by whether the site would drain into the 
               Pueblo reservoir (abovePueblo), below the Pueblo 
               but above John Martin (between), or below John Martin 
               reservoir (belowJM). Altogether, 7 of the sites were 
               above Pueble reservoir, 106 were between the reservoirs,
               and 28 were below John Martin.")
# xtmp
sink("ReservoirByYear.tex")
print(xtmp, caption.placement="top")
sink()

### Table of Fountain Creek, by year:
# as.matrix(with(sites, table(FTN)))
tmp <- rbind(as.matrix(with(SampOccasions, table(year, FTN))),
             as.matrix(t(with(SampOccasions, table(FTN)))))
xtmp <- xtable(tmp, caption="Locations of sampling occasions each year, 
               broken down by whether the site was located in Fountain 
               Creek or not. NO means the site was not on Fountain Creek 
               and YES means the site was on Fountain Creek or drained
               into Fountain Creek. Altogether, 31 of the sites were on 
               Fountain Creek or drained into it, and 110 sites were not on
               Fountain Creek.")
# xtmp
sink("FtnCrkByYear.tex")
print(xtmp, caption.placement="top")
sink()

### Table of Purgatoire, by year:
# as.matrix(with(sites, table(PURG)))
tmp <- rbind(as.matrix(with(SampOccasions, table(year, PURG))),
             as.matrix(t(with(SampOccasions, table(PURG)))))
xtmp <- xtable(tmp, caption="Locations of sampling occasions each year, 
               broken down by whether the site was located in the Purgatoire 
              River or not. NO means the site was not in the Purgatoire River 
               and YES means the site was in the Purgatoire River or drained
               into it. Altogether, 8 of the sites were on 
               the Purgatoire River or drained into it, and 133 sites were not on
               the Purgatoire River.")
# xtmp
sink("PURGByYear.tex")
print(xtmp, caption.placement="top")
sink()

### Table of Main Stem, by year:
# as.matrix(with(sites, table(MAIN)))
tmp <- rbind(as.matrix(with(SampOccasions, table(year, MAIN))),
             as.matrix(t(with(SampOccasions, table(MAIN)))))
xtmp <- xtable(tmp, caption="Locations of sampling occasions each year, 
               broken down by whether the site was located on the main stem of 
              the Arkansas River or not. NO means the site was not on the 
              main stem of the Arkansas River 
               and YES means the site was on the 
              main stem of the Arkansas River. Altogether, 23  sites were on 
               the main stem of the Arkansas River or drained into it, 
and 118 sites were not on
               the main stem of the Arkansas River.")
# xtmp
sink("MainByYear.tex")
print(xtmp, caption.placement="top")
sink()


### Table of Fish Species Detections per Stream Type:
dat.tmp <- subset(samp_dat, SpeciesCode %in% natives)
dat.tmp$SpeciesCode <- as.character(dat.tmp$SpeciesCode)
tmp <- as.matrix(with(dat.tmp, table(SpeciesCode, UNCONNECT)))
# tmp
# how many CONNECT and UNCONNECT passes were there?
# dim(subset(surveys, UNCONNECT==0))  # 331
# dim(subset(surveys, UNCONNECT==1))  # 74
xtmp <- xtable(tmp, caption="The number of \textit{passes} in which each 
               native fish species was detected, broken down by whether the survey 
               occurred on a connected (CONNECT) or unconnected (UNCONNECT) \
               stream. There were 331 passes in the connected streams and 
               74 passes in the unconnected streams.")
sink("FishDetectionsByStreamType.tex")
print(xtmp, caption.placement="top")
sink()
#






### MAPS ---------------------------------------------
setwd("~/testRenvProject//ArcGIS_files/AsShapefiles")

## Read in the shapefiles:
# CO.rg <- readOGR(".", "Colorado_Boundary")
CO.rg <- read_sf(".", "Colorado_Boundary")
# elev1980.rg <- readOGR(".", "elev1980")
elev1980.rg <- read_sf(".", "elev1980")
# HUC.rg <- readOGR(".", "HUC4_AR")
HUC.rg <- read_sf(".", "HUC4_AR")
# I25.rg <- readOGR(".", "I25_Projected")
I25.rg <- read_sf(".", "I25_Projected")
# cities.rg <- readOGR(".", "ArkCities")
cities.rg <- read_sf(".", "ArkCities")
# noNullstreams.rg <- readOGR(".", "FlowlinesArk_NoNull")
noNullstreams.rg <- read_sf(".", "FlowlinesArk_NoNull")
# streams.rg <- readOGR(".", "SurveyDesignStreams")
streams.rg <- read_sf(".", "SurveyDesignStreams")
# streamAsPts.rg <- readOGR(".", "SurveyDesignStreamsAsPoints")
streamAsPts.rg <- read_sf(".", "SurveyDesignStreamsAsPoints")

# Check that all projections match
# print(proj4string(CO.rg))
# st_crs(CO.rg)$proj4string


# print(proj4string(elev1980.rg))
# st_crs(elev1980.rg)$proj4string

# elev1980.nad83 <- spTransform(elev1980.rg, CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

elev1980.nad83 <- st_transform(elev1980.rg, crs="+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# # print(proj4string(elev1980.nad83))
# st_crs(elev1980.nad83)$proj4string
# # print(proj4string(HUC.rg))
# st_crs(HUC.rg)$proj4string
# # print(proj4string(I25.rg))
# st_crs(I25.rg)$proj4string
# # print(proj4string(noNullstreams.rg))
# st_crs(noNullstreams.rg)$proj4string
# # print(proj4string(streams.rg))
# st_crs(streams.rg)$proj4string
# # print(proj4string(streamAsPts.rg))
# st_crs(streamAsPts.rg)$proj4string
# # print(proj4string(cities.rg))
# st_crs(cities.rg)$proj4string

# Map of sampling locations. Color-coded by ifelse PASSNO > 1
# generate a simple map showing all  layers
sites.sp <- sites # SampOccasions
# sites.sp$colPASSES <- ifelse(sites.sp$nSurveys > 1, "black", "gray")

# These functions need library(sp)
coordinates(sites.sp) <- c("UTMX", "UTMY")
proj4string(sites.sp) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
# summary(sites.sp)

#plot(streams.rg, col="blue", lwd=2.0)
#SG: had slighty change path
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/Sites.png", 
    width=1024, height=768)
par(mar=c(0.01, 0.01, 0.01, 0.01))
# plot(HUC.rg, col=gray(0.9), border=gray(0.9), mar=c(0.1, 0.1, 0.1, 0.1))
plot(st_geometry(HUC.rg), col=gray(0.9), border=gray(0.9), mar=c(0.1, 0.1, 0.1, 0.1))

# plot(elev1980.nad83, add=T, col=gray(0.7), border=gray(0.7) )
plot(st_geometry(elev1980.nad83), add=T, col=gray(0.7), border=gray(0.7) )

# plot(CO.rg, border="black", add=T, lwd=3)
plot(st_geometry(CO.rg), border="black", add=T, lwd=3)

#lines(I25.rg, col="red", cex=0.8)
# lines(noNullstreams.rg, col=gray(0.5), lwd=0.8)
plot(st_geometry(streams.rg), add=T, col="dodgerblue3", lwd=2.0)

# add cities (circled)
# points(cities.rg, cex=2.8, col="limegreen", pch=6, lwd=5)
# add labels (using trial and error for placement)
#text(cities.rg, labels=as.character(cities.rg$NAME), col="black",
#     cex=1.7, font=2, offset=0.5, adj=c(1,2))
# Add the sampling locations:
points(subset(sites.sp, nSurveys = 1), col="black", pch=17, cex=1)
points(subset(sites.sp, nSurveys > 1), col="black", pch=8, cex=1.2)
# Add scale bar
#SG: this comes from sdmTools archived package
Scalebar(646000, 4105000, 100000, unit = "km", 
         scale = 0.01, t.cex = 1.5)
# Add title
mtext("Arkansas River Basin", line=-2, cex=2)
dev.off()
#
rm(sites.sp)



# Map of elevations in the stream


# pts.1 <- list("sp.points",
#               SpatialPoints(
#                 cbind(
#                   gwss.mc.out.1.sig[,1],
#                   gwss.mc.out.1.sig[,2])
#               ), cex=2, pch="+", col="black")
# map.layout.1 <-list(pts.1)


huc <- list(HUC.rg, fill = gray(0.9))
elev <- as.list(st_set_geometry(elev1980.nad83, NULL))#, fill=gray(0.7))

# library(plyr)
# dlply(elev,1,c)

# distance bar
l3 = list("SpatialPolygonsRescale", layout.scale.bar(), 
          offset = c(646000, 4104000), scale = 100000, fill=c("transparent","black"))
# zero of distance bar
l4 = list("sp.text", c(646000, 4116000), "0", cex=0.8)
# 500 of distance bar
l4b = list("sp.text", c(696000, 4116000), "500", cex=0.8)
# 1000 km of distance bar
l5 = list("sp.text", c(746000, 4116000), "1000 km", cex=0.8)

png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/elevations.png", 
    width=1024, height=768)

# pts.1 <- list(l3, l4, l4b, l5,elev)

# spplot(as_Spatial(streamAsPts.rg["elevation"]),
#        sp.layout=list(l3,l4,l4b,l5), # list(l3, l4, l4b, l5,elev), #, huc, elev), 
#        main=list(label="Elevation (m)", cex=1),
#        colorkey=T, cex=1, par.settings=list(fontsize=list(text=36)))


# plot(st_geometry(elev1980.nad83))
# points(streamAsPts.rg["elevation"])
# 
# 
# # using the sf package
# plot(elev1980.nad83["geometry"])
# plot(streamAsPts.rg["elevation"], key.pos = 4)
# plot(elev1980.nad83["geometry"], add=T)
# # add scale bar with text 
# Scalebar(646000, 4000000, 100000, unit = "km", 
#          scale = 0.01, t.cex = 1.5)
# 
# plot(st_geometry(elev1980.nad83), add=T)
# plot(l3, add=T)

# using ggplot and sf
ggplot(data=elev1980.nad83["geometry"]) + geom_sf(fill = "white") +
  geom_point(data=streamAsPts.rg, aes(y=UTMY, x=UTMX, color=elevation)) +
  ggspatial::annotation_scale(location = 'br', height=unit(1, "cm"), dist_unit = "km", dist = 200) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        # panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title = element_text(size=30),
        legend.text = element_text(size = 30),
        plot.title = element_text(hjust = 0.5, size=30),
        legend.key.size = unit(2, 'cm')) +
  scale_colour_gradientn(colors=rainbow(5), trans = "reverse") +
  labs(title = "Elevation (m)", x = "", y = "", color = "Elevation")


# Scalebar(646000, 4000000, 100000, unit = "km", 
#          scale = 0.01, t.cex = 1.5)


dev.off()

# Map of stream sizes
all_sites.sp <- all_sites
all_sites.sp$PURG <- as.factor(ifelse(all_sites.sp$PURG==0, "NO", "YES"))
all_sites.sp$FTN <- as.factor(ifelse(all_sites.sp$FTN==0, "NO", "YES"))
all_sites.sp$MAIN <- as.factor(ifelse(all_sites.sp$MAIN==0, "NO", "YES"))
all_sites.sp$UNCONNECT <- as.factor(ifelse(all_sites.sp$UNCONNECT==0, 
                                           "CONNECTED", "UNCONNECTED"))
coordinates(all_sites.sp) <- c("UTMX", "UTMY")
proj4string(all_sites.sp) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
all_sites.sp$streamF <- as.factor(all_sites.sp$StreamSize)

png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/streamsize.png",
    width=1024, height=768)
spplot(all_sites.sp["streamF"],
       # sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Strahler-type stream size", cex=1),
       col.regions=rev(heat.colors(6)),
       #      col.regions=brewer.pal(4, "YlGn"), cex=1, 
       key.space="right", #auto.key=list(cex=3),
       par.settings=list(fontsize=list(text=36)))
dev.off()


## Maps of land cover covariates:
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/crops.png",
    width=1024, height=768)
spplot(all_sites.sp["CROPS"],
       # sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Proportion cropland", cex=1),
       colorkey=T, cex=1, par.settings=list(fontsize=list(text=36)))
dev.off()
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/wetlands.png",
    width=1024, height=768)
spplot(all_sites.sp["WTLNDS"],
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Proportion wetlands", cex=1),
       colorkey=T, cex=1, par.settings=list(fontsize=list(text=36)))
dev.off()
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/dvlpd.png",
    width=1024, height=768)
spplot(all_sites.sp["DVLPD"],
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Proportion developed", cex=1),
       colorkey=T, cex=1, par.settings=list(fontsize=list(text=36)))
dev.off()

## Map of reservoir locations
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/reservoirs.png",
    width=1024, height=768)
spplot(all_sites.sp["RESERVOIR"],
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Reservoir drainages", cex=1),
       col.regions=rev(heat.colors(6)),
       key.space="right",
       par.settings=list(fontsize=list(text=36)))
dev.off()
## Map of Fountain Creek
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/fountain.png",
    width=1024, height=768)
spplot(all_sites.sp["FTN"],
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Fountain Creek", cex=1),
       col.regions=rev(heat.colors(6)),
       key.space="right",
       par.settings=list(fontsize=list(text=36)))
dev.off()
## Map of Purgatoire R
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/purgatoire.png",
    width=1024, height=768)
spplot(all_sites.sp["PURG"],
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Purgatoire River", cex=1),
       col.regions=rev(heat.colors(6)),
       key.space="right",
       par.settings=list(fontsize=list(text=36)))
dev.off()
## Map of Main Stem
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/mainstem.png",
    width=1024, height=768)
spplot(all_sites.sp["MAIN"],
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Main stem of Arkansas River", cex=1),
       col.regions=rev(heat.colors(6)),
       key.space="right",
       par.settings=list(fontsize=list(text=36)))
dev.off()
## Map of CONNECT vs UNCONNECT streams
png("~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/CONNECT.png",
    width=1024, height=768)
spplot(all_sites.sp["UNCONNECT"], bty="n",
       sp.layout=list(l3, l4, l4b, l5, elev),
       main=list(label="Connected and unconnected streams", cex=1),
       col.regions=heat.colors(6),
       key.space="right",
       par.settings=list(fontsize=list(text=36)))
dev.off()
rm(all_sites.sp)



### FIGURES ----------------------------------------------
setwd("~/testRenvProject//Output_Files/2_OrganizeNewData/")

# Histogram of day-of-year when sampling took place:
hist(SampOccasions$yday, xlab="Day-of-year") #  xaxt="n", 
png("EDA/DayOfYearHist.png", width=1024, height=768)
par(mar=c(6, 7, 4, 4))
hist(SampOccasions$yday, xlab="", xaxt="n", las=1, 
     main="Day-of-year that sampling occurred", ylab="",
     cex.lab=3, cex.main=3, cex.axis=3, col="gray")
axis(1, at=seq(100, 350, by=50), cex.axis=3, line=0,
     labels=c("Apr-10", "May-30", "Jul=19",
              "Sep=7", "Oct-27", "Dec-16"))
mtext("No. of sampling occasions", 2, line=4.5, cex=3)
mtext("Day-of-year", 1, line=4, cex=3)
dev.off()
#






# Histogram of nPasses for each survey:
dim(subset(SampOccasions, nSurveys==1))  # 84 with only one pass
plot(as.factor(SampOccasions$nSurveys), xlab="Number of passes") #  xaxt="n", 
png("EDA/nPassesHist.png", width=1024, height=768)
par(mar=c(5, 9, 6, 4))
plot(as.factor(SampOccasions$nSurveys), xlab="", las=1, 
     main="Number of passes per sampling occasion", ylab="", 
     cex.lab=3, cex.main=3, cex.axis=3, cex=3)
dev.off()

# Plot of method used:
plot(as.factor(surveys$METHOD), xlab="Method") #  xaxt="n", 
# dim(subset(surveys, METHOD=="D"))  # 15
# dim(subset(surveys, METHOD=="E"))  # 280
# dim(subset(surveys, METHOD=="S"))  # 125
png("EDA/MethodPlot.png", width=1024, height=768)
par(mar=c(5, 9, 6, 4))
#hist(sites$yday, xlab="Day-of-year") #  xaxt="n", 
plot(as.factor(surveys$METHOD), xlab="", las=1, 
     main="METHOD USED", ylab="", 
     cex.lab=3, cex.main=3, cex.axis=3, cex=3)
mtext("No. of sampling occasions", 2, line=6, cex=3)
dev.off()


### NSPP and MAXCT info ------------------------------
# names(samp_dat)
nspp_df <- unique(subset(samp_dat, 
                         select=-c(sameStatus, Catch, METHOD, PASSNO)))
# nrow=1125
# length(unique(as.character(nspp_df$mySurveyID))) # 201.  good.  matches SampOccasions
# How many sampling occasions without any fish?
# nrow(subset(nspp_df, SpeciesCode=="---"))  # 1...
nspp.table <- as.data.frame(table(nspp_df$mySurveyID))
# names(nspp.table) <- c("mySurveyID", "NSPP")
# dim(nspp.table) # 379. Has extra rows b/c "mySurveyID" was a factor.
# Merge with SampOccasions to only include included occasions:
SampOccasions_wNSPP <- unique( merge(SampOccasions, nspp.table) )
# dim(SampOccasions_wNSPP)
# summary(SampOccasions_wNSPP)
# nrow(subset(SampOccasions_wNSPP, NSPP==1))  # 32.
# nrow=201.  good.  matches SampOccasions
png("EDA/NSPPhist.png", width=1024, height=768)
par(mar=c(6, 7, 4, 4))
##SG: this column NSPP doesn't exist, we don't know why
# hist(SampOccasions_wNSPP$NSPP, xlab="", xaxt="n", las=1, breaks=16,
#      main="Number of species detected per sampling occasion", ylab="",
#      cex.lab=3, cex.main=3, cex.axis=3, col="gray")
# axis(1, at=seq(0, 15, by=5), cex.axis=3, line=0)
# mtext("No. of sampling occasions", 2, line=4.5, cex=3)
#mtext("No. of fish detected", 1, line=4, cex=3)
dev.off()
#
#





## Plot NSPP on map to look for patterns:
# NOT PRINTED- not informative and not pretty
# Doesn't work
# bubble(SampOccasions_wNSPP, "NSPP", maxsize=2.5,  fill=FALSE, add=T)
#

## And add some stats on counts
tmp <-   subset(samp_dat, Catch == "---")
sampNULL <- subset(samp_dat, Catch == "---")$mySurveyID
# samp_dat[samp_dat$mySurveyID %in% sampNULL, c(1:3, 7:10)]
tmp <-   subset(samp_dat, Catch != "---")
maxct <- unique(aggregate(as.numeric(as.character(tmp$Catch)), by=list(as.factor(tmp$mySurveyID)), sum))
# head(maxct)
names(maxct) <- c("mySurveyID", "TOTALCT")
# summary(maxct)
# length(unique(SampOccasions$mySurveyID))  # 201
sites_wCT <- merge(SampOccasions, maxct) # 201-1=200.

png("EDA/TOTALCThist.png", width=1024, height=768)
hist(sites_wCT$TOTALCT, xlab="", ylab="", las=1, xaxt="n", 
     main="Number of fish counted on each sampling occasion", breaks=20,
     cex.lab=3, cex.main=3, cex.axis=3, col="gray")
axis(1, at=seq(0, 3500, by=500), cex.axis=3, line=0)
mtext("No. of sampling occasions", 2, line=4.5, cex=3)
dev.off()





## Plot TOTALCT on map to look for patterns:
# NOT PRINTED- not informative and not pretty
coordinates(sites_wCT) <- c("UTMX", "UTMY")
proj4string(sites_wCT) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#pdf("/Users/broms/Dropbox/Fish_PostDoc/ArkansasR/write_up/NSPP.pdf")
bubble(sites_wCT, "TOTALCT", maxsize=2.5,  fill=FALSE, na.rm=T)
#dev.off()
### nDetects per fish species -----------------
# names(samp_dat)
nDetects <- unique(subset(samp_dat, 
                          select=-c(METHOD, PASSNO, 
                                    sameStatus, Catch)))
# table(nDetects$SpeciesCode)
# dim(unique(cbind(nDetects$UTMX, nDetects$UTMY)))  #  141
#
# And look at nDetects per *site*
#SG: columns watername and WaterNameDateX Not found
nDetects <- unique(subset(samp_dat, 
                          select=-c(METHOD, PASSNO, sameStatus, Catch,
                                    mySurveyID, SampleDate, year, YEAR,
                                    month, yday, YDAY, YDAY2))) #, WaterName, WaterNameDateX
# table(nDetects$SpeciesCode)
# dim(unique(cbind(nDetects$UTMX, nDetects$UTMY)))  #  143
#






### nPasses per samp occasion ------------------
#
orig_dat <- merge(snapped_sites, all_sites)
samp_dat$WaterName <-  orig_dat$WaterName[match(samp_dat$mySurveyID, orig_dat$mySurveyID)]
samp_dat$X <-  as.numeric(with(samp_dat, gsub("\\_.*", "", mySurveyID) ))
samp_dat$WaterNameDateX <-  with(samp_dat, paste(WaterName, "_", SampleDate, 
                                                 "!", UTMX, sep=""))
#
nPasses <- as.matrix(with(samp_dat, table(list(WaterNameDateX, sameStatus))))
# head(nPasses)
WaterNameDateX <- row.names(nPasses)
WaterName <- gsub("\\_.*", "", WaterNameDateX)
# head(WaterName)
UTMX <- gsub(".*\\!", "", WaterNameDateX)
# head(UTMX)
Date.tmp <- gsub(".*\\_", "", WaterNameDateX)
# head(Date.tmp)
Date <- gsub("\\!.*", "", Date.tmp)
# head(Date)
# head(cbind(WaterName, Date, UTMX))
row.names(nPasses) <- NULL
names(nPasses) <- NULL
nPasses.df <- as.data.frame(nPasses[1:nrow(nPasses), 1:8])
# head(nPasses.df)
nPasses.out <- data.frame(WaterName, Date, UTMX, nPasses.df)
# head(nPasses.out)
# dim(nPasses.out)  # 197...  Hmmm- DOES NOT matches sampling occasions.

## Add original surveyID back in for easier clean-up and comparison
# sort(names(orig_dat))
# names(nPasses.out)
# nPasses.out$UTMX <- as.numeric(as.character(nPasses.out$UTMX))
# tmp <- unique(orig_dat[, c("WaterName", "SampleDate", "UTMX", "UTMY", "SurveyID")])
# nPass.wSurveyID <- unique(merge(nPasses.out, tmp,
#                                 by.x=c("WaterName", "Date", "UTMX"),
#                                 by.y=c("WaterName", "SampleDate", "UTMX"), all.x=T))
# head(nPass.wSurveyID)
#

# getwd()
# setwd("/Volumes/CPW_Work/Optimum Sampling/Ark_Sampling_bwa_test/EDA")
write.csv(nPasses.out, "~/testRenvProject//Output_Files/2_OrganizeNewData/EDA/nPassesPerSampOcc.csv", row.names = F)
#

# end of file. -----------------------------
