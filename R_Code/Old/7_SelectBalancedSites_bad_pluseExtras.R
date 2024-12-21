#### SelectBalancedSites.R

###   This script imports data and selects spatially-balanced sampling locations
###  10 sites per year, for 3 years: 2016, 2017, 2018

## INPUTS:
## ArkData.RData


##OUTPUTS:
##ArkBalanceSites.RData

## Set directory:
setwd("/Volumes/CPW_Work/Optimum Sampling/Ark_Optimal_bwa/Output_Files/SelectBalancedSites")  # Fill in as appropriate
# setwd("~/ArkModels")  # for server


# remotes::install_version("spsurvey", version = "3.3", repos = "http://cran.us.r-project.org")





## LIBRARIES:
# library(fields)  # for cover.designs() function
# library(spsurvey)  # for GRTS survey design
# NOTE: spsurvey is different as of version 5.0. Try older version 4.1.4


# install.packages("/Volumes/CPW_Work/Old_R_packages/spsurvey_4.1.4.tar.gz", repos = NULL, type="source")

# install old version of spsurvey
# old dependencies for spsurvey
# install.packages("Hmisc")
# geos_loc <- "/Volumes/CPW_Work/Old_R_packages/rgeos_0.6-4.tar.gz"
# install.packages("/Volumes/CPW_Work/Old_R_packages/rgeos_0.6-4.tar.gz", repos=NULL, type="source")
# # 


# spsurvey_loc <- "https://cran.r-project.org/src/contrib/Archive/spsurvey/spsurvey_4.1.2.tar.gz"
# install.packages(spsurvey_loc, repos=NULL, type="source")


library(spsurvey)

# library(foreign) # read .dbf files
library(sp)
library(sf)



### Import the data and source functions -----------------------
load("/Volumes/CPW_Work/Optimum Sampling/Ark_Optimal_bwa/Output_Files/ArkData.RData")
source("/Volumes/CPW_Work/Optimum Sampling/Ark_Optimal_bwa/Source_Files/ArkFunctions.R")
# load("~/ArkCrossVal/ArkData.RData")  # for server
# source("~/ArkCrossVal/ArkFunctions.R")  # for server

## Use GRTS because of Basin characterics ------------
# B/c of different structure of Ark River basin, need GRTS to get
# enough sampling on the main stem of the Ark
# Use Unstratified, unequal prob's with panel design (Last Section)

### Unstratified design, unequal prob's among stream sizes 
### With Panel structure because survey is over time 

## Convert predSiteIDs object to a shapefile

for(i in 1:length(pred_sites$StreamSize)){
  if(pred_sites$StreamSize[i]==1){
    pred_sites$StreamSize2[i] <- "one"
  }
  if(pred_sites$StreamSize[i]==2){
    pred_sites$StreamSize2[i] <- "two"
  }
  if(pred_sites$StreamSize[i]==3){
    pred_sites$StreamSize2[i] <- "three"
  }
  if(pred_sites$StreamSize[i]==4){
    pred_sites$StreamSize2[i] <- "four"
  }
}



predSites.sp <- pred_sites
coordinates(predSites.sp) <- c("UTMX", "UTMY")
proj4string(predSites.sp) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

sp2shape(sp.obj=predSites.sp, shpfilename="predSites_shp") # function isn't available int he spsurvey package anymore.  Need to us sf

# make sf object from the above data
sf_obj <- st_as_sf(predSites.sp)

# write the shapefile
# st_write(sf_obj, "predSites_shp.shp")

# Read the shapefile
# att <- read.dbf("predSites_shp.dbf")
# study_area <- st_read("predSites_shp.dbf")

# set.seed(2016)

Paneldsgn <- list(None=list(panel=c(Year1=10, Year2=10, Year3=10), # number of points that you want in the strata,
                            seltype="Unequal",
                            caty.n=c("1"=4*3,
                                     "2"=1*3,
                                     "3"=2*3,
                                     "4"=3*3),
                            over=10))
#grts stands for generalized random tesselation stratified sampling
Panelsites <- grts(design=Paneldsgn,
                     DesignID="UNEQUAL",
                     type.frame="finite",
                     src.frame="shapefile",
                    in.shape="predSites_shp.shp",
                    att.frame=att,
                    mdcaty="SIZE",
                    shapefile=FALSE)

# New function from spsurvey - grts
# first argument is the sampling frame (sf object)
# n_base the number of sites in the base(main) sample: total

# unequal inclusion probabilites based on "StreamSize" 1, 2, 3, 4, changed numbers to chr above to get it to work with grts
caty.n=c(one=12, two=3, three=6, four=9)

sites_grts <- grts(sf_obj, n_base=30, 
                   caty_var="StreamSize2",
                   caty_n=caty.n,
                   n_over=10)
sites_grts
plot(sites_grts, sf_obj, key.width = lcm(3))
 
summary(sites_grts)
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n"))
points(Panelsites@coords[, 1], Panelsites@coords[, 2], 
       col=as.factor(Panelsites$panel), pch=16, cex=1.5)
BalancedSites <- subset(Panelsites, panel != "OverSamp")
BalancedSites <- Panelsites
save(BalancedSites, file="ArkBalancedSites.RData")

summary(BalancedSites@data$pointid)
summary(AllSites$pointid)
str(BalancedSites)
#

# end of file. ---------------------------------------------------------