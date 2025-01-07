#### SelectBalancedSites.R

###   This script imports data and selects spatially-balanced sampling locations
###  10 sites per year, for 3 years: 2016, 2017, 2018

## INPUTS:
## ArkData.RData


##OUTPUTS:
##ArkBalanceSites.RData

## Set directory:
setwd(file.path(here(), "Output_Files/7_SelectBalancedSites"))  # Fill in as appropriate
# setwd("~/ArkModels")  # for server


#remotes::install_version("spsurvey", version = "3.3", repos = "http://cran.us.r-project.org")
#SG: this is just bc i downloaded the package locally
#still need dependencies below though
#install.packages("C:/Users/graffs.NATURENET.001/Downloads/spsurvey_3.3.tar.gz", repos = NULL, type = "source")

#detach("package:foreign", unload=TRUE)

## LIBRARIES:
# library(fields)  # for cover.designs() function

# NOTE: AS of 26Sept2024 spsurvey version is 5.5.1. 
# spsurvey changed functions in 5.0.
# Need older package.  ::::  Try older version 4.1.2



# install old version of spsurvey
# old dependencies for spsurvey
#install.packages("Hmisc")
#"https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_3.17-4.tar.gz"
# geos_loc <- "/Volumes/CPW_Work/Old_R_packages/rgeos_0.6-4.tar.gz"
# install.packages("/Volumes/CPW_Work/Old_R_packages/rgeos_0.6-4.tar.gz", repos=NULL, type="source")
# # SG spot: 
#install.packages("C:/Users/graffs.NATURENET.001/Downloads/rgeos_0.6-4.tar.gz", repos = NULL, type = "source")
#on R version 3.5, might need older version of rtools as well in order to compile packages
#https://cran.r-project.org/bin/windows/Rtools/history.html

# spsurvey_loc <- "https://cran.r-project.org/src/contrib/Archive/spsurvey/spsurvey_4.1.2.tar.gz"
# install.packages(spsurvey_loc, repos=NULL, type="source")
#SG: for Sp dependency:
# install.packages("C:/Users/graffs.NATURENET.001/Downloads/sp_1.2-3.tar.gz", repos = NULL, type = "source")

library(spsurvey)
library(foreign) # read .dbf files
library(sp)
library(sf)



### Import the data and source functions -----------------------
#SG: not loading .rdata bc it doesn't work as well in Rmarkdown. if needed .rdata files are instead called explicitly in Rmarkdown file
#base::load("~/testRenvProject/Output_Files/2_OrganizeNewData/ArkData.Rdata")
# source("~/testRenvProject/Source_Files/ArkFunctions.R")


## Use GRTS because of Basin characterics ------------
# B/c of different structure of Ark River basin, need GRTS to get
# enough sampling on the main stem of the Ark
# Use Unstratified, unequal prob's with panel design (Last Section)

# NOTE: This is what is being used
### Unstratified design, unequal prob's among stream sizes 
### With Panel structure because survey is over time 

## Convert predSiteIDs object to a shapefile


predSites.sp <- pred_sites
coordinates(predSites.sp) <- c("UTMX", "UTMY")
proj4string(predSites.sp) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# this doesn't work use sf
# sp2shape(sp.obj=predSites.sp, shpfilename="predSites_shp") # function isn't available int he spsurvey package anymore.  Need to us sf

# make sf object from the above data
sf_obj <- st_as_sf(predSites.sp)

# write the shapefile
st_write(sf_obj, "predSites_shp.shp")

# Read the shapefile for att in the grts function
att <- foreign::read.dbf("predSites_shp.dbf")
# study_area <- st_read("predSites_shp.dbf")

set.seed(2016)

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

 
# summary(Panelsites)
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n"))
points(Panelsites@coords[, 1], Panelsites@coords[, 2], 
       col=as.factor(Panelsites$panel), pch=16, cex=1.5)
BalancedSites <- subset(Panelsites, panel != "OverSamp")
BalancedSites <- Panelsites
#saveRDS(BalancedSites, file="ArkBalancedSites.rds")
#this is saved in master RUnscript
#save(BalancedSites, file=file.path(here(), "Output_Files/7_SelectBalancedSites/ArkBalancedSites.RData"))
#BalancedSites <- readRDS("Source_Files/ArkBalancedSites.rds")
# summary(BalancedSites@data$pointid)
# summary(AllSites$pointid)
# str(BalancedSites)
#

# end of file. ---------------------------------------------------------
