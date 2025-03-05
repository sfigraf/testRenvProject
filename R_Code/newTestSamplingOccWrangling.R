### My work to get Grant's data to work
library(tidyverse)
library(leaflet)
library(sf)
library(here)
setwd(here())
#bring in data
TestSampingOccasions_Snapped_GRANT <- read_csv("Input_Files/TestSampingOccasions_Snapped_GRANT.csv")

###Deselecting the UTMX and UTMY coords and replacing them with the new correct Snapped Coordindates
TestSampingOccasions_Snapped_GRANT1 <- TestSampingOccasions_Snapped_GRANT %>%
  select(-UTMX, -UTMY) %>%
  rename(UTMX = SNAPUTM_X, 
         UTMY= SNAPUTM_Y)


# ##subject new Grant data to same processes as old snapped Sites ---------------

TestSampingOccasions_Snapped_GRANT1[TestSampingOccasions_Snapped_GRANT1$SurveyID==43945, "UTMX"] <- 620287
TestSampingOccasions_Snapped_GRANT1[TestSampingOccasions_Snapped_GRANT1$SurveyID==43945, "UTMY"] <- 4167592
TestSampingOccasions_Snapped_GRANT1[TestSampingOccasions_Snapped_GRANT1$SurveyID==45532, "UTMX"] <- 620287
TestSampingOccasions_Snapped_GRANT1[TestSampingOccasions_Snapped_GRANT1$SurveyID==45532, "UTMY"] <- 4167592
#
#new column for SpeciesCod is 21, so it's just easiest to find that column by name
names(TestSampingOccasions_Snapped_GRANT1)[which(names(TestSampingOccasions_Snapped_GRANT1) == "SpeciesCod")] <- "SpeciesCode"


TestSampingOccasions_Snapped_GRANT1 <- subset(TestSampingOccasions_Snapped_GRANT1, !(SurveyID %in% c(25784, 23897))) 

# bring in allSites data, same as script 2 --------------------------------
## Read in ArcGIS Stream Network info 
streams <- foreign::read.dbf("Input_Files/SurveyDesignStreams.dbf")
# str(streams)
## Change INT/PER to more proper Connect/Unconnect
streams$StreamType <- ifelse(streams$StreamType=="INT", "UNCONNECT", "CONNECT")

stream_pts <- foreign::read.dbf("Input_Files/SurveyDesignStreamsAsPoints.dbf")
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


# MERGE -------------------------------------------------------------------

# base::merge() finds columns in common by running this line: intersect(names(TestSampingOccasions_Snapped_GRANT1), names(all_sites))
# in the new Grant Data, there is an ObjectID column from GIS that shouldn't be there; specify which columns to merge on
# or better yet, make this data have the same columns as snapped_sites (some columns might need to be renamed)
samp_datGrant <- base::merge(TestSampingOccasions_Snapped_GRANT1, all_sites, by = c("UTMX","UTMY"))


# QAQC --------------------------------------------------------------------
###check to make sure grant Data lines up with snapped Sites and stream_points
###visually check by mapping, using leaflet

#mapping prep
#stream network
arkStreamNetwork <- read_sf("ArcGIS_files/AsShapefiles/arkansasStreamNetwork.shp")
#remove z element, not needed on 2d maps
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")
#gets utms as columns
arkStreamNetwork1 <- arkStreamNetwork1 %>%
  mutate(UTMX = st_coordinates(arkStreamNetwork1)[row_number(),1], 
         UTMY = st_coordinates(arkStreamNetwork1)[row_number(),2])
#define latLongcrs that's necessary for mapping with leaflet
latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326
#transform the stream Network to new crs
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) 

#make spatial object of UTMs before they were snapped
#has slightly more rows than TestSampingOccasions_Snapped_GRANT1 bc snapped_sites wrangling has not happened on this orginal file
UnSnapped_GRANTSF <- st_transform(st_as_sf(TestSampingOccasions_Snapped_GRANT, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS) 

#make spatial object of UTMs after they were snapped
#has slightly more rows than TestSampingOccasions_Snapped_GRANT1 bc snapped_sites wrangling has not happened on this orginal file
Snapped_GRANTSF <- st_transform(st_as_sf(TestSampingOccasions_Snapped_GRANT, coords = c("SNAPUTM_X", "SNAPUTM_Y"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                                                     latLongCRS)
#create map using leaflet package
#stream network is first data defined
# to save map, click "export" in the viewer and "export as webpage"
leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines(
    popup = paste(
      "WaterBody:", arkStreamNetwork2$NAME1, "<br>"
    )
  ) %>%

  addAwesomeMarkers(data = Snapped_GRANTSF, 
                    group = "New Sites Grant Snapped", 
                    popup = paste(
                      "Grant Snapped Site <br>",
                      "WaterBody:", Snapped_GRANTSF$WaterName, "<br>", 
                      "UTMX:", Snapped_GRANTSF$SNAPUTM_X, "<br>", 
                      "UTMY:", Snapped_GRANTSF$SNAPUTM_Y, "<br>"),
                    clusterOptions = markerClusterOptions()
  ) %>%
  addAwesomeMarkers(data = st_transform(st_as_sf(stream_pts, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), latLongCRS),
                    group = "Stream Points", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      iconColor = 'black',
                      library = 'ion',
                      #iconHeight = 20,
                      markerColor = "green"
                    ), 
                    popup = paste(
                      "Stream Points <br>",
                      "UTMX:", stream_pts$UTMX, "<br>", 
                      "UTMY:", stream_pts$UTMY, "<br>"),
                    clusterOptions = markerClusterOptions()
  ) %>%
  addAwesomeMarkers(data = st_transform(UnSnapped_GRANTSF, latLongCRS),
                    group = "New Sites Grant UNSnapped",
                    icon = leaflet::awesomeIcons(
                      library = 'ion',
                      markerColor = "red"
                    ),
                    popup = paste(
                      "Grant UNSnapped Site <br>",
                      "WaterBody:", UnSnapped_GRANTSF$WaterName, "<br>",
                      "UTMX:", UnSnapped_GRANTSF$UTMX, "<br>",
                      "UTMY:", UnSnapped_GRANTSF$UTMY, "<br>"),
                    clusterOptions = markerClusterOptions()) %>%
  addLayersControl(overlayGroups = c( "New Sites Grant Snapped", "New Sites Grant UNSnapped", "Stream Points")) %>%
  hideGroup(c("Stream Points", "New Sites Grant UNSnapped"))


