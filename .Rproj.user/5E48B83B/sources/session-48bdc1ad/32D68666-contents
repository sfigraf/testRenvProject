

datasub <- data[,c("WaterID", "WaterName", "Catch", "mySurveyID")]

snapped_sites <- foreign::read.dbf("Input_Files/snappedSurveys.dbf")
snapped_sitessub <- snapped_sites[,c("WaterID", "WaterName", "Catch", "mySurveyID", "UTMX", "UTMY")]

difs <- left_join(snapped_sites, data,  by = c("Protocol", "Gear", "SampleDate", "Location", "Elevation", "UTMZone", "UTMX", "UTMY", "Basin", "AvgWidth", "SiteType", 
                                              "SurveyID", "Status", "NewStatus",
                  "sameStatus", "Catch", "Source", "mySurveyID"))

library(sf)
latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326

arkStreamNetwork <- read_sf("ArcGIS_files/arkansasStreamNetwork.shp")
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")
arkStreamNetwork1 <- arkStreamNetwork1 %>%
  mutate(UTMX = st_coordinates(arkStreamNetwork1)[row_number(),1], 
         UTMY = st_coordinates(arkStreamNetwork1)[row_number(),2])
sitesSF <- st_as_sf(data, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE)
#st_crs(sitesSF)
#trasnformedSites <- sf::st_transform(sitesSF, latLongCRS)
#stattions needed to calculate movements and distance moved
sitesAndNetwork <- sf::st_join(sitesSF, arkStreamNetwork1, st_nearest_feature)
###trying snapping to allSites
all_sitesNADSF <- st_as_sf(all_sites, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE)

sitesAllSitesPOints <- sf::st_join(sitesSF, all_sitesNADSF, st_nearest_feature)

leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines() %>%
  addAwesomeMarkers(data = st_transform(sitesSF, latLongCRS), 
             icon = leaflet::awesomeIcons(
               icon = 'add',
               iconColor = 'black',
               library = 'ion',
               #iconHeight = 20,
               markerColor = "purple"
             ), 
             popup = paste(
               "Survey Data Sites <br>",
               "WaterBody:", sitesSF$WaterName, "<br>", 
               "UTMX:", sitesSF$UTMX, "<br>", 
               "UTMY:", sitesSF$UTMY, "<br>"),
             clusterOptions = markerClusterOptions()) %>%
  addAwesomeMarkers(data = st_transform(all_sitesNADSF, latLongCRS), 
             popup = paste(
               "ALL SITES <br>",
               "WaterBody:", all_sitesNADSF$GNIS_Name, "<br>", 
               "UTMX:", all_sitesNADSF$UTMX, "<br>", 
               "UTMY:", all_sitesNADSF$UTMY, "<br>"),
             clusterOptions = markerClusterOptions())

####
#use the new UTM from allSites as point 
sitesAllSitesPOints1 <- sitesAllSitesPOints %>%
  as.data.frame() %>%
  st_as_sf(coords = c("UTMX.y", "UTMY.y"), crs = st_crs(arkStreamNetwork), remove = FALSE)

##check if it got it right?
leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines() %>%
  addAwesomeMarkers(data = st_transform(sitesAllSitesPOints1, latLongCRS), 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      iconColor = 'black',
                      library = 'ion',
                      #iconHeight = 20,
                      markerColor = "purple"
                    ), 
                    popup = paste(
                      "Joined Sites <br>",
                      "WaterBody:", sitesAllSitesPOints1$WaterName, "<br>", 
                      "UTMX:", sitesAllSitesPOints1$UTMX.y, "<br>", 
                      "UTMY:", sitesAllSitesPOints1$UTMY.y, "<br>"),
                    clusterOptions = markerClusterOptions())

x <- sitesAndNetwork %>%
  select(WaterName, NAME1) %>%
  mutate(UTMX = st_coordinates(sitesAndNetwork)[row_number(),1], 
         UTMY = st_coordinates(sitesAndNetwork)[row_number(),2], )

#install.packages("leaflet")
library(leaflet)
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) %>%
  mutate(Lat = st_coordinates(arkStreamNetwork2)[row_number(),1], 
         Long = st_coordinates(arkStreamNetwork2)[row_number(),2])
leaflet(st_transform(sitesSF, latLongCRS)) %>%
  addTiles() %>%
  addMarkers(clusterOptions = markerClusterOptions()) %>%
  addPolygons(data = arkStreamNetwork2, lat = ~arkStreamNetwork2$Lat, lng = ~arkStreamNetwork2$Long)
  #addPolylines(st_transform(arkStreamNetwork1, latLongCRS))

leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines() %>%
  addMarkers(data = st_transform(sitesSF, latLongCRS), 
             popup = paste(
               "WaterBody:", sitesSF$WaterName, "<br>", 
               "UTMX:", sitesSF$UTMX, "<br>", 
               "UTMY:", sitesSF$UTMY, "<br>"),
             clusterOptions = markerClusterOptions())


###snappedSitesMapping
snapped_sites1 <- st_transform(st_as_sf(snapped_sites, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS)

leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines(
    popup = paste(
      "WaterBody:", arkStreamNetwork2$NAME1, "<br>"
    )
  ) %>%
  addMarkers(data = snapped_sites1, 
             popup = paste(
               "WaterBody:", snapped_sites1$WaterName, "<br>", 
               "UTMX:", snapped_sites1$UTMX, "<br>", 
               "UTMY:", snapped_sites1$UTMY, "<br>"
               ), 
             clusterOptions = markerClusterOptions())


#########snapped sites data wrangling
sitesAllSitesPOints2 <- as.data.frame(sitesAllSitesPOints1) %>%
  rename(UTMX = UTMX.y, 
         UTMY = UTMY.y)
newNames <- c('SurveyPurpose', 'Protocol', 'Gear', 'TotalEffort', 'EffortMetric', 'WaterID', 'WaterName', 'SampleDate', 'Location', 
              'Elevation', 'UTMZone', 'UTMX', 'UTMY', 'Basin', 'StationLength', 'AvgWidth', 'StationCode', 'SiteType', 'SurveyID', 
              'SpeciesCode', 'Status', 'NewStatus', 'sameStatus', 'Catch', 'Source', 'mySurveyID')

newSnappedSites <- sitesAllSitesPOints2 %>%
  select(all_of(newNames))
#survey data snapped to correct location in the river (
#one that works is snapped data from july 6 2016
#new snapped surveys that brian added doesn't work
#original 2016 snapped surveys
#columns 20
#
snapped_sites <- foreign::read.dbf("Input_Files/snappedSurveys.dbf")

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

####All sites wrangling

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


##streamNetworkandAllSites 
all_sites1 <- st_transform(st_as_sf(all_sites, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS)
stream_pts1 <- st_transform(st_as_sf(stream_pts, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                           latLongCRS)

leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines(
    popup = paste(
      "WaterBody:", arkStreamNetwork2$NAME1, "<br>"
    )
  ) %>%
  addMarkers(data = all_sites1, 
             popup = paste(
               "WaterBody:", all_sites1$GNIS_Name, "<br>", 
               "UTMX:", all_sites1$UTMX, "<br>", 
               "UTMY:", all_sites1$UTMY, "<br>"
             ), 
             clusterOptions = markerClusterOptions())
  

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

