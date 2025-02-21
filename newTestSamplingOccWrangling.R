
Ark_Data_TEST_bwa <- read_csv("Ark_Data_TEST_bwa.csv")

TestSampingOccasions_Snapped_GRANT <- read_csv("TestSampingOccasions_Snapped_GRANT.csv")

#mapping prep
#stream network
arkStreamNetwork <- read_sf("ArcGIS_files/arkansasStreamNetwork.shp")
#remove z element, not needed on 2d maps
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")
#gets utms as columns
arkStreamNetwork1 <- arkStreamNetwork1 %>%
  mutate(UTMX = st_coordinates(arkStreamNetwork1)[row_number(),1], 
         UTMY = st_coordinates(arkStreamNetwork1)[row_number(),2])
latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) 

# stream points
stream_pts <- foreign::read.dbf("Input_Files/SurveyDesignStreamsAsPoints.dbf")


UnSnapped_GRANTSF <- st_transform(st_as_sf(TestSampingOccasions_Snapped_GRANT, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS) #%>%
  #distinct(UTMX, UTMY, .keep_all = T)

Snapped_GRANTSF <- st_transform(st_as_sf(TestSampingOccasions_Snapped_GRANT, coords = c("SNAPUTM_X", "SNAPUTM_Y"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                                                     latLongCRS) #%>%
  #distinct(UTMX, UTMY, .keep_all = T)

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

####Merging with allSites

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
# dim(pred_sites)
# 479 sites to choose from.


### STEP 2: Add SampOccasion-level info and covariates --------

# Import the snapped locations to get correct covariate info:
#survey data snapped to correct location in the river (
#one that works is snapped data from july 6 2016
#new snapped surveys that brian added doesn't work
#original 2016 snapped surveys
#columns 20
#
#snapped_sites <- foreign::read.dbf("Input_Files/snappedSurveys.dbf")
#TestSampingOccasions_Snapped_GRANT
## After comparing data to stream network, need to edit one site's location:
TestSampingOccasions_Snapped_GRANT[TestSampingOccasions_Snapped_GRANT$SurveyID==43945, "SNAPUTM_X"] <- 620287
TestSampingOccasions_Snapped_GRANT[TestSampingOccasions_Snapped_GRANT$SurveyID==43945, "SNAPUTM_Y"] <- 4167592
TestSampingOccasions_Snapped_GRANT[TestSampingOccasions_Snapped_GRANT$SurveyID==45532, "SNAPUTM_X"] <- 620287
TestSampingOccasions_Snapped_GRANT[TestSampingOccasions_Snapped_GRANT$SurveyID==45532, "SNAPUTM_Y"] <- 4167592
#
#names(snapped_sites)
# survey ID changed to Species Code
#SG: 
names(TestSampingOccasions_Snapped_GRANT)[which(names(TestSampingOccasions_Snapped_GRANT) == "SpeciesCod")] <- "SpeciesCode"
#names(TestSampingOccasions_Snapped_GRANT)[20] <- "SpeciesCode"

## Remove some sites that are still hanging around and shouldn't be:
## Sites are repeated counts of other surveys, but slightly diff coords:
#SG: changedSurveyID here to Species Code, but also code 25784 doesn't seem to even be in here
TestSampingOccasions_Snapped_GRANT <- subset(TestSampingOccasions_Snapped_GRANT, !(SurveyID %in% c(25784, 23897))) 

TestSampingOccasions_Snapped_GRANT1 <- TestSampingOccasions_Snapped_GRANT %>%
  select(-UTMX, -UTMY) %>%
  rename(UTMX = SNAPUTM_X, 
         UTMY= SNAPUTM_Y)

#left_join works, merge() doesnt...
samp_datGrant <- base::merge(all_sites, TestSampingOccasions_Snapped_GRANT1, by = c("UTMX","UTMY" ))
samp_datGrant <- merge(TestSampingOccasions_Snapped_GRANT, all_sites)
#autmatically joins on this:
intersect(names(all_sites), names(TestSampingOccasions_Snapped_GRANT))
