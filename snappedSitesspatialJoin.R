####SCRIPT TO "SNAP" DATA to Stream Network
library(tidyverse)
library(sf)
library(leaflet)
##### bring in site/survey data from Andy/Brian
# this is what gets saved in script 1
data <- read.csv("Input_Files/Cleaned_Data_include_dat_27Sept2024.csv", stringsAsFactors=F)

#bring in Stream network: right now used to get correct coordinate system and check to make sure data is being snapped correctly
arkStreamNetwork <- read_sf("ArcGIS_files/arkansasStreamNetwork.shp")
#remove z element, not needed on maps and causes a little errors 
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")
#gets utms as columns
arkStreamNetwork1 <- arkStreamNetwork1 %>%
  mutate(UTMX = st_coordinates(arkStreamNetwork1)[row_number(),1], 
         UTMY = st_coordinates(arkStreamNetwork1)[row_number(),2])
#change data to sf object in preparation for spatial join with same crs as streamNetwork
surveySitesSF <- st_as_sf(data, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE)


### get sites as point data, aka all sites, aka joined streams and streampoints file
#could maybe just do stream points
stream_pts <- foreign::read.dbf("Input_Files/SurveyDesignStreamsAsPoints.dbf")
#convert to sf object
stream_ptsNADSF <- st_as_sf(stream_pts, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE)

#spatially join to nearest feature. Currently keeps geometry of left dataset but contains UTMs as columns of both DFs
sitesAllSitesPOints <- sf::st_join(surveySitesSF, stream_ptsNADSF, st_nearest_feature)
#make new df of data with desired UTM coords from streamPOints df
sitesAllSitesPOints1 <- sitesAllSitesPOints %>%
  as.data.frame() %>%
  st_as_sf(coords = c("UTMX.y", "UTMY.y"), crs = st_crs(arkStreamNetwork), remove = FALSE)

###check how it lines up with snapped Sites 
snapped_sites <- foreign::read.dbf("Input_Files/snappedSurveys.dbf")
#snapped sites has 300 less entries, idk why, q for brian/ryan

write.csv(as.data.frame(sitesAllSitesPOints1), "sitesAllSitesPOints1.csv")
write.csv(as.data.frame(snapped_sites), "snapped_sites.csv")

###visually check

latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) %>%
  mutate(Lat = st_coordinates(arkStreamNetwork2)[row_number(),1], 
         Long = st_coordinates(arkStreamNetwork2)[row_number(),2])

#snapped sites ready to map

##as Joined
#these are the differences in UTMs between OG data and new data
differences <- anti_join(as.data.frame(sitesAllSitesPOints1), snapped_sites, by = c("UTMX.y" = "UTMX", "UTMY.y" = "UTMY"))
differencesSF <- st_as_sf(differences)


###snappedSitesMapping
snapped_sites1 <- st_transform(st_as_sf(snapped_sites, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS) %>%
  distinct(UTMX, UTMY, .keep_all = T)

leaflet(arkStreamNetwork2) %>%
  addTiles() %>%
  addPolylines(
    popup = paste(
      "WaterBody:", arkStreamNetwork2$NAME1, "<br>"
    )
  ) %>%
  addAwesomeMarkers(data = st_transform(sitesAllSitesPOints1, latLongCRS), 
                    group = "Survey Sites (Post-Snap)",
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      iconColor = 'black',
                      library = 'ion',
                      #iconHeight = 20,
                      markerColor = "purple"
                    ), 
                    popup = paste(
                      "Survey Data Sites (Post-Snap) <br>",
                      "WaterBody:", sitesAllSitesPOints1$WaterName, "<br>", 
                      "UTMX:", sitesAllSitesPOints1$UTMX.y, "<br>", 
                      "UTMY:", sitesAllSitesPOints1$UTMY.y, "<br>"),
                    clusterOptions = markerClusterOptions()) %>%
  addAwesomeMarkers(data = snapped_sites1, 
                    group = "Snapped Sites", 
                    popup = paste(
                      "SNAPPED SITES <br>",
                      "WaterBody:", snapped_sites1$WaterName, "<br>", 
                      "UTMX:", snapped_sites1$UTMX, "<br>", 
                      "UTMY:", snapped_sites1$UTMY, "<br>")
                    #clusterOptions = markerClusterOptions()
  ) %>%
  addAwesomeMarkers(data = st_transform(differencesSF, latLongCRS),
                    group = "Differences", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      iconColor = 'black',
                      library = 'ion',
                      #iconHeight = 20,
                      markerColor = "orange"
                    ), 
                    popup = paste(
                      "Survey Data Sites: NOT in Snapped Sites <br>",
                      "WaterBody:", differencesSF$WaterName, "<br>", 
                      "UTMX:", differencesSF$UTMX.y, "<br>", 
                      "UTMY:", differencesSF$UTMY.y, "<br>"),
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
  addAwesomeMarkers(data = st_transform(surveySitesSF, latLongCRS), 
                    group = "Survey Sites (Pre-Snap)", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      iconColor = 'black',
                      library = 'ion',
                      #iconHeight = 20,
                      markerColor = "yellow"
                    ), 
                    popup = paste(
                      "Survey Data Sites (Pre-Snap) <br>",
                      "WaterBody:", surveySitesSF$WaterName, "<br>", 
                      "UTMX:", surveySitesSF$UTMX, "<br>", 
                      "UTMY:", surveySitesSF$UTMY, "<br>"),
                    clusterOptions = markerClusterOptions()) %>%
  addLayersControl(overlayGroups = c("Survey Sites (Post-Snap)", "Survey Sites (Pre-Snap)", "Snapped Sites", "Differences", "Stream Points")) %>%
  hideGroup(c("Differences", "Stream Points", "Survey Sites (Pre-Snap)"))
  


###making new df similar to snapped sites
# cat(names(snapped_sites), sep = "', '")
sitesAllSitesPOints2 <- as.data.frame(sitesAllSitesPOints1) %>%
  rename(UTMX = UTMX.y, 
         UTMY = UTMY.y)
newNames <- c('SurveyPurpose', 'Protocol', 'Gear', 'TotalEffort', 'EffortMetric', 'WaterID', 'WaterName', 'SampleDate', 'Location', 
              'Elevation', 'UTMZone', 'UTMX', 'UTMY', 'Basin', 'StationLength', 'AvgWidth', 'StationCode', 'SiteType', 'SurveyID', 
              'SpeciesCode', 'Status', 'NewStatus', 'sameStatus', 'Catch', 'Source', 'mySurveyID')

newSnappedSites <- sitesAllSitesPOints2 %>%
  select(all_of(newNames))

#####
#this is allsites after it has been merged and wrnalged with streams and stream points in script 2

allSitesWrangled <- all_sites
#original 2016 snapped surveys
#columns 20
#

#this stuff is all the snapped sites data wrangling
## After comparing data to stream network, need to edit one site's location:
newSnappedSites[newSnappedSites$SurveyID==43945, "UTMX"] <- 620287
newSnappedSites[newSnappedSites$SurveyID==43945, "UTMY"] <- 4167592
newSnappedSites[newSnappedSites$SurveyID==45532, "UTMX"] <- 620287
newSnappedSites[newSnappedSites$SurveyID==45532, "UTMY"] <- 4167592
#
#names(newSnappedSites)
# survey ID changed to Species Code
#SG: 
#
names(newSnappedSites)[which(names(newSnappedSites) == "SpeciesCod")] <- "SpeciesCode"

## Remove some sites that are still hanging around and shouldn't be:
## Sites are repeated counts of other surveys, but slightly diff coords:
#SG: changedSurveyID here to Species Code, but also code 25784 doesn't seem to even be in here
newSnappedSites <- subset(newSnappedSites, !(SurveyID %in% c(25784, 23897)))
# dim(newSnappedSites)  # 3197 obs.
#

## IF YOU GET TO END AND NEED TO REMOVE MORE SAMP OCCASIONS -------------
# Run the "Ark_OrgSampOccForExport.R" on the new/updated data spreadsheet.
# tmp1 <- include_dat[, c("mySurveyID", "Catch", 
#                         "SpeciesCodee", "sameStatus")]
# tmp2 <- newSnappedSites[, c("WaterName", "SampleDate", "UTMX", "UTMY",
#                           "Source", "mySurveyID")]
# tmp2$mySurveyID <- as.character(tmp2$mySurveyID)
# tmp <- merge(tmp1, tmp2, all=T)
# tmp <- unique(tmp)
# tmp$mySurveyID <- as.factor(tmp$mySurveyID)
# summary(tmp)
# dim(tmp)
# dim(include_dat)  # should match, then:
# newSnappedSites <- tmp
# -----------------------------------

# Add covariate info associated with each sampled site: 
# names(newSnappedSites)
# names(all_sites)

##with new joining, seems to work!!
#these are the columns that automatically are merged
#intersect(names(all_sites), names(newSnappedSites))

samp_datNew <- base::merge(all_sites, newSnappedSites)
samp_datNew <- merge(newSnappedSites, all_sites)


##checking for difs
difs <- anti_join(samp_datNew, samp_dat1, join_by(UTMX, UTMY))

sampdatNew2 <- anti_join(samp_datNew, difs)

samp_dat1 <- samp_dat %>%
  rename(StationCode = StationCod, 
         StationLength = StationLen, 
         EffortMetric = EffortMetr, 
         TotalEffort = TotalEffor, 
         SurveyPurpose = SurveyPurp
         )
difs <- left_join(samp_datNew, samp_dat1, join_by(UTMX, UTMY, mySurveyID)) #, SurveyPurpose, Protocol, Gear, TotalEffort, EffortMetric
difs <- difs %>%
  select(order(colnames(.))) %>%
  arrange(mySurveyID)

names(difs)
