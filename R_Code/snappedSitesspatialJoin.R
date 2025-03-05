####SCRIPT TO "SNAP" DATA to Stream Network
library(tidyverse)
library(sf)
library(leaflet)
library(here)

##### bring in site/survey data from Andy/Brian
# this is what gets saved in script 1
setwd(here())
data <- read.csv("Input_Files/Cleaned_Data_include_dat_27Sept2024.csv", stringsAsFactors=F)

#if you want to append more data to the file, you can save it as a csv here 
# then read it in again in the above line with the new name changed. make sure columns stay the same
write.csv(data, "newData.csv")

#bring in Stream network
# right now used to get correct coordinate system and check to make sure data is being snapped correctly
arkStreamNetwork <- read_sf("ArcGIS_files/AsShapefiles/arkansasStreamNetwork.shp")
#remove z element, not needed on maps and causes a little errors 
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")

#change data to sf object in preparation for spatial join with same crs as streamNetwork
surveySitesSF <- st_as_sf(data, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE)


### bring in stream points data: a point file made from stream network to reliably snap sites
stream_pts <- foreign::read.dbf("Input_Files/SurveyDesignStreamsAsPoints.dbf")
#convert to sf object in order to spatially join with survey data, using the same CRS as the stream network
stream_ptsNADSF <- st_as_sf(stream_pts, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE)

#spatially join DFs to nearest feature. Currently keeps geometry of left dataset but contains UTMs as columns of both DFs
newSitesSnapped0 <- sf::st_join(surveySitesSF, stream_ptsNADSF, st_nearest_feature)
#make new df of data with desired UTM coords from stream_points
newSitesSnapped1 <- newSitesSnapped0 %>%
  #first change back to normal DF
  as.data.frame() %>%
  #then change back to sf object with desired coordinates and desired CRS
  st_as_sf(coords = c("UTMX.y", "UTMY.y"), crs = st_crs(arkStreamNetwork), remove = FALSE)


#give hte new data the same columns as the original snaped_surveys
#read in surveys
snapped_sites <- foreign::read.dbf("Input_Files/snappedSurveys.dbf")

newSitesSnapped2 <- as.data.frame(newSitesSnapped1) %>%
  rename(UTMX = UTMX.y, 
         UTMY = UTMY.y)
#these names are column names of snapped_sites, with a few of them changed for how snappedSurveys came out of GIS (ie SurveyPurp to SurveyPurpose)
#easiest to find with this line: cat(names(snapped_sites), sep = "', '")
newNames <- c('SurveyPurpose', 'Protocol', 'Gear', 'TotalEffort', 'EffortMetric', 'WaterID', 'WaterName', 'SampleDate', 'Location', 
              'Elevation', 'UTMZone', 'UTMX', 'UTMY', 'Basin', 'StationLength', 'AvgWidth', 'StationCode', 'SiteType', 'SurveyID', 
              'SpeciesCode', 'Status', 'NewStatus', 'sameStatus', 'Catch', 'Source', 'mySurveyID')

##THIS is the new equivalent to snapped_sites
#plug in newSnappedSites back to script 2 as if it was snapped_sites that just got read into the script as a dbf file
newSnappedSites <- newSitesSnapped2 %>%
  select(all_of(newNames))



# QAQC  -------------------------------------------------------------------



###check how it lines up with snapped Sites and stream_points
###visually check by mapping, using leaflet

#define latLongcrs that's necessary for mapping with leaflet
latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326
#transform the stream Network to new crs
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) 

##as Joined
#create "differences" file between newSitesSnapped and previous snapped sites
#snapped sites has 300 less entries from stuff done by original author
differences <- anti_join(as.data.frame(newSitesSnapped1), snapped_sites, by = c("UTMX.y" = "UTMX", "UTMY.y" = "UTMY"))
#make differences a sf object for mapping later
differencesSF <- st_as_sf(differences)

#change snapped sites to latLong crs
#also take away duplicate entries, because we just need to know if hte new snapped sites are located in the same spot as the old ones,
#not how many surveys there are at a specific location
snapped_sites1 <- st_transform(st_as_sf(snapped_sites, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS) %>%
  distinct(UTMX, UTMY, .keep_all = T)

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
  #new data is added with different colors, groups, and custom labels/popups
  #using newSitesSnapped1 because it's still an sf object
  #you can plot non-sf objects but you just need to specify lat/long coords
  addAwesomeMarkers(data = st_transform(newSitesSnapped1, latLongCRS), 
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
                      "WaterBody:", newSitesSnapped1$WaterName, "<br>", 
                      "UTMX:", newSitesSnapped1$UTMX.y, "<br>", 
                      "UTMY:", newSitesSnapped1$UTMY.y, "<br>"),
                    clusterOptions = markerClusterOptions()) %>%
  addAwesomeMarkers(data = snapped_sites1, 
                    group = "Snapped Sites", 
                    popup = paste(
                      "SNAPPED SITES <br>",
                      "WaterBody:", snapped_sites1$WaterName, "<br>", 
                      "UTMX:", snapped_sites1$UTMX, "<br>", 
                      "UTMY:", snapped_sites1$UTMY, "<br>")
                    #can make this a cluster as well but not really needed since there's just 1 entry per UTM combo for this data
                    #clusterOptions = markerClusterOptions()
  ) %>%
  addAwesomeMarkers(data = st_transform(differencesSF, latLongCRS),
                    group = "Differences", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      iconColor = 'black',
                      library = 'ion',
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
  #makes layer control for user
  addLayersControl(overlayGroups = c("Survey Sites (Post-Snap)", "Survey Sites (Pre-Snap)", "Snapped Sites", "Differences", "Stream Points")) %>%
  #start with some layers off
  hideGroup(c("Differences", "Stream Points", "Survey Sites (Pre-Snap)"))

