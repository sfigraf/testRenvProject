##### bring in site/survey data from bios
# this is what gets saved in script 1
data <- read.csv("Input_Files/Cleaned_Data_include_dat_27Sept2024.csv", stringsAsFactors=F)

#bring in Stream network: right now used to get correct coordinate system and check to make sure data is being snapped correctly
arkStreamNetwork <- read_sf("ArcGIS_files/arkansasStreamNetwork.shp")
#remove z element, not needed on 2d maps
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

###snappedSitesMapping
snapped_sites1 <- st_transform(st_as_sf(snapped_sites, coords = c("UTMX", "UTMY"), crs = st_crs(arkStreamNetwork), remove = FALSE), 
                               latLongCRS)
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
                      "Survey Data Sites <br>",
                      "WaterBody:", sitesAllSitesPOints1$WaterName, "<br>", 
                      "UTMX:", sitesAllSitesPOints1$UTMX.y, "<br>", 
                      "UTMY:", sitesAllSitesPOints1$UTMY.y, "<br>"),
                    clusterOptions = markerClusterOptions()) %>%
  addAwesomeMarkers(data = snapped_sites1, 
                    popup = paste(
                      "SNAPPED SITES <br>",
                      "WaterBody:", snapped_sites1$WaterName, "<br>", 
                      "UTMX:", snapped_sites1$UTMX, "<br>", 
                      "UTMY:", snapped_sites1$UTMY, "<br>")
                    #clusterOptions = markerClusterOptions()
                    )
