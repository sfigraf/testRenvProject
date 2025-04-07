
# plot sites selected as best.
library(here)
library(maptools) # for readShapePoints() and to write a kml file
# library(rgdal) # for readOGR function; no longer available 2024
library(sf) # use instead of rgdal
library(sp)
library(ggplot2)
# plotKML is no loger available, it went away with rgdal
# library(plotKML) # for "easier" conversion to kml file
# library(devtools)
# install_github("envirometrix/plotKML") # install github version that ingegrates with sf
#

# Get sites "SelectOASDsites.R" output:
# load("~/Desktop/Fish_PostDoc/ArkansasR/ArkSelectSites/ArkAllSiteSelections.RData")
# setwd("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/")
load(paste0(here(),"/Output_Files/8_SelectOASDsites/ArkAllSiteSelections.RData"))
# setwd("/Volumes/CPW/CPW_Work/Optimum_Sampling/Ark_Optimal_Final/Output_Files/9_ArkPlotFutureSites")

# read_sf

## AllSites == all possible sampling sites = OrigSites + pred_sites = 620
## all_sites == more continous view of the basin. 47848 sites

## Turn optimalSites into Spatial object:
optimal1 <- AllSites[OASDsites[1, ], ] # newsites2016
coordinates(optimal1) <- ~ UTMX + UTMY
proj4string(optimal1) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# SP@proj4string

optimal2 <- AllSites[OASDsites[2, ], ] # newsites2017
coordinates(optimal2) <- ~ UTMX + UTMY
proj4string(optimal2) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# SP@proj4string

# for some reason, OASDsites[3, ] didn't save in workspace.  Add here:
min.qInd <- which(qAll[[yr]] == min(qAll[[yr]]))[1]
OASD.qInd[yr] <- min.qInd
OASDq[yr] <- qAll[[yr]][min.qInd]
OASDsites[yr, ] <- SitesAll[[yr]][min.qInd, ]
extraSitesIndex <- c(extraSitesIndex, OASDsites[yr, ])

optimal3 <- AllSites[OASDsites[3, ], ] # newsites2018
coordinates(optimal3) <- ~ UTMX + UTMY
proj4string(optimal3) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# SP@proj4string


# Plot in R to compare to balanced sites:
palette(rainbow(6))
par(mfrow=c(1,1))
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n"))
#points(R, cex=1.5)
points(BalancedSites@coords[, 1], BalancedSites@coords[, 2], 
       col=as.numeric(as.factor(BalancedSites$panel)), pch=16, cex=1.5)
points(optimal1, col=4, pch=8, cex=1.5)
points(optimal2, col=5, pch=8, cex=1.5)
points(optimal3, col=6, pch=8, cex=1.5)

## Order sites of importance/value for 2016 --------------------

# Because Ark sampling has been opportunistic, sampling the balanced sites is most
#  important first. 
# Check out the details of the optimal 2016 sites:
# as.data.frame(optimal1)[, c(1, 6:8, 12, 16:19, 25)]

# Run loop to see what happens when values are added or removed
tmp2016 <- extraSitesIndex[c(31:40)]
tmpOther <- extraSitesIndex[1:30]  # balanced sites from 2017, 2018
tmp.q <- vector(length=length(tmp2016))
for (i in 1:length(tmp2016)) {
print(i)
tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                 AllSitesOccuProbs, AllSurveysDetectProbs,
                 nAllSites, 
                 c(tmp2016[i], tmpOther),
                 nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
print(tmp.q[i])
}
# if q goes up, then site i is least important.
sample(c(1:8, 10), 1)  # site 6 = importance 1





### Export plots for Procedure Manual ---------------------------
palette("default")
pdf(paste0(here(),"/Output_Files/9_ArkPlotFutureSites/ArkFutureSites.pdf"), width=5)# height=5)
par(mar=c(0.1, 0.1, 2.1, 0.1), mfrow=c(3, 1))
# plot(SP, pch=16, cex=0.5, col="gray", main="2016")
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n", main="Year 1"))
# points(subset(balanced.spat, year== "2016"), col=1, pch=16, cex=1.5)
with(subset(BalancedSites@data, panel=="Year1"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
points(optimal1, bg=5, pch=23, cex=1.2, col="black")
# dev.off()
# pdf("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/FutureSites2017.pdf", height=5)
par(mar=c(0.1, 0.1, 2.1, 0.1))
# plot(SP, pch=16, cex=0.5, col="gray", main="2016")
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n", main="Year 2"))
# points(subset(balanced.spat, year== "2017"), col=1, pch=16, cex=1.5)
with(subset(BalancedSites@data, panel=="Year2"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
points(optimal2, bg=5, pch=23, cex=1.2, col="black")
# dev.off()
# pdf("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/FutureSites2018.pdf", height=5)
par(mar=c(0.1, 0.1, 2.1, 0.1))
# plot(SP, pch=16, cex=0.5, col="gray", main="2016")
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n", main="Year 3"))
# points(subset(balanced.spat, year== "2018"), col=1, pch=16, cex=1.5)
with(subset(BalancedSites@data, panel=="Year3"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
points(optimal3, bg=5, pch=23, cex=1.2, col="black")
dev.off()





####
## Output results to Excel ---------------------
####


# CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# 
allSampSites <- rbind(AllSites[which(AllSites$pointid %in% BalancedSites$pointid), ], # bSites
                      AllSites[OASDsites[1, ], ],
                      AllSites[OASDsites[2, ], ],
                      AllSites[OASDsites[3, ], ])
allSampSites$SiteSource <- c(rep("balanced", 41), rep("optimal", 30))# rep(c("balanced", "optimal"), each=30)
allSampSites$year <- c(rep(c(1, 2, 3), each=10), rep("OVERSAMPLE", 11), # messed with oversample number
                       rep(c(1, 2, 3), each=10))
write.csv(allSampSites, 
          paste0(here(),"/Output_Files/9_ArkPlotFutureSites/Ark_FutureSamplingLocations.csv"),
          row.names = FALSE)







#####################################################################################################


library(leaflet)
#original Coord system for inital sf object making
nativeCordinateSystem <- "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
#defining latLong crs for leaflet plotting
latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326
### Year1
ARKbalanced1 <- subset(BalancedSites, panel=="Year1")

ARKbalanced1_sf <-  st_as_sf(ARKbalanced1, coords = c("xcoord", "ycoord"), remove = FALSE)
st_crs(ARKbalanced1_sf) <- nativeCordinateSystem


#transform to new crs for plotting with leaflet
ARKbalanced1_sf <- st_transform(ARKbalanced1_sf, latLongCRS) 
### Year2
ARKbalanced2 <- subset(BalancedSites, panel=="Year2")# subset(balanced.spat, year==2017)
ARKbalanced2_sf <-  st_as_sf(ARKbalanced2, coords = c("xcoord", "ycoord"), remove = FALSE)
st_crs(ARKbalanced2_sf) <- nativeCordinateSystem
#transform to new crs for plotting with leaflet
ARKbalanced2_sf <- st_transform(ARKbalanced2_sf, latLongCRS) 

### Year 3
ARKbalanced3 <- subset(BalancedSites, panel=="Year3")# subset(balanced.spat, year==2017)
ARKbalanced3_sf <-  st_as_sf(ARKbalanced3, coords = c("xcoord", "ycoord"), remove = FALSE)
st_crs(ARKbalanced3_sf) <- nativeCordinateSystem
#transform to new crs for plotting with leaflet
ARKbalanced3_sf <- st_transform(ARKbalanced3_sf, latLongCRS) 

#optimal sites
##2016
optimal1_sf <-  st_as_sf(optimal1)
#transform to new crs for plotting with leaflet
optimal1_sf <- st_transform(optimal1_sf, latLongCRS) 

##2017
optimal2_sf <-  st_as_sf(optimal2)
#transform to new crs for plotting with leaflet
optimal2_sf <- st_transform(optimal2_sf, latLongCRS) 

##2018
optimal3_sf <-  st_as_sf(optimal3)
#transform to new crs for plotting with leaflet
optimal3_sf <- st_transform(optimal3_sf, latLongCRS) 

##stream network 
arkStreamNetwork <- read_sf(file.path(paste0(here(), "/ArcGIS_files/AsShapefiles/arkansasStreamNetwork.shp")))
#remove z element, not needed on maps and causes a little errors 
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) 

map <- leaflet() %>%
  addTiles(options = providerTileOptions(maxZoom = 100), group = "OSM") %>%
  addProviderTiles(providers$Esri.WorldImagery,
                   options = providerTileOptions(maxZoom = 100), 
                   group = "Satellite"
  ) %>%
  addPolylines(data = arkStreamNetwork2,
               group = "Stream Network", 
               popup = paste(
                 "WaterBody:", arkStreamNetwork2$NAME1, "<br>"
               )
  ) %>%
  addAwesomeMarkers(data = ARKbalanced1_sf, 
                    group = "ARKbalanced year 1", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "blue"
                    ), 
                    popup = paste(
                      "ARKbalanced year 1", "<br>", 
                      "Point ID: ", ARKbalanced1_sf$pointid, "<br>", 
                      "UTMX:", ARKbalanced1_sf$xcoord, "<br>",
                      "UTMY:", ARKbalanced1_sf$ycoord, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = ARKbalanced2_sf, 
                    group = "ARKbalanced year 2", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "orange"
                    ), 
                    popup = paste(
                      "ARKbalanced year 2", "<br>", 
                      "Point ID: ", ARKbalanced2_sf$pointid, "<br>", 
                      "UTMX:", ARKbalanced2_sf$xcoord, "<br>",
                      "UTMY:", ARKbalanced2_sf$ycoord, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = ARKbalanced3_sf, 
                    group = "ARKbalanced year 3", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "green"
                    ), 
                    popup = paste(
                      "ARKbalanced year 3", "<br>", 
                      "Point ID: ", ARKbalanced3_sf$pointid, "<br>", 
                      "UTMX:", ARKbalanced3_sf$xcoord, "<br>",
                      "UTMY:", ARKbalanced3_sf$ycoord, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = optimal1_sf, 
                    group = "Optimal year 1", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "cadetblue"
                    ), 
                    popup = paste(
                      "Optimal year 1", "<br>", 
                      "Site ID: ", optimal1_sf$SiteID, "<br>", 
                      "OASD ID:", optimal1_sf$OASDid, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = optimal2_sf, 
                    group = "Optimal year 2", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "red"
                    ), 
                    popup = paste(
                      "Optimal year 2", "<br>", 
                      "Site ID: ", optimal2_sf$SiteID, "<br>", 
                      "OASD ID:", optimal2_sf$OASDid, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = optimal3_sf, 
                    group = "Optimal year 3", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "purple"
                    ), 
                    popup = paste(
                      "Optimal year 3", "<br>", 
                      "Site ID: ", optimal3_sf$SiteID, "<br>", 
                      "OASD ID:", optimal3_sf$OASDid, "<br>"
                    )
  ) %>%
  addLayersControl(overlayGroups = c("Stream Network", "ARKbalanced Year1", "ARKbalanced Year2", "ARKbalanced Year3", "Optimal Year1", "Optimal Year2", "Optimal Year3"), 
                   baseGroups = c("OSM", "Satellite"))




# library(htmlwidgets)
saveWidget(map, file=paste0(here(),"/Output_Files/9_ArkPlotFutureSites/map_of_sites.html"))








# end of file ------------------------------------
