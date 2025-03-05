
# plot sites selected as best.
library(maptools) # for readShapePoints() and to write a kml file
#install.packages("C:/Users/graffs.NATURENET.001/Downloads/maptools_0.8-40.tar.gz", repos = NULL, type = "source")
# library(rgdal) # for readOGR function; no longer available 2024
library(sf) # use instead of rgdal
library(sp)
library(ggplot2)
# plotKML is no loger available, it went away with rgdal
#library(plotKML) # for "easier" conversion to kml file
# library(devtools)
# install_github("envirometrix/plotKML") # install github version that ingegrates with sf
#

# Get sites "SelectOASDsites.R" output:
#base::load("~/Desktop/Fish_PostDoc/ArkansasR/ArkSelectSites/ArkAllSiteSelections.RData")
# setwd("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/")
#SG: not loading .rdata bc it doesn't work as well in Rmarkdown. if all scripts ran in succession, shouldn't be a need to call .Rdata files
base::load("~/testRenvProject/Output_Files/8_SelectOASDsites/ArkAllSiteSelections.RData")
setwd(file.path(here(), "Output_Files/9_ArkPlotFutureSites"))

# read_sf

## AllSites == all possible sampling sites = OrigSites + pred_sites = 620
## all_sites == more continous view of the basin. 47848 sites

## Turn optimalSites into Spatial object:
optimal2016 <- AllSites[OASDsites[1, ], ] # newsites2016
coordinates(optimal2016) <- ~ UTMX + UTMY
proj4string(optimal2016) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# SP@proj4string

optimal2017 <- AllSites[OASDsites[2, ], ] # newsites2017
coordinates(optimal2017) <- ~ UTMX + UTMY
proj4string(optimal2017) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# SP@proj4string

# for some reason, OASDsites[3, ] didn't save in workspace.  Add here:
min.qInd <- which(qAll[[yr]] == min(qAll[[yr]]))[1]
OASD.qInd[yr] <- min.qInd
OASDq[yr] <- qAll[[yr]][min.qInd]
OASDsites[yr, ] <- SitesAll[[yr]][min.qInd, ]
extraSitesIndex <- c(extraSitesIndex, OASDsites[yr, ])

optimal2018 <- AllSites[OASDsites[3, ], ] # newsites2018
coordinates(optimal2018) <- ~ UTMX + UTMY
proj4string(optimal2018) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# SP@proj4string


# Plot in R to compare to balanced sites:
palette(rainbow(6))
par(mfrow=c(1,1))
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n"))
#points(R, cex=1.5)
points(BalancedSites@coords[, 1], BalancedSites@coords[, 2], 
       col=as.numeric(as.factor(BalancedSites$panel)), pch=16, cex=1.5)
points(optimal2016, col=4, pch=8, cex=1.5)
points(optimal2017, col=5, pch=8, cex=1.5)
points(optimal2018, col=6, pch=8, cex=1.5)

## Order sites of importance/value for 2016 --------------------

# Because Ark sampling has been opportunistic, sampling the balanced sites is most
#  important first. 
# Check out the details of the optimal 2016 sites:
# as.data.frame(optimal2016)[, c(1, 6:8, 12, 16:19, 25)]

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
tmpOther <- c(tmpOther, tmp2016[6])
for (i in c(1:5, 7:10)) {
  print(i)
  tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                               AllSitesOccuProbs, AllSurveysDetectProbs,
                               nAllSites, 
                               c(tmp2016[i], tmpOther),
                               nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
  print(tmp.q[i])
}
# can pick anything but site 9. Try site 5.
tmpOther <- c(tmpOther, tmp2016[5])
for (i in c(1:4, 7:10)) {
  print(i)
  tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                               AllSitesOccuProbs, AllSurveysDetectProbs,
                               nAllSites, 
                               c(tmp2016[i], tmpOther),
                               nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
  print(tmp.q[i])
}
# can pick anything but site 9. Try site 10.
tmpOther <- c(tmpOther, tmp2016[10])
for (i in c(1:4, 8:9)) {
  print(i)
  tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                               AllSitesOccuProbs, AllSurveysDetectProbs,
                               nAllSites, 
                               c(tmp2016[i], tmpOther),
                               nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
  print(tmp.q[i])
}
# can pick anything but site 9. Try site 8.
tmpOther <- c(tmpOther, tmp2016[8])
for (i in c(1:4, 9)) {
  print(i)
  tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                               AllSitesOccuProbs, AllSurveysDetectProbs,
                               nAllSites, 
                               c(tmp2016[i], tmpOther),
                               nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
  print(tmp.q[i])
}
# can pick anything but site 9. Try site 2.
tmpOther <- c(tmpOther, tmp2016[2])
for (i in c(1, 3:4, 9)) {
  print(i)
  tmp.q[i] <- sum(apply(calc_q(repSppIndices, mod, 
                               AllSitesOccuProbs, AllSurveysDetectProbs,
                               nAllSites, 
                               c(tmp2016[i], tmpOther),
                               nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
  print(tmp.q[i])
}




### Export plots for Procedure Manual ---------------------------
palette("default")
pdf("ArkFutureSites.pdf", width=5)# height=5)
par(mar=c(0.1, 0.1, 2.1, 0.1), mfrow=c(3, 1))
# plot(SP, pch=16, cex=0.5, col="gray", main="2016")
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n", main="2016"))
# points(subset(balanced.spat, year== "2016"), col=1, pch=16, cex=1.5)
with(subset(BalancedSites@data, panel=="Year1"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
points(optimal2016, bg=5, pch=23, cex=1.2, col="black")
# dev.off()
# pdf("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/FutureSites2017.pdf", height=5)
par(mar=c(0.1, 0.1, 2.1, 0.1))
# plot(SP, pch=16, cex=0.5, col="gray", main="2016")
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n", main="2017"))
# points(subset(balanced.spat, year== "2017"), col=1, pch=16, cex=1.5)
with(subset(BalancedSites@data, panel=="Year2"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
points(optimal2017, bg=5, pch=23, cex=1.2, col="black")
# dev.off()
# pdf("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/FutureSites2018.pdf", height=5)
par(mar=c(0.1, 0.1, 2.1, 0.1))
# plot(SP, pch=16, cex=0.5, col="gray", main="2016")
with(all_sites, plot(UTMX, UTMY, pch=16, col="gray", yaxt="n", xaxt="n", main="2018"))
# points(subset(balanced.spat, year== "2018"), col=1, pch=16, cex=1.5)
with(subset(BalancedSites@data, panel=="Year3"), 
     points(xcoord, ycoord, col=1, pch=16, cex=1.5))
points(optimal2018, bg=5, pch=23, cex=1.2, col="black")
dev.off()


## Output results to Excel and google Earth ---------------------

# CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")# 
allSampSites <- rbind(AllSites[which(AllSites$pointid %in% BalancedSites$pointid), ], # bSites
                      AllSites[OASDsites[1, ], ],
                      AllSites[OASDsites[2, ], ],
                      AllSites[OASDsites[3, ], ])
allSampSites$SiteSource <- c(rep("balanced", 40), rep("optimal", 30))# rep(c("balanced", "optimal"), each=30)
allSampSites$year <- c(rep(c(2016, 2017, 2018), each=10), rep("OVERSAMPLE", 10), 
                       rep(c(2016, 2017, 2018), each=10))
write.csv(allSampSites, 
          "Ark_FutureSamplingLocations2Nstarts_scriptonly.csv",
          row.names = FALSE)


############ Good Up to Here ################
library(leaflet)
#original Coord system for inital sf object making
nativeCordinateSystem <- "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
#defining latLong crs for leaflet plotting
latLongCRS <- st_crs("+proj=longlat +datum=WGS84 +no_defs") #should be same as +init=epsg:4326
###2016
ARKbalanced2016 <- subset(BalancedSites, panel=="Year1")

ARKbalanced2016_sf <-  st_as_sf(ARKbalanced2016, coords = c("xcoord", "ycoord"), remove = FALSE)
st_crs(ARKbalanced2016_sf) <- nativeCordinateSystem


#transform to new crs for plotting with leaflet
ARKbalanced2016_sf <- st_transform(ARKbalanced2016_sf, latLongCRS) 
###2017
ARKbalanced2017 <- subset(BalancedSites, panel=="Year2")# subset(balanced.spat, year==2017)
ARKbalanced2017_sf <-  st_as_sf(ARKbalanced2017, coords = c("xcoord", "ycoord"), remove = FALSE)
st_crs(ARKbalanced2017_sf) <- nativeCordinateSystem
#transform to new crs for plotting with leaflet
ARKbalanced2017_sf <- st_transform(ARKbalanced2017_sf, latLongCRS) 

###2018
ARKbalanced2018 <- subset(BalancedSites, panel=="Year3")# subset(balanced.spat, year==2017)
ARKbalanced2018_sf <-  st_as_sf(ARKbalanced2018, coords = c("xcoord", "ycoord"), remove = FALSE)
st_crs(ARKbalanced2018_sf) <- nativeCordinateSystem
#transform to new crs for plotting with leaflet
ARKbalanced2018_sf <- st_transform(ARKbalanced2018_sf, latLongCRS) 

#optimal sites
##2016
optimal2016_sf <-  st_as_sf(optimal2016)
#transform to new crs for plotting with leaflet
optimal2016_sf <- st_transform(optimal2016_sf, latLongCRS) 

##2017
optimal2017_sf <-  st_as_sf(optimal2017)
#transform to new crs for plotting with leaflet
optimal2017_sf <- st_transform(optimal2017_sf, latLongCRS) 

##2018
optimal2018_sf <-  st_as_sf(optimal2018)
#transform to new crs for plotting with leaflet
optimal2018_sf <- st_transform(optimal2018_sf, latLongCRS) 

##stream network 
arkStreamNetwork <- read_sf(file.path(paste0(here(), "/ArcGIS_files/AsShapefiles/arkansasStreamNetwork.shp")))
#remove z element, not needed on maps and causes a little errors 
arkStreamNetwork1 <- st_zm(arkStreamNetwork, drop = TRUE, what = "ZM")
arkStreamNetwork2 <- st_transform(arkStreamNetwork1, latLongCRS) 

leaflet() %>%
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
  addAwesomeMarkers(data = ARKbalanced2016_sf, 
                    group = "ARKbalanced 2016", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "blue"
                    ), 
                    popup = paste(
                      "ARKbalanced 2016", "<br>", 
                      "Point ID: ", ARKbalanced2016_sf$pointid, "<br>", 
                      "UTMX:", ARKbalanced2016_sf$xcoord, "<br>",
                      "UTMY:", ARKbalanced2016_sf$ycoord, "<br>"
                    )
                    ) %>%
  addAwesomeMarkers(data = ARKbalanced2017_sf, 
                    group = "ARKbalanced 2017", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "orange"
                    ), 
                    popup = paste(
                      "ARKbalanced 2017", "<br>", 
                      "Point ID: ", ARKbalanced2017_sf$pointid, "<br>", 
                      "UTMX:", ARKbalanced2017_sf$xcoord, "<br>",
                      "UTMY:", ARKbalanced2017_sf$ycoord, "<br>"
                    )
                    ) %>%
  addAwesomeMarkers(data = ARKbalanced2018_sf, 
                    group = "ARKbalanced 2018", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "green"
                    ), 
                    popup = paste(
                      "ARKbalanced 2018", "<br>", 
                      "Point ID: ", ARKbalanced2018_sf$pointid, "<br>", 
                      "UTMX:", ARKbalanced2018_sf$xcoord, "<br>",
                      "UTMY:", ARKbalanced2018_sf$ycoord, "<br>"
                    )
                    ) %>%
  addAwesomeMarkers(data = optimal2016_sf, 
                    group = "Optimal 2016", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "cadetblue"
                    ), 
                    popup = paste(
                      "Optimal 2016", "<br>", 
                      "Site ID: ", optimal2016_sf$SiteID, "<br>", 
                      "OASD ID:", optimal2016_sf$OASDid, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = optimal2017_sf, 
                    group = "Optimal 2017", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "red"
                    ), 
                    popup = paste(
                      "Optimal 2017", "<br>", 
                      "Site ID: ", optimal2017_sf$SiteID, "<br>", 
                      "OASD ID:", optimal2017_sf$OASDid, "<br>"
                    )
  ) %>%
  addAwesomeMarkers(data = optimal2018_sf, 
                    group = "Optimal 2018", 
                    icon = leaflet::awesomeIcons(
                      icon = 'add',
                      library = 'ion',
                      markerColor = "purple"
                    ), 
                    popup = paste(
                      "Optimal 2018", "<br>", 
                      "Site ID: ", optimal2018_sf$SiteID, "<br>", 
                      "OASD ID:", optimal2018_sf$OASDid, "<br>"
                    )
  ) %>%
  addLayersControl(overlayGroups = c("Stream Network", "ARKbalanced 2016", "ARKbalanced 2017", "ARKbalanced 2018", "Optimal 2016", "Optimal 2017", "Optimal 2018"), 
                   baseGroups = c("OSM", "Satellite"))
  
## Check-out balanced sites in google Earth:
# setwd("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/")
# ARKbalanced2016 <- subset(BalancedSites, panel=="Year1")# subset(balanced.spat, year==2016)
# ARKbalanced2016 <- SpatialPointsDataFrame(ARKbalanced2016@coords, ARKbalanced2016@data,
#   proj4string=CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# 
# ggplot(ARKbalanced2016, aes(x=data$coords[,1],y=data$coords[,2])) + 
#   geom_polygon() +
#   geom_path(color="gray")
# 
# 
# plotKML(ARKbalanced2016["geometry"], 
#         points_names="2016", 
#         colour_scale="#00FF00",  # lime
#         balloon=TRUE)
# ARKbalanced2017 <- subset(BalancedSites, panel=="Year2")# subset(balanced.spat, year==2017)
# ARKbalanced2017 <- SpatialPointsDataFrame(ARKbalanced2017@coords, ARKbalanced2017@data,
#    proj4string=CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# plotKML(ARKbalanced2017, 
#         points_names="2017", 
#         colour_scale="#00FFFF", # aqua
#         balloon=TRUE)
# ARKbalanced2018 <- subset(BalancedSites, panel=="Year3")# subset(balanced.spat, year==2017)
# ARKbalanced2018 <- SpatialPointsDataFrame(ARKbalanced2018@coords, ARKbalanced2018@data,
#   proj4string=CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
# plotKML(ARKbalanced2018, 
#         points_names="2018", 
#         colour_scale="#FF00FF",
#         balloon=TRUE)
# 
# ## And add optimal sites to G. Earth:
# plotKML(optimal2016, file.name = "ARKoptimal2016.kml",
#         points_names="2016", 
#         colour_scale="#008000",  # green (=darker lime)
#         balloon=TRUE)
# plotKML(optimal2017, file.name = "ARKoptimal2017.kml",
#         points_names="2017", 
#         colour_scale="#008080",  # teal (=darker aqua)
#         balloon=TRUE)
# plotKML(optimal2018, file.name = "ARKoptimal2018.kml",
#         points_names="2018", 
#         colour_scale="#800080",  # purple (=darker fuchsia)
#         balloon=TRUE)
## change point and label sizes inside of google earth.
## (R to kml functions don't work that well)




### ------------------------------
## Sites are sorta redundant... What happens if you scale back on redundant sites?

# OASDsites[newsites2016, ]
# round(iDist(OASDsites[newsites2016, c("x_coord", "y_coord")]))
# # rm site 10- very close to site 3 and 7
# OASDsites[newsites2016[-10], ]
# ## q value with all 10 optimal sites:
# new.q <- calc_q(repSppIndices, mod, occu.mat, det.mat,
#                 nPredSites, c(extraSitesIndex, newsites2016),
#                 nSurveys, nMCMC=nMCMC)
# sum( apply(new.q, 1, sum) )  # 258.4468
# # q value without site 10
# new.q <- calc_q(repSppIndices, mod, occu.mat, det.mat,
#                 nPredSites, c(extraSitesIndex, newsites2016[-c(1, 10)], 330, 300),
#                 nSurveys, nMCMC=nMCMC)
# sum( apply(new.q, 1, sum) )  # 259.1267





### Plot of Random Sampling vs OASD   ----------------------------

# and for a random design search
# set.seed(2016)
# RandomSites2016 <- sample(1:nAllSites, nAdd, replace=T)
# random2016 <- sum(apply(calc_q(repSppIndices, mod, 
#                              AllSitesOccuProbs, AllSurveysDetectProbs,
#                              nAllSites, 
#                              c(extraSitesIndex[1:30], RandomSites2016),
#                              nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# random2016  # 725.1816
# (random2016 / length(repSppIndices) ) / nAllSites  # 0.12996
# 
# set.seed(2017)
# RandomSites2017 <- sample(1:nAllSites, nAdd, replace=T)
# random2017 <- sum(apply(calc_q(repSppIndices, mod, 
#                                AllSitesOccuProbs, AllSurveysDetectProbs,
#                                nAllSites, 
#                                c(extraSitesIndex[1:40], RandomSites2017),
#                                nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# random2017  # 710.1885
# (random2017 / length(repSppIndices) ) / nAllSites  # 0.127
# ###
# random2017b <- sum(apply(calc_q(repSppIndices, mod, 
#                                AllSitesOccuProbs, AllSurveysDetectProbs,
#                                nAllSites, 
#                                c(extraSitesIndex[1:30], RandomSites2016, RandomSites2017),
#                                nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# random2017b  # 716.9112
# (random2017b / length(repSppIndices) ) / nAllSites  # 0.1285
# ###
# set.seed(2018)
# RandomSites2018 <- sample(1:nAllSites, nAdd, replace=T)
# random2018 <- sum(apply(calc_q(repSppIndices, mod, 
#                                AllSitesOccuProbs, AllSurveysDetectProbs,
#                                nAllSites, 
#                                c(extraSitesIndex[1:50], RandomSites2018),
#                                nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# random2018  # 703.5089
# (random2018 / length(repSppIndices) ) / nAllSites  # 0.126
# random2018b <- sum(apply(calc_q(repSppIndices, mod, 
#                                AllSitesOccuProbs, AllSurveysDetectProbs,
#                                nAllSites, 
#                                c(extraSitesIndex[1:30], RandomSites2016, 
#                                  RandomSites2017, RandomSites2018),
#                                nSurveys=nSurveys, nMCMC=nMCMC), 1, sum))
# random2018b  # 710.34
# (random2018b / length(repSppIndices) ) / nAllSites  # 0.127
# 
# # setwd("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/FutureSites/")
# 
# pdf("CompareQ.pdf", height=5)
# par(mar=c(4, 4, 3, 1))
# plot(OASDq, xaxt="n", type="o", 
#      xlab="Number of additional sampling occasions", 
#      ylab="q(d)", las=1, ylim=c(690, 732), lwd=2,
#      main="q under different sampling schemes")
# points(c(random2016, random2017b, random2018b), 
#        type="b", lty=2, pch=8, lwd=3)
# mtext(c(10, 20, 30), 1, at=1:3, line=1)
# abline(h=high.q, col="gray", lwd=3, lty=4)
# legend(1, 700, c("Simple Random Sampling", "Optimal Sampling"),
#        lty=2:1, pch=c(8, 1), bty="n", lwd=2)
# dev.off()
#

# end of file ------------------------------------
