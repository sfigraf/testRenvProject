
#  Process results for best fitting model
# (create tables, figures to show output)
library(mcmcplots)  # to check convergence
library(sf)
library(sp)
library(boot)  # for inv.logit function
library(lattice)  # for trellis.par.set() fxn to change plot colors
library(xtable)  # easily export table to Latex
library(RColorBrewer)  # to change palette for better readability
# display.brewer.all()

load("~/CPWOptimalSampling//Output_Files/4_RunModels/ArkInSampleModels.RData")
setwd("~/CPWOptimalSampling//Output_Files/")
mod.no <- 25
mod <- out[[mod.no]]
nPsiCoef <- nPsiList[[mod.no]]
nDetCoef <- nDetList[[mod.no]]

### Check results -----------------
# str(mod, 1)
# mcmcplot(mod, c("muPsiCoef", "sigmaPsiCoef"))
# mcmcplot(mod, "psiCoef")
# mcmcplot(mod, c("muDetCoef", "sigmaDetCoef"))
# mcmcplot(mod, "detCoef")

### Set-up for maps ------------------
pred_sites$X2 <- pred_sites$X ^ 2 
pred_sites$Y2 <- pred_sites$Y ^ 2
pred_sites$aboveP <- ifelse(pred_sites$RESERVOIR=="abovePueblo", 1, 0)
pred_sites$belowJM <- ifelse(pred_sites$RESERVOIR=="belowJM", 1, 0)
# psiX.mat <- cbind(1, pred_sites[, c("Y", "Y2", "ELEV", "ELEV2",
#                                     "UNCONNECT", "aboveP", "belowJM")])
psiX.mat <- cbind(1, pred_sites[, c("SIZE")])

## Create a more continuous map:
all_sites$X2 <- all_sites$X ^ 2 
all_sites$Y2 <- all_sites$Y ^ 2
all_sites$aboveP <- ifelse(all_sites$RESERVOIR=="abovePueblo", 1, 0)
all_sites$belowJM <- ifelse(all_sites$RESERVOIR=="belowJM", 1, 0)
# psiX.mat <- cbind(1, all_sites[, c("Y", "Y2", "ELEV", "ELEV2",
#                                     "UNCONNECT", "aboveP", "belowJM")])
psiX.mat <- cbind(1, all_sites[, c("SIZE")])


## Representative Species??
## ONLY for maps, NOT for design criterion
spp.ind <- c(which(spp %in% c("ARD", "LND", "SMM")), 17)  #17=RCS
spp.out <- c("ARD", "LND", "SMM", "RCS")
fullSppName <- c("Arkansas Darter", "Longnose Dace", 
                 "Suckermouth Minnow", "River Carpsucker")

### Create occupancy predictions for all natives ---------------------
# occuMap.se <- occuMap <- matrix(NA, nrow=nrow(pred_sites), ncol=nPossSpp)
occuMap.se <- occuMap <- matrix(NA, nrow=nrow(all_sites), ncol=nPossSpp)

for (j in 1:nPossSpp) {
  occu.mat <- as.matrix(psiX.mat) %*% 
    t(mod$BUGSoutput$sims.list$psiCoef[, j, ])
  occuMap[, j] <- apply(inv.logit(occu.mat), 1, mean)
  occuMap.se[, j] <- apply(inv.logit(occu.mat), 1, sd)
}
# predmaps <- data.frame(id=pred_sites$pointid, 
#                        x_coord=pred_sites$UTMX, 
#                        y_coord=pred_sites$UTMY, 
#                        occuMap, occuMap.se )
predmaps <- data.frame(id=all_sites$pointid, 
                       SIZE=all_sites$SIZE,
                       x_coord=all_sites$UTMX, 
                       y_coord=all_sites$UTMY, 
                       occuMap, occuMap.se )
predmaps <- predmaps[order(predmaps$SIZE), ]
names(predmaps)[5:(length(natives)+4)] <- paste(natives, "_Occupancy", sep="")
names(predmaps)[21] <- "RCS_Occupancy"
names(predmaps)[22:(21+length(natives))] <- paste(natives, "sd", sep="")
names(predmaps)[38] <- "RCSsd"
# names(predmaps)
coordinates(predmaps) <- c("x_coord", "y_coord")
proj4string(predmaps) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
### Create different OCCUPANCY maps ----------------

# Fig of Occu Probs 3 rep spp + 1 undetected for main text:
# set_col_regions(rev(grey.colors(100)))
my.palette <- brewer.pal(n = 7, name = "YlOrRd")
my.palette <- c("#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
                "#e31a1c", "#bd0026", "#800026", "#4c0016", "#33000f")
# colors from  http://colorbrewer2.org/ and http://www.color-hex.com/

png("~/CPWOptimalSampling/Output_Files/5_ModelResults/OccuProb_maps.png", 
    width=1024, height=1000)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, c("ARD_Occupancy", "LND_Occupancy", "SMM_Occupancy", 
                   "RCS_Occupancy"), 
       colorkey=T, as.table=T, cex=1, layout=c(2, 2), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[1:7], cuts = 6)  # cuts must be n - 1
dev.off()
#

# Maps for a specific species:
# jpeg("~/Dropbox/Fish_PostDoc/MultiSpp/write_up/LAC_maps.jpg", height=700)
# spplot(predmaps, c("LAC", "LACsd"), colorkey=T, as.table=T, cex=1)
# dev.off()
#

# Maps for ALL species!!
#set_col_regions(bpy.colors(100))
png("~/CPWOptimalSampling//Output_Files/5_ModelResults/OccuProb_mapsALL_1.png",
    width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, 3:8, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette, cuts=9)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/OccuProb_mapsALL_2.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, 9:14, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette, cuts=9)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/OccuProb_mapsALL_3.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, 15:19, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette, cuts=9)
dev.off()

# SD Maps for ALL species!!
png("~/CPWOptimalSampling//Output_Files/5_ModelResults/OccuProbSD_mapsALL_1.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, 20:25, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette, cuts=9)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/OccuProbSD_mapsALL_2.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, 26:31, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette, cuts=9)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/OccuProbSD_mapsALL_3.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps, 32:36, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette, cuts=9)
dev.off()

### Create DETECTION predictions for all natives ----------------
# need to make some survey-level assumptions:
# use year =2010
# Method = electrofishing
# YDAY = median(YDAY) = 0 because variable is scaled
# Total Ct = mean( log( Total Ct) ).  because makes sense.
# detXmat <- cbind(1, YDAY=0, YDAY2=0, SEINE=0, PASSNO=1,
#                  y2009=0, y2010=0, y2011=0, y2012=0, y2013=0,
#                  y2014=0, y2015=0, SIZE=pred_sites$SIZE)
detXmat <- cbind(1, SEINE=0, 
 #                TotalCt = mean(surveys$logTotalCt), 
                 SIZE=pred_sites$SIZE)
detMap.se <- detMap <- matrix(NA, nrow=nrow(pred_sites), ncol=nPossSpp)
detXmat <- cbind(1, SEINE=0, 
  #               TotalCt = mean(surveys$logTotalCt), 
                 SIZE=all_sites$SIZE)
detMap.se <- detMap <- matrix(NA, nrow=nrow(all_sites), ncol=nPossSpp)
for(j in 1:nPossSpp){
  det.mat <- as.matrix(detXmat) %*% 
    t(mod$BUGSoutput$sims.list$detCoef[, j, ])
  detMap[, j] <- apply(inv.logit(det.mat), 1, mean)
  detMap.se[, j] <- apply(inv.logit(det.mat), 1, sd)
}  # loop takes awhile

# predmaps.det <- data.frame(id=pred_sites$pointid, 
#                            x_coord=pred_sites$UTMX, 
#                            y_coord=pred_sites$UTMY,
#                            detMap, detMap.se )
predmaps.det <- data.frame(id=all_sites$pointid, 
                           SIZE=all_sites$SIZE,
                           x_coord=all_sites$UTMX, 
                           y_coord=all_sites$UTMY,
                           detMap, detMap.se )
predmaps.det <- predmaps.det[order(predmaps.det$SIZE), ]
names(predmaps.det)[5:(length(natives)+4)] <- paste(natives, "_Detection", sep="")
names(predmaps.det)[21] <- c("RCS_Detection")
names(predmaps.det)[22:(21+length(natives))] <- paste(natives, "sd", sep="")
names(predmaps.det)[38] <- "RCSsd"
# names(predmaps.det)
coordinates(predmaps.det) <- c("x_coord", "y_coord")
proj4string(predmaps.det) <- CRS("+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


### Create different DETECTION maps ------------------
# Fig of Det Probs  5 rep spp + 1 undetected for main text:
png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProb_maps.png", 
    width=850, height=1024)
spplot(predmaps.det, c("ARD_Detection", "LND_Detection", 
                       "SMM_Detection", "RCS_Detection"), 
       colorkey=T, as.table=T, cex=1, layout=c(2, 2), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(2, 3, 5, 7, 9)], cuts=4)
dev.off()

# Maps for a specific species:
# jpeg("~/Dropbox/Fish_PostDoc/MultiSpp/write_up/LAC_maps.jpg", height=700)
# spplot(predmaps.det, c("LAC", "LACsd"), colorkey=T, as.table=T, cex=1)
# dev.off()
#
# Maps for ALL species!!
png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProb_mapsALL_1.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps.det, 3:8, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(3, 5, 7, 8, 10)], cuts=4)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProb_mapsALL_2.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps.det, 9:14, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(1, 2, 3, 5, 7, 8, 10)], cuts=6)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProb_mapsALL_3.png", width=850, height=1024)
trellis.par.set(strip.background=list(col="gray80"))
spplot(predmaps.det, 15:19, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(3, 5, 7, 8, 10)], cuts=4)
dev.off()

# SD Maps for ALL species!!
png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProbSD_mapsALL_1.png", width=850, height=1024)
spplot(predmaps.det, 20:25, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(1, 2, 3, 5, 7)], cuts=4)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProbSD_mapsALL_2.png", width=850, height=1024)
spplot(predmaps.det, 26:31, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(1, 2, 3, 5, 7, 8)], cuts=4)
dev.off()

png("~/CPWOptimalSampling//Output_Files/5_ModelResults/DetProbSD_mapsALL_3.png", width=850, height=1024)
spplot(predmaps.det, 32:36, colorkey=T, as.table=T, cex=0.8, layout=c(2, 3), 
       par.settings=list(fontsize=list(text=26)),
       col.regions = my.palette[c(1, 2, 3, 5, 7, 8, 10)], cuts=6)
dev.off()

### Compare Occu maps to Count maps ------------------------------
#png("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/ArkModelResults/Detects1.png", width=850, height=1024)
# png("~/Dropbox/Fish_PostDoc/ArkansasR/write_up/ArkModelResults/Detects3.png", 
#     width=1024, height=1400)
# par(mfrow=c(3, 2), mar=c(0.5, 1, 3, 1))
# for(i in 13:16){
# spp.ind <- i  # which(natives == "SMM")) ## pick a Species
# sumDetects <- apply(Y[spp.ind, , ], 1, sum, na.rm=T)
# #nSurveys <- apply(Y[1, , ], 1, function(Y) sum(!is.na(Y)))
# # pch.type <- ifelse(sumDetects > 0, 1, 16)
# # cex.type <- ifelse(sumDetects > 0, 4, 1.5)
# pch.type <- ifelse(sumDetects > 0, 16, 16)
# cex.type <- ifelse(sumDetects > 0, 3, 0)
# col.type <- ifelse(sumDetects > 0, "black", "navy")
# plot(all_sites$UTMX, all_sites$UTMY, pch=16, cex=0.5, col="gray", 
#      xaxt="n", yaxt="n")
# with(countmap, points(sites$UTMX, sites$UTMY, 
#                       pch=pch.type, cex=cex.type, col=col.type))
# title(paste(natives[spp.ind], "detected", sep=" "), cex.main=2)
# }
# dev.off()


### Create tables of parameter estimates for each native -----------------
# OccNames <- c("Intercept", "Y", "Y2", "CROPS", "DVLPD",
#               "WTLNDS", "SIZE")
# DetNames <- c("Intercept", "YDAY", "YDAY2", "SEINE", "PASSNO", 
#               "y2009","y2010", "y2011", "y2012", 
#               "y2013", "y2014", "y2015", "SIZE")
OccNames <- c("(Intercept)", "\\texttt{size}")
DetNames <- c("(Intercept)", "\\texttt{seine}", "\\texttt{total ct}", "\\texttt{size}")







### Print spp output tables to LaTeX ---------------------
nativeNames <- c("Arkansas Darter", "Black Bullhead", "Central Stoneroller",
                 "Channel Catfish", "Fathead Minnow", "Flathead Chub",
                 "Green Sunfish", "Longnose Dace", "Orangespotted Sunfish",
                 "Plains Killifish", "Plains Minnow", "Red Shiner", 
                 "Sand Shiner", "Southern Redbelly Dace", "Suckermouth Minnow",
                 "White Sucker", "River Carpsucker")
natives[17] <- "RCS"
for(n in 1:length(natives)){
  out.file <- paste("~/CPWOptimalSampling//Output_Files/5_ModelResults/", natives[n], "table.tex", sep="")

  psiCoef <- round(apply(mod$BUGSoutput$sims.list$psiCoef[, n, ],
                         2, mean), 2)  # n.saved X nPossSpp X nPsiCoef  
  psiCoef.sd <- round(apply(mod$BUGSoutput$sims.list$psiCoef[, n, ],
                            2, sd), 2)  # n.saved X nPossSpp X nPsiCoef  
  detCoef <- round(apply(mod$BUGSoutput$sims.list$detCoef[, n, ],
                         2, mean), 2)  # n.saved X nPossSpp X nPsiCoef  
  detCoef.sd <- round(apply(mod$BUGSoutput$sims.list$detCoef[, n, ],
                            2, sd), 2)  # n.saved X nPossSpp X nPsiCoef    

  sink(out.file, append=FALSE)
  cat("\\begin{table}[ht]
\\centering
\\caption{", nativeNames[n], " (", natives[n], ") parameter estimates and standard deviations.} 
\\vspace{0.2cm}
\\begin{tabular}{ lrr }
      \\toprule 
      Occupancy: & &  \\\\
      \\hline 
      & Mean & (SD)  \\\\  ")
  for(i in 1:nPsiCoef){
    cat(OccNames[i], " & ", psiCoef[i], " & (",  psiCoef.sd[i], 
        ")  \\\\ \n", sep="")
  }     
  cat("\\hline 
      \\noalign{\\smallskip}
      Detection: & &  \\\\
      \\hline 
      & Mean & (SD)  \\\\  \n")
  for(i in 1:nDetCoef){
    cat(DetNames[i], " & ", detCoef[i], " & (",  detCoef.sd[i], 
        ") \\\\ \n", sep="")
  }  
  cat("
      \\bottomrule
      \\end{tabular} 
      \\end{table}")
  sink() 
}
#


### Create table of detection prob's ----------------------
natives[17] <- "RCS"
# detXmat <- cbind(1, YDAY=0, YDAY2=0, SEINE=0, PASSNO=1,
#                  y2009=0, y2010=0, y2011=0, y2012=0, y2013=0,
#                  y2014=0, y2015=0, SIZE=1:4)
detXmat <- cbind(1, SEINE=0, SIZE=1:4) #TOTAL_CT=mean(surveys$logTotalCt), SIZE=1:4)
detMap.se <- detMap <- matrix(NA, nrow=4, ncol=length(natives))
for(j in 1:(length(natives))){
  det.mat <- as.matrix(detXmat) %*% 
    t(mod$BUGSoutput$sims.list$detCoef[, j, ])
  detMap[, j] <- apply(inv.logit(det.mat), 1, mean)
  detMap.se[, j] <- apply(inv.logit(det.mat), 1, sd)
}  # loop takes awhile

detOut <- data.frame(SPP=natives, round(t(detMap), 2))#, detMap.se )
names(detOut)[2:5] <- c("size1", "size2", "size3", "size4")
xtmp <- xtable(detOut,
               caption="An example of the detection probabilities associated with 
                one electrofishing pass.
               Streams of size 1 are lower order streams with lower flows, and 
               stream size 4 represents the main stem of the Arkansas River.",
               label="ExampleDetectProbs", auto=T)
sink("~/CPWOptimalSampling//Output_Files/5_ModelResults/ExampleDetectProbsEF.tex")
print(xtmp, include.rownames=F, caption.placement="top") 
sink()

## Repeat above for both SEINE and EF.


# end of file.


