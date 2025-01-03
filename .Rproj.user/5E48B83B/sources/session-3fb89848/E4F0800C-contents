start_time <- Sys.time()
gcinfo(TRUE)
source("global.R")
end_time <- Sys.time()

print(paste("Reading in libraries took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

start_time <- Sys.time()
print("Script 1")

source("R_Code/1_Ark_OrgSampOccForExport.R")
end_time <- Sys.time()
print(paste("Script 1 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))


start_time <- Sys.time()
print("Script 2")
setwd("~/testRenvProject/")
source("R_Code/2_Ark_OrganizeDataNew.R")
save.image("~/testRenvProject/Output_Files/2_OrganizeNewData/ArkData.Rdata")

end_time <- Sys.time()


print(paste("Script 2 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))



# start_time <- Sys.time()
# print("Script 3")
# source("R_Code/3_RunCrossVal.R")
# end_time <- Sys.time()

#print(paste("Script 3 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))


start_time <- Sys.time()
print("Script 4")
setwd("~/testRenvProject/")
source("R_Code/4_RunModels.R")
save.image("~/testRenvProject/Output_Files/4_RunModels/ArkInSampleModels.RData") # save often in case it crashes

end_time <- Sys.time()

print(paste("Script 4 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

start_time <- Sys.time()
print("Script 5")
setwd("~/testRenvProject/")
source("R_Code/5_ModelResults.R")
end_time <- Sys.time()

print(paste("Script 5 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

start_time <- Sys.time()
print("Script 6")
setwd("~/testRenvProject/")
source("R_Code/6_BestModResidPlots.R")
end_time <- Sys.time()

print(paste("Script 6 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

start_time <- Sys.time()
print("Script 7")
setwd("~/testRenvProject/")
source("R_Code/7_SelectBalancedSites.R")
save(BalancedSites, file="~/testRenvProject/Output_Files/7_SelectBalancedSites/ArkBalancedSites.RData")

end_time <- Sys.time()

print(paste("Script 7 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

start_time <- Sys.time()
print("Script 8")
setwd("~/testRenvProject/")
rm(list = ls())
base::load("~/testRenvProject/Output_Files/2_OrganizeNewData/ArkData.Rdata")
base::load("~/testRenvProject/Output_Files/4_RunModels/ArkInSampleModels.RData")
base::load("~/testRenvProject/Output_Files/7_SelectBalancedSites/ArkBalancedSites.RData") # base::load("ArkBalancedSites.RData")

source("R_Code/8_SelectOASDsites.R")
end_time <- Sys.time()

print(paste("Script 8 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

start_time <- Sys.time()
print("Script 9")
setwd("~/testRenvProject/")
source("R_Code/9_ArkPlotFutureSites.R")
end_time <- Sys.time()
print(paste("Script 9 took", round(difftime(end_time, start_time, units = "mins"),2), "minutes."))

