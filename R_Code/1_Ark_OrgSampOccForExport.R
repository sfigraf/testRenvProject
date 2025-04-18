
## The purpose of the file is to clean-up the sampling occasions 
##  spreadsheet. 

# INPUTS:
# Set working directory to Input_Files Location
# 2024 Data Cleaned_Data_include_dat_25Sept2024.csv

###renv code:
# renv::init()
# renv::snapshot()
# OUTPUTS:
# IncludeSampOccasions.dbf of cleaned locations
# IncludeSampOccasions.dbf goes into ArcGIS to get the snapped locations.

# Libraries:
#library(foreign) # to read and write .dbf files

# Set the directory:
setwd(here())

data <- read.csv("Input_Files/Cleaned_Data_include_dat_27Sept2024.csv", stringsAsFactors=F)


### Write for ArcGIS and later Org ------------------------------
# Write .dbf file of the include_SampOcc object
write.dbf(data, "Output_Files/1_OrgSampOccForExport/IncludedSampOccasions.dbf")
#
# write.csv(data, "/Volumes/CPW_Work/Optimum Sampling/Ark_Optimal_bwa/Output_Files/Cleaned_Data_include_dat_25Sept2024.csv", row.names=F)

# end of file. -----------------------------

#base::load("~/testRenvProject/Output_Files/2_OrganizeNewData/ArkData.Rdata")
# base::save.image("ArkData111.Rdata")
