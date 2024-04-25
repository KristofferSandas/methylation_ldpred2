# Script used to extract the illumina manifests for EPIC1 and
# 450K from the comeback package.

# Load comeback system data
sysdata_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                     "For_Kristoffer/comeback_package/comeback/R/",
                     "sysdata.rda", sep="")

load(sysdata_path)

# Get the 450K data
I450K_Manifest<- init_data$I450K_Manifest
str(I450K_Manifest)
head(I450K_Manifest)

# Get the EPIC data
EPIC_Manifest<- init_data$EPIC_Manifest
str(EPIC_Manifest)
head(EPIC_Manifest)

save(I450K_Manifest, file = "illumina_manifests/I450K_Manifest.RData")
save(EPIC_Manifest, file = "illumina_manifests/EPIC_Manifest.RData")

# Clean up 
rm(sysdata_path, init_data, I450K_Manifest, EPIC_Manifest)
