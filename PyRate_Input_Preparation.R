library(dplyr)
library(googlesheets4)
source("/Applications/PyRate/pyrate_utilities.R")

# read in the ASH Taxonomy List
gs4_deauth()
prinput <- read_sheet("https://docs.google.com/spreadsheets/d/1rZQaO_Qhdmfo5Ck2LkPD9UrBGibjukFwVHA9mM3-PoE/edit#gid=0")

# filter it down to just Varanus and just species-level fossils
spinput <- filter(prinput, accepted_rank == "species" & genus == "Varanus")

# remove duplicate records (same species at the same site)
redinput <- spinput[,c("spp_extant", "genus", "accepted_name", 
                       "accepted_rank", "max_ma", "min_ma", "fossil_site", "source")]
unique.input <- distinct(redinput)

# write this new input data to a file
write.csv(unique.input, file="~/Documents/GitHub/MonitorDiversification/Data/FossilRecords_Varanus.csv", row.names = FALSE)

# designate the extant taxa
exinput <- filter(unique.input, spp_extant == "extant")
extants <- unique(exinput$accepted_name)

# create PyRate input files
extract.ages.pbdb(file="~/Documents/GitHub/MonitorDiversification/Data/FossilRecords_Varanus.csv",
                  extant_species = extants, replicates = 100) # replicates = 10
