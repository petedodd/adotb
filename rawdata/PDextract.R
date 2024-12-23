## data on age-specific biths for pregnancy analysis:
## restricting to relevant piece of large (~3Gb) WPP2024
library(here)
library(data.table)

## PD <- fread("~/Downloads/WPP2024_Fertility_by_Age1.csv.gz") # zipped full data
PD <- fread("~/Downloads/WPP2024_Fertility_by_Age5.csv.gz")
print(object.size(PD), units = "auto") # ~3Gb for 1 year, 1Gb for 5 yr
PD <- PD[Time == 2019] # only need 2019


PD[,unique(ISO3_code)]
PD <- PD[ISO3_code != ""] # only need specific countries
print(object.size(PD), units = "auto") # better


save(PD, file = here("rawdata/PD.Rdata"))

