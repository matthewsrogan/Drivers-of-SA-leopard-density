#####
#ROGAN ET AL R SCRIPT
#####

#Code to model leopard density at 27 study sites using multi-session spatial capture-recapture 
#models with inhomogenous density.

##### STEP 1: PREPARE DATA #####

### PREP WORKSPACE
require(tidyverse)
require(secr)
require(rgdal)

setwd(".") #set working director as repository

#READ IN DATA

#create list of trap objects. Coordinates are in a custom equidistant conic projection
trap.files <- list.files("./data", pattern = "ms.chtrap")#set of txt files, one for each session (i.e., survey)

#convert trap files to list of trap objects
trap.objs <- lapply(trap.files, function(file){ 
  read.traps(file, 
             detector = "count", 
             binary.usage = F)}) #usage is measured in number of trap days

#Read in captures
captures <- read.table("./data/ms.chcapt.txt", #text file with the session, individual ID, station, occasion, and partially observed sex of each capture
                       header = F, 
                       col.names = c("Session", "Individual", "Occasion", "StationID", "Sex"))
captures$Sex <- as.factor(captures$Sex) #convert to factor: 1 = female, 2 = male

#Read in session covariates
Sesh_covs <- read_csv("./data/Session_covariates.csv",
                      col_types = cols(Biome = col_character())) #read in table with session-level covariates for each session (i.e., survey)
Sesh_covs_scl <- Sesh_covs %>% 
  mutate_at(vars(log_PA:Eff_size), ~scale(.x)[,1])

### MAKE CAPTURE HISTORY OBJECT
#combine captures and traps into capthist
ms.ch <- make.capthist(captures, 
                       trap.objs, 
                       fmt = "trapID", 
                       noccasions = 1, 
                       covnames = "Sex",
                       bysession = T, 
                       sortrows = T, 
                       noncapt = "9999") #used for sites with 0 captures
verify(ms.ch)

### MAKE MASKS

#Read in bounding polygon for masks.
poly <- rgdal::readOGR("./GIS_layers/bounding_polygon.shp") #GIS layer of study region with "urban" and "water" land uses excluded
sp::proj4string(poly) #projection info of equidistant conic

#make masks for each trap object within 15 km buffer around stations with 500m cells
mask15km_list <- map(trap.objs, ~make.mask(.x, 
                                           buffer = 15000, 
                                           spacing = 500, 
                                           type = "trapbuffer", 
                                           poly = poly, 
                                           keep.poly = F,
                                           check.poly = F))#If true, will return warning for one station on edge of Loskop Dam

##### END OF STEP 1 #####

