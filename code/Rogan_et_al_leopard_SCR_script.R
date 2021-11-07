#####
###ROGAN ET AL R SCRIPT###
#Model leopard density at 27 study sites using multi-session spatial capture-recapture models with inhomogenous density


###MAJOR STEPS: 
#1) prep workspace and data [Starts line 15]
#2) fit session model [Starts line 64] 
#3) extract covariate values within constrained state space [Starts line 137]
#4) fit univariate models to derive starting values [Starts line 316] 
#5) fit models and examine results [Starts line 474]


#####
###STEP 1: PREP WORKSPACE
library(tidyverse)
library(secr)

setwd("./Archived_data")#navigate to folder with archived data

#READ IN DATA

#create list of trap objects. Coordinates are in a custom equidistant conic projection
trap.files <- list.files(".", pattern = "ms.chtrap")#set of txt files, one for each session (i.e., survey)

#convert trap files to list of trap objects
trap.objs <- lapply(trap.files, function(file){ 
  read.traps(file, 
             detector = "count", 
             binary.usage = F)}) #usage is measured in number of trap days

#Read in captures
captures <- read.table("ms.chcapt.txt", #text file with the session, individual ID, station, occasion, and partially observed sex of each capture
                       header = F, 
                       col.names = c("Session", "Individual", "Occasion", "StationID", "Sex"))
captures$Sex <- as.factor(captures$Sex) #convert to factor: 1 = female, 2 = male

#Read in session covariates
Sesh_covs <- read_csv("Session_covariates.csv",
                      col_types = cols(Biome = col_character())) #read in table with session-level covariates for each session (i.e., survey)
Sesh_covs_scl <- Sesh_covs %>% 
  mutate_at(vars(log_PA:Eff_size), ~scale(.x)[,1])


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

#Read in bounding polygon for masks.
poly <- rgdal::readOGR("./GIS_layers/bounding_polygon.shp") #GIS layer of study region with "urban" and "water" land uses excluded
sp::proj4string(poly) #projection info of equidistant conic

###END OF STEP 1
#####

#####
###STEP 2: FIT SESSION MODEL WITH FULL STATE SPACE
#make masks for each trap object within 15 km buffer around stations with 500m cells
mask15km_list <- map(trap.objs, ~make.mask(.x, 
                                           buffer = 15000, 
                                           spacing = 500, 
                                           type = "trapbuffer", 
                                           poly = poly, 
                                           keep.poly = F,
                                           check.poly = F))#If true, will return warning for one station on edge of Loskop Dam


#starting values for D~session model derived from Rogan et al 2019
D0 <- log(6.47*1e-4)#Density intercept (i.e., first session) based on Atherstone 2016 density estimate from Rogan et al. 2019

##starting values for the effect of each site relative to the intercept
#create a function to approximate a beta coefficent based on D per 100 sq. km. and D intercept in log-transformed units per ha)
beta_fn <- function(D, d0) return(log(D*1e-4) - d0) 

#use 0.5 leopards per 100 sq. km. for sites with one or no captures, 1 for sites with a few individuals but too few to model effectively, else the estimate from Rogan et al 2019
D_est <- c(0.5, 1.7, 0.5, 9.5, 3.6, 8.4, 6.8, 6.9, 4.3, 9.5, 2.7, 
           9.8, 3.7, 1, 8.7, 8.2, 1.7, 0.5, 7.5, 7.1, 6.9, 6.2,
           3.1, 3, 3.7, 1)

##starting values for the detection function
#approximate medians for each biome from Rogan et al 2019
det_est <- c(log(0.04), #first lambda0 (log transformed)
             0, 
             log(0.03) - log(0.04), 
             #then sigma (log transformed)
             log(1800), 
             log(2200) - log(1800), 
             log(2600) - log(1800))

#combine into one vector of starting values
Beta0 <- c(D0, beta_fn(D_est, D0), det_est)

###FIT MODEL
#Note, model fitting may take 24 hours or longer. Adjust number of cores (argument "ncores") as needed.
session_mod <- secr.fit(capthist = ms.ch, 
                        model = list(D ~ session, #session as categorical predictor
                                     lambda0 ~ as.factor(Biome), 
                                     sigma ~ as.factor(Biome)), 
                        mask = mask15km_list, 
                        detectfn = "HHN", #Hazard half-normal
                        binomN = 0, #poisson process
                        start = Beta0, 
                        sessioncov = as.data.frame(Sesh_covs_scl),
                        details = list(fastproximity = FALSE), 
                        method = "Nelder-Mead", 
                        trace = F, 
                        ncores = 6, 
                        control = list(maxit = 19999))

#check convergence
session_mod$fit$convergence

#sigma estimates:
exp(coef(session_mod)[31, 1]) #Biome 1 - Forest
exp(coef(session_mod)[31, 1] + coef(session_mod)[32, 1]) #Biome 10 - grasslands
exp(coef(session_mod)[31, 1] + coef(session_mod)[33, 1]) #Biome 7 - savannahs

#calculate sigma_hat
round(mean(exp(coef(session_mod)[31, 1]),
     exp(coef(session_mod)[31, 1] + coef(session_mod)[32, 1]),
     exp(coef(session_mod)[31, 1] + coef(session_mod)[33, 1])))

#save progress
save.image("Session_mod_image.rdata")

###END OF STEP 2
#####

#####
###STEP 3: EXTRACT COVARIATE VALUES
###Create sigma masks
masks_sigma <- lapply(trap.objs, function(traps){
  msk <- make.mask(traps, 
                   buffer = 2.447*sigma_hat, 
                   spacing = 500, 
                   type = "trapbuffer",
                   poly = bb_poly, 
                   keep.poly = F, 
                   check.poly = F)
  return(msk)
})

names(masks_sigma) <- Sesh_covs$SessionID

#run session model with constrained state space to compare with 15km mask
sigma_session_mod <- secr.fit(ms.ch, 
                              model = list(D ~ session, lambda0 ~ as.factor(Biome), sigma ~ as.factor(Biome)),
                              mask = masks_sigma_covs, 
                              detectfn = "HHN", 
                              binomN = 0, 
                              start = session_mod, 
                              sessioncov = as.data.frame(Sesh_covs_scl), 
                              details = list(fastproximity = FALSE), 
                              method = "Nelder-Mead", 
                              trace = F, 
                              ncores = 6, 
                              control = list(maxit = 19999))

#compare parameter values to session_mod
par_comp <- full_join(as_tibble(coef(session_mod), rownames = "par"), 
                      as_tibble(coef(sigma_session_mod), rownames = "par"),
                      by = "par", 
                      suffix = c("_full", "_sigma")) %>%
  mutate_at(vars(beta_full:ucl_sigma), ~{round(.x, 2)}) %>%
  mutate(out_full = paste0(beta_full, " (", SE.beta_full, ")"), 
         out_sigma = paste0(beta_sigma, " (", SE.beta_sigma, ")"),
         prop_diff = round((beta_sigma - beta_full)/beta_full, 2))

par_comp %>% dplyr::select(par, out_full, out_sigma, prop_diff) %>% 
  write_csv("full vs sigma state space comparison.csv")

###extract raster values
#Load spatial packages
library(raster)
library(sf)
library(nngeo)

##create functions for calculating distance weighted means
#function for calculating weights of each raster pixel
kernel_weight <- function(distance, sigma){
  exp(-distance^2/(2*sigma^2))/sum(exp(-distance^2/(2*sigma^2)))
}

#function to calculate raster values within 2.447*sigma of each point and weight them using kernel_weight()
dw_covs <- function(mask_pts, r, sigma, proj){#mask_pts is spatial object of habitat cell centroids, 
  #r is the predictor variable raster, sigma is the distance weighting scale parameter, proj is the projection of the masks
  
  #convert raster to pts, clip extent, and transform to the projection of the habitat mask
  if(str_detect(names(r), "A20") == T) {r[r<0] = NA} #set NA values in NDVI layers
  r_pts <- rasterToPoints(crop(r, st_combine(mask_pts) %>% 
                                 st_convex_hull() %>% 
                                 st_buffer(dist = 2.447*sigma) %>% 
                                 st_transform(st_crs(r)) %>%  
                                 as_Spatial()), spatial = T) %>% 
    st_as_sf() %>% 
    st_transform(proj)
  
  r_pts <- mutate(r_pts, x = st_coordinates(r_pts)[,1], y = st_coordinates(r_pts)[,2]) %>% 
    st_set_geometry(NULL)
  
  st_geometry(mask_pts) <- NULL
  Var <- as.numeric(NA)
  
  #for each hab pt, create dwm. 
  for (j in 1:nrow(mask_pts)){
    pt <- mask_pts[j, c("x", "y")]
    dist_vect <- sqrt((pt$x - r_pts$x)^2 + (pt$y - r_pts$y)^2)
    Var[j] = sum(r_pts[which(dist_vect < 2.447*sigma), 1] *
                   kernel_weight(dist_vect[dist_vect < (2.447 * sigma)], sigma))
  }
  return(Var)
}

###Read in covariate layers
#NDVI - folder with layers for the entire study period
#extract dates of each layer
ndvi_fol <- "C:/Users/Matthew Rogan/Dropbox/PhD/GIS/Env_Variables/NDVI_WGS84"
ndvi_lyrs <- str_sub(list.files(path=ndvi_fol, pattern = ".tif"), 2, 8) %>% unique()
ndvi_days<-as.integer(ndvi_lyrs) #set of julian days ndvi layers were collected
ndvi_dates<-as.POSIXlt(ndvi_lyrs, format="%Y%j")

#HFI
#read in hfp, already transformed to cec when previously clipped to southern African
hfi <- raster("./GIS_layers/hfp_SA_cec.tif")
hfi[hfi == 128] = NA

#Edge_dist
#read in unprotected areas and transform to cec
edge <- st_read("./GIS_layers/unprotected.shp") %>% 
  st_transform(cec)

#dist2river
#read in SA rivers
rivers <- st_read("./GIS_layers/Rivers_SA_WGS84.shp") %>%
  filter(CLASS == "Perennial") %>% 
  st_combine() %>% 
  st_transform(cec)

###extract covariate values for each survey
msk.sigma.covs <- list() #list of spatial points of sigma masks with covariate values

for(i in 1:nrow(Sesh_covs)){
  sesh <- Sesh_covs$SessionID[i]
  midpoint <- Sesh_covs$Midpoint[i] #use midpoint to identify appropriate NDVI layer
  
  #Find the NDVI lyr closest to midpoint and load it to object ndvi
  ndvi_day <- ndvi_days[which.min(abs(ndvi_days-(as.integer(format(ymd(midpoint), "%Y%j")) - 8)))]
  ndvi <- raster(paste0(ndvi_fol, "/A", ndvi_day, ".tif"))/10000
  ndvi[ndvi < 0] = NA #exclude water bodies
  
  #convert masts to spatial points
  msk_pts <- masks_sigma[[i]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = cec, remove = F)
  
  #sample covariates
  msk.sigma.covs[[i]] <- msk_pts %>% mutate(NDVI = dw_covs(msk_pts, ndvi, sigma_hat, cec),
                                    HFI = dw_covs(msk_pts, hfi, sigma_hat, cec),
                                    Edge = simplify(st_nn(msk_pts, edge, k=1, returnDist = T)[[2]]),
                                    River = simplify(st_nn(msk_pts, rivers, k=1, returnDist = T)[[2]]),
                                    Edge_log = log(replace(Edge, Edge < 1, 1)),
                                    River_log = log(replace(River, River < 1, 1)),
                                    Prtcd = if_else(Edge == 0, 0, 1),
                                    SessionID = Sesh_covs$SessionID,) 
}

###scale covariates
#combine list of mask pts into data frame with continuous predictors scaled to have mean 0 and std dev 1
hab_covs_sigma <- do.call(bind_rows, 
                          lapply(msk.sigma.covs2, function(x) {
                            x %>% st_set_geometry(NULL)
                            })) %>% 
  dplyr::select(-x, -y, -Edge, -River) %>%
  mutate_at(vars(c("NDVI", "HFI", "River_log")), function(x) scale(x)[,1]) %>%
  rename(River = River_log)

#scale the edge variable using only values inside protected areas (ie Prtcd == 1)
log_dists <- filter(hab_covs_sigma, Prtcd == 1) %>% pull(Edge_log)
hab_covs_sigma <- hab_covs_sigma %>%
  mutate(Edge = if_else(Prtcd == 0, 0, (Edge_log - mean(log_dists))/sd(log_dists)))

###check correlations
cor_df <- dplyr::select(Sesh_covs_scl, SessionID, Prey, Domestic, log_PA, Mgmt, Eff_size, Edge_core) %>%
  full_join(hab_covs_sigma, by = "SessionID") %>% dplyr::select(-SessionID)
cor_tbl <- round(cor(cor_df), digits = 2)
write_csv(cor_tbl, "covariate_correlations.csv")

###add covariates to list of mask objects
masks_sigma_covs <- lapply(seq_along(masks_sigma), function(i){
  sesh = names(masks_sigma)[[i]]
  msk = masks_sigma[[i]]
  covariates(msk) <- filter(hab_covs_sigma, SessionID == sesh) %>%
    as.data.frame()
  return(msk)
})

names(masks_sigma_covs) <- names(masks_sigma)

#save progress as image and just objects needed for subsequent steps
save.image("covariate_processing_image.rdata")
save(ms.ch, masks_sigma_covs, Sesh_covs_scl, session_mod, sigma_hat, 
     file = "model_fitting_inputs.rdata")

#clear workspace
rm(list = ls())

###END OF STEP 3
#####

#####
###STEP 4: DERIVE STARTING VALUES FOR MODELS
#extract starting values for univariate models from session_mod
#calculate mean density (on log scale) among sessions to use as intercept for density
d0 <- coef(session_mod)$beta[1]
d_betas <- c(d0, coef(session_mod)$beta[2:27] + d0)
d_hat <- mean(d_betas) #geometric mean density across sites

#extract detection starting values
det_starts <- c(coef(session_mod)[c("lambda0", 
                                  "lambda0.as.factor(Biome)10", 
                                  "lambda0.as.factor(Biome)7"), 1],
                      log(sigma_hat))

univar_starts = round(c(d_hat, #density intercept
                        0, #starting value for predictor effect
                        det_starts), 
                      digits = 3)

#create second set of start values for NDVI which has linear and quadratic effect and Edge:Prtcd which has Prtcd main effect
univar_starts2 = round(c(d_hat, 0, 0, det_starts), 3)

###FIT UNIVARIATE MODELS


#Prey
univar_mods$prey <- secr.fit(ms.ch, 
                             model = list(D ~ Prey, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#NDVI
univar_mods$ndvi <- secr.fit(ms.ch, 
                             model = list(D ~ NDVI, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts2, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#River distance
univar_mods$river <- secr.fit(ms.ch, 
                             model = list(D ~ River, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#PA size
univar_mods$prey <- secr.fit(ms.ch, 
                             model = list(D ~ log_PA, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#Management type
univar_mods$mgmt <- secr.fit(ms.ch, 
                             model = list(D ~ mgmt, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#Domestic
univar_mods$dom <- secr.fit(ms.ch, 
                             model = list(D ~ Domestic, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#HFI
univar_mods$hfi <- secr.fit(ms.ch, 
                             model = list(D ~ HFI, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#Protected
univar_mods$prtcd <- secr.fit(ms.ch, 
                             model = list(D ~ Prtcd, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))


#Edge distance - include Prtcd main effect as well.
univar_mods$edge <- secr.fit(ms.ch, 
                             model = list(D ~ Prtcd + Prtcd:Edge, lambda0 ~ as.factor(Biome), sigma ~ 1), 
                             mask = masks_sigma_covs,
                             detectfn = "HHN", 
                             sessioncov = as.data.frame(Sesh_covs_scl), 
                             trace = T, 
                             start = univar_starts2, 
                             binomN = 0, 
                             ncores = 6, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#save univariate models
write_rds(univar_mod, "Univariate_models.rds")

###END STEP 4
#####

#####
###STEP 5: FIT MODELS
#specify competing models
models <- list()
models$m1 <- list(D ~ Prey + NDVI + I(NDVI^2) + River, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m2 <- list(D ~ NDVI + I(NDVI^2) + River, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m3 <- list(D ~ Prtcd, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m4 <- list(D ~ log_PA + Mgmt + Prtcd, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m5 <- list(D ~ log_PA, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m6 <- list(D ~ HFI + Domestic + Prtcd + Prtcd:Edge, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m7 <- list(D ~ Eff_size, lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m8 <- list(D ~ log_PA + Mgmt + Domestic + HFI + Prtcd + Prtcd:Edge, 
                  lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m9 <- list(D ~ Domestic + HFI + Prtcd + Prtcd:Edge + Mgmt:HFI + Mgmt:Prtcd + Mgmt:Prtcd:Edge,
                  lambda0 ~ as.factor(Biome):h2, sigma ~ h2)
models$m10 <- list(D ~ log_PA + Mgmt + Domestic + NDVI + I(NDVI^2) + River, 
                   lambda0 ~ as.factor(Biome):h2, sigma ~ h2)

#Extract starting values for detecton parameters from mean of univariate models
pars.list <- lapply(univar_mods, function(mod) return(mod$fit$par))
d0 <- mean(sapply(pars.list, first))

lam0 <- mean(sapply(pars.list, function(x) x["lambda0"]))
lam0.10 <- mean(sapply(pars.list, function(x) x["lambda0.as.factor(Biome)10"]))
lam0.7 <- mean(sapply(pars.list, function(x) x["lambda0.as.factor(Biome)7"]))
sigma <- mean(sapply(pars.list, function(x) x["sigma"]))

#create sex specific starting values
sigma.fem <- sigma - 0.2
sigma.male <- 0.4 
pmix <- logit(0.4)
det_starts <- c(lam0, 0, lam0.10, lam0.7, 0, lam0.10, lam0.7, sigma.fem, sigma.male, pmix)

#Specify start values for each model
start_list <- list(m1 = c(d0, univar_mods$prey$fit$par[2], univar_mods$ndvi$fit$par[2:3],
                          univar_mods$river$fit$par[2], det_starts),
                   m2 = c(d0, univar_mods$ndvi$fit$par[2:3],
                          univar_mods$river$fit$par[2], det_starts),
                   m3 = c(d0, univar_mods$prtcd$fit$par[2], det_starts),
                   m4 = c(d0, univar_mods$size$fit$par[2], univar_mods$mgmt$fit$par[2],
                          univar_mods$prtcd$fit$par[2], det_starts),
                   m5 = c(d0, univar_mods$size$fit$par[2], det_starts),
                   m6 = c(d0, univar_mods$hfi$fit$par[2], univar_mods$dom$fit$par[2],
                          univar_mods$edge$fit$par[2:3], det_starts),
                   m7 = c(d0, 0, det_starts),
                   m8 = c(d0, univar_mods$size$fit$par[2], univar_mods$mgmt$fit$par[2],
                          univar_mods$dom$fit$par[2], univar_mods$hfi$fit$par[2],
                          univar_mods$edge$fit$par[2:3], det_starts),
                   m9 = c(d0, univar_mods$dom$fit$par[2], univar_mods$hfi$fit$par[2],
                          univar_mods$edge$fit$par[2:3], rep(0, 3), det_starts),
                   m10 = c(d0, univar_mods$size$fit$par[2], univar_mods$mgmt$fit$par[2],
                           univar_mods$dom$fit$par[2], univar_mods$ndvi$fit$par[2:3],
                           univar_mods$river$fit$par[2], det_starts))

###FIT MODELS
fit_models <- map2(models, start_list, function(mod, strt){
  secr.fit(ms.ch, 
           model = mod, 
           mask = masks_sigma_covs, 
           ncores = 6, 
           hcov = "Sex",
           detectfn = "HHN", 
           sessioncov = as.data.frame(Sesh_covs_scl), 
           trace = F, 
           details = list(fastproximity = FALSE),
           method = "Nelder-Mead", 
           control = list(maxit = 19999), 
           start = strt)
}) 

#refit with updated start values
mods <- map(models, function(mod){
  secr.fit(ms.ch,
           model = mod$model, 
           mask = masks_sigma_covs, 
           ncores = 6, 
           hcov = "Sex",
           detectfn = "HHN", 
           sessioncov = as.data.frame(Sesh_covs_scl), 
           trace = F, 
           details = list(fastproximity = FALSE),
           method = "Nelder-Mead", 
           control = list(maxit = 19999), 
           start = mod)
})

###Look at results
secrmods <- secrlist(mods)

#check convergence
lapply(mods, function(mod) print(mod$fit$convergence))

#check compatibility
AICcompatible(mods)

#check support
AIC(mods, criterion = "AICc", dmax = 7)

#save AIC table
AIC_tbl <- as.tibble(AIC(mods, criterion = "AICc", dmax = 7), rownames = "Hyp") %>% 
  mutate(model = str_sub(model, 1, -46),
         dAICc = round(dAICc, 2), Weight = round(AICcwt, 2)) %>% 
  dplyr::select(-detectfn, -logLik, -AICc) %>%
  write_csv("AICc_table.csv")

#Extract predictor effects
coef_list <- list()
for(m in names(mods)){
  coef_list[[m]] <- as_tibble(coef(mods[[m]]), rownames = "par") %>% 
    mutate(Hyp = m)
}

par_vals <- do.call(bind_rows, coef_list) %>%
  mutate_if(is.numeric, ~round(.x, 2)) %>%
  mutate(estimate = paste0(beta, " (", SE.beta, ")")) %>%
  dplyr::select(Hyp, par, estimate) %>%
  group_by(Hyp, par) %>%
  spread(key = par, value = estimate, fill = NA) %>%
  left_join(dplyr::select(AIC_tbl, Hyp, dAICc, Weight), by = "Hyp") %>%
  arrange(dAICc) %>%
  rename_at(vars(D.Domestic:D.River), ~{str_sub(., 3, -1)}) %>%
  mutate_at(vars(Domestic:River), ~{replace_na(., "")}) 

#write predictor coefficients to table
dplyr::select(par_vals, Hyp:River, dAICc, Weight) %>%
  write_csv("Predictor_effects_table.csv")

#Extract detection parameters
var_names <- c("lam0", "lam0_1fem", "lam0_10fem", "lam0_7fem", "lam0_1male", "lam0_10male",
               "lam0_7male", "sigma", "sigma_male", "pmix_male")
det_list <- list()
for(m in names(mods)){
  det_list[[m]] <- as_tibble(coef(mods[[m]])) %>% slice((nrow(.)-9):nrow(.)) %>%
    mutate(Hyp = m, par = var_names)
}

det_pars <- do.call(bind_rows, det_list) %>%
  dplyr::select(Hyp, par, beta, SE.beta) %>%
  group_by(Hyp) %>% 
  mutate(val = paste0(beta, " (", SE.beta,))
spread(key = par, value = beta) %>%
  rowwise() %>% 
  mutate(lam0_1fem = exp(lam0 + lam0_1fem),
         lam0_10fem = exp(lam0 + lam0_10fem),
         lam0_7fem = exp(lam0 + lam0_7fem),
         lam0_1male = exp(lam0 + lam0_1male),
         lam0_10male = exp(lam0 + lam0_10male),
         lam0_7male = exp(lam0 + lam0_7male),
         sigma_fem = exp(sigma), sigma_male = exp(sigma + sigma_male),
         pmix_male = invlogit(pmix_male)) %>%
  dplyr::select(-lam0, -sigma) %>%
  mutate_at(vars(lam0_10fem:lam0_7male), ~{round(.x, 3)}) %>%
  mutate(pmix_male = round(pmix_male, 2)) %>%
  mutate_at(vars(sigma_male:sigma_fem), ~{round(.x, 0)}) %>%
  left_join(dplyr::select(AIC_tbl, Hyp, dAICc, Weight)) %>%
  arrange(dAICc) %>% write_csv("Detection_parameters_table.csv")

###END STEP 5
#####


