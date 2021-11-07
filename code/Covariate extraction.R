#####
#ROGAN ET AL R SCRIPT
#####

#Code to model leopard density at 27 study sites using multi-session spatial capture-recapture 
#models with inhomogenous density.

#run scripts './code/data_prep.R' and './code/Preliminary_model.R' prior to running this script.

#See Rogan et al. Appendix S1 for details of predictor covariates

##### STEP 3: EXTRACT COVARIATE VALUES #####

###CREATE SIGMA MASKS
#makes constrained masks based on a buffer of 2.447*sigma_hat
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

### FIT SESSION MODEL
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
  write_csv("./outputs/full vs sigma state space comparison.csv")

###EXTRACT RASTER VALUES
#Load spatial packages
library(raster)
library(sf)
library(nngeo)

###Read in covariate layers
#get equidistant conic projecting info from the poly object
cec <- crs(poly)

#NDVI - folder with layers for the entire study period
#NB: NDVI layers should be sourced directly from a NASA remote sensing data repository

#extract dates of each layer
ndvi_fol <- "" #specify path to folder with NDVI layers
ndvi_lyrs <- str_sub(list.files(path=ndvi_fol, pattern = ".tif"), 2, 8) %>% unique()
ndvi_days<-as.integer(ndvi_lyrs) #set of julian days ndvi layers were collected
ndvi_dates<-as.POSIXlt(ndvi_lyrs, format="%Y%j")

#HFI
#read in hfp, already transformed to cec when previously clipped to southern African
hfi <- raster("./data/GIS_layers/hfp_SA_cec.tif")

# set no data value to NA
hfi[hfi == 128] = NA

#Edge_dist
#read in unprotected areas and transform to cec
edge <- st_read("./GIS_layers/unprotected.shp") %>% 
  st_transform(cec)

#dist2river
#read in SA rivers
#NB: Rivers layer should be sourced directly from http://www.dwa.gov.za/iwqs/gis_data/river/rivs500k.html
rivers <- st_read("./data/GIS_layers/Rivers_SA_WGS84.shp") %>%
  filter(CLASS == "Perennial") %>% 
  st_combine() %>% 
  st_transform(cec)

###extract covariate values for each survey
msk.sigma.covs <- list() #empty list to populate with spatial points of sigma masks with covariate values

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
save.image("./outputs/covariate_processing_image.rdata")

save(ms.ch, masks_sigma_covs, Sesh_covs_scl, session_mod, sigma_hat, 
     file = "./outputs/model_fitting_inputs.rdata")

### clear workspace
rm(list = ls())

##### END OF STEP 3 #####
