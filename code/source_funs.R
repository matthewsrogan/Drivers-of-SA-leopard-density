#####
### Create function for specifying starting values of beta coefficients
#####

#create a function to approximate a beta coefficent based on D per 100 sq. km. 
#and D intercept in log-transformed units per ha)

beta_fn <- function(D, d0) return(log(D*1e-4) - d0) 

#########################################################################################

#####
### Create functions for calculating distance-weighted mean predictor values from rasters
#####

#####
### kernel_weight
#####

#function for calculating weights of each raster pixel
kernel_weight <- function(distance, sigma){
  exp(-distance^2/(2*sigma^2))/sum(exp(-distance^2/(2*sigma^2)))
}


#####
### dw_covs
#####

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