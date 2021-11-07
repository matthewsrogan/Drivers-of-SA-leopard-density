#####
#ROGAN ET AL R SCRIPT
#####

#Code to model leopard density at 27 study sites using multi-session spatial capture-recapture 
#models with inhomogenous density.

#run script './code/data_prep.R' prior to running this script.

#Source custom functions
source(".code/source_funs.R")

##### STEP 2: FIT PRELIMINARY MODEL #####

### SPECIFY STARTING VALUES

#starting values for D~session model derived from Rogan et al 2019
#Density intercept (i.e., first session) based on Atherstone 2016 density estimate from Rogan et al. 2019
D0 <- log(6.47*1e-4) # convert denstiy per 11 sq km to density per hectare with log link

# starting values for the effect of each site relative to the intercept
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

#combine into one vector of starting values using beta_fn() from source_funs.R
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
save.image("./outputs/Session_mod_image.rdata")

##### END OF STEP 2 #####