#####
#ROGAN ET AL R SCRIPT
#####

#Code to model leopard density at 27 study sites using multi-session spatial capture-recapture 
#models with inhomogenous density.

##### STEP 4: MODEL FITTING #####
require(tidyverse)
require(secr)

setwd(".") #set working director as repository

#set number of available cores for fitting models in parallel
n_cores = 6

### LOAD INPUTS
load("./outputs/model_fitting_inputs.rdata")

### DERIVE STARTING VALUES FOR MODELS
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
                             ncores = n_cores, 
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
                             ncores = n_cores, 
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
                              ncores = n_cores, 
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
                             ncores = n_cores, 
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
                             ncores = n_cores, 
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
                            ncores = n_cores, 
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
                            ncores = n_cores, 
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
                              ncores = n_cores, 
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
                             ncores = n_cores, 
                             details = list(fastproximity = FALSE),
                             method = "Nelder-Mead", 
                             control = list(maxit = 19999))

#save univariate models
write_rds(univar_mod, "./outputs/Univariate_models.rds")

#####
### SPECIFY MODELS
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

### EXTRACT STARTING VALUES
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
           ncores = n_cores, 
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
           ncores = n_cores, 
           hcov = "Sex",
           detectfn = "HHN", 
           sessioncov = as.data.frame(Sesh_covs_scl), 
           trace = F, 
           details = list(fastproximity = FALSE),
           method = "Nelder-Mead", 
           control = list(maxit = 19999), 
           start = mod)
})

### EVALUATE RESULTS
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
  write_csv("./outputs/AICc_table.csv")

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
  write_csv("./outputs/Predictor_effects_table.csv")

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
  arrange(dAICc) %>% write_csv("./outputs/Detection_parameters_table.csv")

##### END OF STEP 4 #####

