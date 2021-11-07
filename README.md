# Drivers-of-SA-leopard-density
Code and data from Rogan et al. Troubled spots: Human impacts constrain the density of an apex predator inside protected areas.

This repository contains data and a series of scripts to fit multi-session spatial capture-recapture models with inhomogenous density using package 'secr' (Efford 2019).
The models are fit to camera-trap survey data of leopard populations in South African protected areas.

This repositor replicates the code of Rogan et al. Appendix S3 but splits the workflow into a series of distinct scripts to improve readability.
The analytical workflow is divided into four steps: 1) prepping data, 2) Fitting a preliminary model, 3) Extracting covariate values, and 4) Fitting and evaluating models. 
A fifth script, 'source_funs.R' contains three custom functions used in the workflow, one of which pertains to creating starting values for secr models and two that pertain to
calculating distance-weighted mean values from rasters.
