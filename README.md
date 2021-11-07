# Drivers-of-SA-leopard-density
Code and data from Rogan et al. Troubled spots: Human impacts constrain the density of an apex predator inside protected areas.

This repository contains data and a series of scripts to fit multi-session spatial capture-recapture models with inhomogenous density 
using package 'secr' (see https://cran.r-project.org/web/packages/secr/vignettes/secr-datainput.pdf for more information).

The models are fit to camera-trap survey data of leopard populations in 27 South African protected areas.

This repository replicates the code of Rogan et al. Appendix S3 but splits the workflow into a series of distinct scripts to improve readability.
The analytical workflow is divided into four steps: 1) prepping data, 2) Fitting a preliminary model, 3) Extracting covariate values, and 4) Fitting and evaluating models. 
A fifth script, 'source_funs.R' contains three custom functions used in the workflow, one of which pertains to creating starting values for secr models and two that pertain to
calculating distance-weighted mean values from rasters.

The data folder consists of a multi-session capture history file (ms.chcapt.txt), 27 trap files (ms.chtrap[SITE].txt), 
and three spatial layers: a polygon shapefile denoting potential leopard habitat used to define the model statespace,
a polygon shapefile denoting non-protected terrestrial areas of Southern Africa,
and a raster of the 2009 human footprint index acquired from Venter et al. 2016 (https://doi.org/10.1038/sdata.2016.67)
that has been clipped to southern Africa and reprojected to the custom equidistant conic projection used in the analysis.

Data were collected by Panthera's Southern Africa Leopard Monitoring Program and processed by Panthera's Integrated Data Systems team. 
Contact Matt Rogan (mrogan@panthera.org) for information about the data presented here and IDS@panthera.org for information about unprocessed camera-trap survey data.

