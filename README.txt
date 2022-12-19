This repository contains code to reproduce the results of

Tuel and Martius, On the persistence of warm and cold spells in the Northern Hemisphere extratropics:
  regionalisation, synoptic-scale dynamics, and temperature budget. Submitted to Earth System Dynamics.

Contents:

  - Normalise_temperature.R normalises the daily-mean temperature data with a 30-day, 7-year moving window.
  - clustering_model_t2m.R implements the quantile regression with geopotential height as covariate, and
    the regionalisation with the Partitioning Around Medoids (PAM) algorithm.
  - The netCDF file "regionalisation_results_21day_1deg_NH_PAM_best.nc" contains the results of the 
    regionalisation (for DJF/JJA and cold/warm spells).

ERA5 data (2-meter temperature and geopotential height at 500 hPa) can be obtained from the
  Copernicus Climate Data Store (CDS) at https://doi.org/10.24381/cds.adbb2d47.
