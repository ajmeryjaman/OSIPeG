## This R-package constructs valid inference for penalized G-estimation using a decorrelated score approach.

This repository contains the following folders:

Folder | Description
--- | ---
R | Contains the source codes
man | Contains the documentation of each function used

The R folder contains the following files:

File | Description
--- | ---
[osipeg_confint.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/osipeg_confint.R) | Contains the main function osipeg_confint() which implements our method
[optimalW_dantzig.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/optimalW_dantzig.R) | Contains the function optimalW_dantzig() that calculates the sparse weight vector via the Dantzig selector 
[optimalW_lasso.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/optimalW_lasso.R) | Contains the function optimalW_lasso() that calculates the sparse weight vector via the LASSO 

Please see the example given in [osipeg_confint.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/osipeg_confint.R) to generate a longitudinal data set and implement our method. Or, do the following:

#### R commands for installing and using our package

library(devtools) # if already installed, otherwise need to install it first

install_github("ajmeryjaman/OSIPeg")

library(OSIPeg)

?osipeg_confint # To follow the example

