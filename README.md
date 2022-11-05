# Analysis for Geographic Mobility and Fitness of *Streptococcus pneumoniae*
*This code was run with R version 3.6.1*
## Human Mobility Model Analysis
This analysis integrates human mobility data and genomic data to mechanistically understand pneumococcal migration in South Africa.

### Install R packages. 
```install.packages(c("raster","rgdal","data.table","doParallel","ucminf","doMC","Rcpp","RcppEigen","Rfast","abind","ggplot2","fmcmc","coda","dplyr","ape","lubridate","tmaptools","geodist","PBSmapping","stringi","maptools","tidyr","stringr","ComplexHeatmap","circlize"))```

### Process raw input data.  
Raw input data is in this folder:
```./scripts/manage_rawData/```
1) **Landscan Data** (https://landscan.ornl.gov/) <br />
*This script reads in the landscan population data and associated shapefiles from GADM (https://gadm.org/data.html) to create a datatable with all population levels ad the municipality (N=234) level population sizes.*  
Script: ```031122_LandScanMunic.R``` <br />
Output: ```./data/landscan2017/LandScan_PopulationN.RData```  <br />

 Download this folder from FigShare (https://figshare.com/s/675e41ed68ece18c5c61)<br />
  Run this code ```mkdir ./data/landscan2017``` <br />
  Unzip downloaded files and place them in ```./data/landscan2017/```<br />
  ```mkdir ./data/shapefiles```<br />
  Download South Africa shapefiles from GADM (https://gadm.org/data.html)<br />
  Unzip downloaded files and place them in ```./data/shapefiles/```<br />
  
2) **Facebook Data** (https://dataforgood.facebook.com/) <br />
*This script reads in the raw facebook mobility data from the disaster movement range maps (data downloaded from FigShare URL). Alternatively you can download the output file ```mvment_SA.provinces.RData``` from FigShare directly (https://figshare.com/s/7eb72568387c476e62f5)*<br />
Script: ```facebook_rawData.R``` <br />
Output: ```./data/landscan2017/landscan_populations.RData```  <br />
```./data/facebook/mvment_SA.provinces.RData```  <br />
### Create model input files
Scripts to generate input files location: ```./scripts/createinputs/``` <br />
Output RData files location: ```./modelinput_data/``` <br />
1) **Human Mobility** <br />
*This script creates probability of moblity matrices at the province level (9X9) and the municipality level (234X234)*  
Script: ```mobility_probability.R``` <br />
Output: Province Mobility:```./modelinput_data/cdr.mat.one.RData ``` <br />
Municipality Mobility:```./modelinput_data/cdr.mat.town.one.RData ``` <br />
2) **Municipality Population Sizes** <br />
Script: ```mobility_probability.R``` <br />
Output: Municipality Populations:```./modelinput_data/pop_municipality.2017LS.RData```  <br />
3) **Metadata to Generate Likelihood**<br />
*This creates a list of length 9 (number of GPSCs) containing 8 columns with information for each pair including ids, collection location (province numbered 1-9), time in days between pairs, collection year for each pair.*<br />
Script: ```031122_GPSC_metadata.R```<br />
Output: Metadata File: ```./modelinput_data/dat.tmp.allser.RData```<br />
Population size (province): ```./modelinput_data/pop_2019.RData: ```<br />
4) **tMRCA between pairs**<br />
*This creates a list of length 9 (number of GPSCs) with pairwise time to MRCA for each pair (half the time between pairs (years)).*<br />
Script: ```031122_makeGDistMatrix.R```<br />
Output: ```./modelinput_data/tMRCAs.RData```<br />

### Run MCMC Model<br />
*You will have to create an output directory to write the chains to and can then run the code to fit the model*<br />
  ```cd ./MCMC_model/```<br/>
 ```mkdir outputs```<br/>
 **Files**<br/>
Likelihoods Functions: ```./MCMC_model/LikelihoodFunctions/```<br />
Matrix Multiplication Cpp file: ```./MCMC_model/MatrixMultiplication.cpp```<br />
Run Model Files: ```./MCMC_model/RunModel/041122_Pneumo_MCMC_MUNIC.R```<br />
1) Open ```./MCMC_model/RunModel/041122_Pneumo_MCMC_MUNIC.R``` <br />
2) Set the number of iterations ```iters=20000```
3) Set the chain number ```chain=1```
4) Load R on your cluster environment and submit a job with ```Rscript ./MCMC_model/RunModel/041122_Pneumo_MCMC_MUNIC.R```
5) Repeat with at least 3 chains remembering to change the ```chain``` variable as this will prevent your files from overwriting each other.

### Test Model Fit <br />
*Test your Model Fit against the data*
1) Run code: ```mkdir ./MCMC_model/TestFit/Plots``` & ```mkdir ./MCMC_model/TestFit/Data```
2) Open ```051122_MunicFitTest.R```
3) Set the number of iterations to match the number of iterations in your model run under variable ```iters`` 


### Simulations to Test Model Function <br />



### Model Projection Simulations <br />

## Relative Risk Analysis



## Fitness Analysis
