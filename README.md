# Analysis for Geographic Mobility and Fitness of *Streptococcus pneumoniae*

## Human Mobility Model Analysis
This analysis integrates human mobility data and genomic data to mechanistically understand pneumococcal migration in South Africa.
### Process raw input data.  

Raw input data is in this folder:
```./scripts/manage_rawData/```
1) **Landscan Data** (https://landscan.ornl.gov/) <br />
*This script reads in the landscan population data and associated shapefiles from GADM (https://gadm.org/data.html) to create a datatable with all population levels ad the municipality (N=234) level population sizes.*  
Script: ```031122_LandScanMunic.R``` <br />
Output: ```./data/landscan2017/LandScan_PopulationN.RData```  <br />
2) **Facebook Data** <br />(https://dataforgood.facebook.com/) <br />
*This script reads in the raw facebook mobility data from the disaster movement range maps (data downloaded from FigShare URL)*<br />
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
3)**Metadata to Generate Likelihood**<br />
*This creates a list of length 9 (number of GPSCs) containing 8 columns with information for each pair including ids, collection location (province numbered 1-9), time in days between pairs, collection year for each pair.*<br />
Script: ```031122_GPSC_metadata.R```<br />
Output: Metadata File: ```./modelinput_data/dat.tmp.allser.RData```<br />
Population size (province): ```./modelinput_data/pop_2019.RData: ```<br />
4) **tMRCA between pairs**<br />
*This creates a list of length 9 (number of GPSCs) with pairwise time to MRCA for each pair (half the time between pairs (years)).*<br />
Script: ```031122_makeGDistMatrix.R```<br />
Output: ```./modelinput_data/tMRCAs.RData```<br />

## Relative Risk Analysis



## Fitness Analysis
