# Description of how to run the Fitness model 
## Vaccine type and non-vaccine type serotypes

## Penicillin resistance in vaccine type and non-vaccine type serotypes
### Process the raw data
Script: ```./Data/alltogether/111122_fitness.formatting.VaccineStatus_AMR_AMRVaxStat.R```<br />
Output: This will output the files with the numbers of vaccine type and non-vaccine type, AMR, and whether the VT or NVT had penicillin resistance. <br />

### Format the data input file for the model
Script: ```./Data/DataCreators/data_creator_fitness_model_AMRVT_together.R```<br />
Output: ```./Data/alltogether/Data_model_061622_ref_R.VT_fitness.rds```<br />

## Set up to run the model
1) mkdir ```RunModel```
2) 
