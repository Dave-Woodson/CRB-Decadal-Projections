# CRB-Decadal-Projections

Read Me

For the decadal flow projections and MTOM simulations, see below steps for running code. 

1. For decadal forecasts, first run the .Rmd file "midTermForecast-v2.07-DPLE---MultiYears---Blind---RFonly.Rmd" in the /Code/Multi-year Projections/ folder. 
	a) Ensure that the associated library "midTermForecast_Library_v1.4.R" is loaded
	b) The required input data for this program is in the /Data/Input Data/ folder (including subdirectories)
	c) You should be able to run this program as an R notebook all at once and it will output the decadal forecasts and skill scores. 
	d) The output should match all the contents of /Data/Output Data and Analysis Results/Multi-year Projections/ 
	e) Aggregated forecasts and skill scores for all mean lengths are located in "RF_forecasts_skillscores_v2.07_blindKFold.Rds" and "RF_forecasts_v2.07_blindKFold.Rds" 

2. For MTOM simulations from the decadal forecasts, follow the below steps:
	a) In the /Code and Data/Code/MTOM simulations/ directory, first run the "WaterYear_Flow_Disaggregation - v1.8.R" file. 
		i)This dissaggregates the decadal forecasts to a monthly, sub-basin resolution.
		ii) You will need to change the directories of the input files (i.e., point to where the decadal forecasts from 1d and 1e are stored)
	b) Open the macros.xlsm file and go to the Developer Tab then open Macros. 
		i) Edit the 'ReplaceDate_RF' macro so that it point to where the disaggregated forecasts created in 2a are stored.
		ii) This macro changes the date format of the forecasts so it can be read by MTOM
	c) The forecasts are now ready for MTOM input. 
	d) After running MTOM via RiverSMART, run the 'mtom_plotting_loop - v1.9.R' code. Make sure the correct path is set for reading in the MTOM simulations and the 'mtom_postProcessing_library - v1.4.R' function library
		i) This code will calculate various performance metrics on the MTOM simulations and output related graphics.
