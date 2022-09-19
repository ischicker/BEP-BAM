# Multivariate postprocessing techniques

The file `PP_example_thesis.R` file contains the code to construct the example that is used in the paper to illustrate the postprocessing techniques.

The files that are used in the simulations can be found in the `multiv_pp-master` folder and is adapted from the code of Lerch et al., 2020, https://doi.org/10.5194/npg-2019-62, where `Settings.R` contains several simulation parameter values.

## Simulations

The `Simulation code` folder contains all the code that is used in simulations. Here the files `ECC_T2M_Emos_subfunctions.R` and `emos_T2M_mean_singleForecast_subfunctions-orig.R` contain several helper functions for univariate postprocessing and computing scores. The actual postprocessing can be done from the `run_setting_archimedean.R` file and saves the results to the `Data/Rdata` folder. It sources all files from the `sourceArchimedean` folder for computing scores, univariate and multivariate postprocessing code. Further processing of the data is done in the `processing code for simulation output` folder, where test statistics are computed and saved in the `DM_TestStatistic_computation_Archimedean.R` file from the `Data/TestStatistic` folder.

## LAEF data

The data from the LAEF system can be imported with the `getData.R` file, where the data is cleaned and transferred to the global environment. The file `process_data.R` provides several initial visualizations of the data. The main file for the multivariate postprocessing is `mvpp_data.R`, which loads in the data after which it is postprocessed and sources from the `source` folder for postprocessing routines. Lastly the results are saved in the `Data/Rdata_LAEF` folder. Further processing of the data is done in the `processing code for simulation output` folder, where test statistics are computed and saved in the `DM_TestStatistic_computation_data.R` file from the `Data/TestStatistic` folder.

## Visualizations

The postprocessed data are visualized in various ways. All files relating to this are present in the `reproducing results and figures` folder, with the following usages:

- `ArchimedeanPlots.R`: depiction of typical two-dimensional copula output to showcase tail behavior.
- `ParamFitting.R`: visualization to show how well the copula parameter can be fitted from the data.
- `TauScatterPlot.R`: visualization of parameter output of copula fitting versus the copula parameter from observations.
- `TimeBoxPlots.R`: creates a figure to illustrate computational differences in postprocessing methods for simulations
- `TimeBoxPlots_data.R`: creates a figure to illustrate computational differences in postprocessing methods for the LAEF data
- `multivariate_Archimedean.R`: creates the box-plots of the Diebold-Mariano test statistics for various settings of the simulation and scoring rules.
- `multivariate_data.R`: creates the box-plots of the Diebold-Mariano test statistics for the LAEF data and various scoring rules.
- `mv_hist_data.R`: creates multivariate rank histograms for the data.
