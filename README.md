This repository contains the scripts and functions developed for the estimation and analysis of MS-VECM (Markov Switching Vector Error Correction Models). This work is part of a Stochastic Processes thesis titled:

    "Macroeconomics factors influencing the Eurostoxx 50 under different regimes".

Core Scripts

    Main.R: It performs the time series analysis on Eurostoxx 50 and macroeconomic factors, leveraging the custom msvecm package.
    Alternative_MSwM_tsdyn.R: A workaround solution designed to handle multivariate time series, bypassing the current architectural limitations of the MSwM package.

Functions & Data Management

    R_functions/:  A dedicated folder containing all custom-built functions required to run the MS-VECM models and perform diagnostic tests.
    import_ecbdata.R: A utility function to import external data.
        Note: For immediate use, the processed dataset is already stored in obj.rData.
    data_simulation.R: A script specifically designed to simulate time series data for model validation and testing.
