Project Structure
This repository contains the tools and scripts developed for the estimation and analysis of MS-VECM (Markov Switching Vector Error Correction Models).

Core Scripts

    Main.R: The core script of the thesis development. It performs the complete time series analysis leveraging the custom msvecm package.
    Alternative_MSwM_tsdyn.R: A workaround solution designed to handle multivariate time series, bypassing the current architectural limitations of the MSwM package.

Functions & Data Management

    R_functions/: A dedicated directory containing all custom-built functions required to execute the MS-VECM workflow and perform statistical validation tests.
    import_ecbdata.R: A utility function to import external data.
        Note: For immediate use, the processed dataset is already stored in obj.rData.
    data_simulation.R: A script specifically designed to simulate time series data for model testing and scenario analysis.
