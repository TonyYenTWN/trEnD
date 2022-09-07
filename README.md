![alt text](https://github.com/TonyYenTWN/trEnD/blob/main/pic/main_logo.jpg?raw=true)
# trEnD: Model for Transition of End-users and the Distribution Power Network
Welcome to trEnD, a c++ based model of the energy transition, focusing on the end-users and the distribution network. The source codes are currently undergoing revision and being tested by real data of Norway in my phD project. Feel free to comment on any issues.

# User Guide
## Preparing the Input .csv Files
An example of the input csv files can be found [here](https://github.com/TonyYenTWN/trEnD/tree/main/csv/input). In case the model is applied for other datasets, the format of the csv files should be kept the same. Below is a list of the input .csv files and their format.

### .csv Files for Power Market
The following files can be found in the folder [power_market](https://github.com/TonyYenTWN/trEnD/tree/main/csv/input/power_market). If not specifically mentioned, the time series extracted from [ENTSO-E Transparency Platform](https://transparency.entsoe.eu/).

- [cbt_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/cbt_forecast_2021.csv) stores the time series of scheduled cross border transmission flows of the modeled bidding zones and their directly connected neighbors. Since there can be two directions of transmission flow between 2 connected bidding zones, each connection contributes to 2 columns of time series.
- [control_reserve_activated_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/control_reserve_activated_2021.csv) stores the time series of activated control reserve in the modeled bidding zones. Since the activated control reserve can be both positive and negative, each bidding zone contributes to 2 columns of time series.
- [da_price_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/da_price_2021.csv) stores the time series of day-ahead prices in the modeled bidding zones and their directly connected neighbors. The prices in Great Britian was extracted from [Elexon Portal](https://www.elexonportal.co.uk/).
- [generation_solar_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_solar_forecast_2021.csv) stores the time series of scheduled solar power generation in the modeled bidding zones and their directly connected neighbors.
- [generation_total_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_total_forecast_2021.csv) stores the time series of  scheduled total generation in the modeled bidding zones and their directly connected neighbors.
- [generation_wind_offshore_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_wind_offshore_forecast_2021.csv) stores the time series of  scheduled offshore wind power generation in the modeled bidding zones and their directly connected neighbors.
- [generation_wind_onshore_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_wind_onshore_forecast_2021.csv) stores the time series of  scheduled onshore wind power generation in the modeled bidding zones and their directly connected neighbors. The data of Great Britian were extracted from [Elexon Portal](https://www.elexonportal.co.uk/).
- [merit_order_curve_q_assimilated_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/merit_order_curve_q_assimilated_2021.csv) stores the merit order curves of the modeled bidding zones and their directly connected neighbors.
- [solar_radiation_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/solar_radiation_2021.csv) stores the time series of solar radiation data from meteorology stations in Norway.


### .csv Files for Power Network
The files can be found in the folder [power_network](https://github.com/TonyYenTWN/trEnD/tree/main/csv/input/power_network).

## Compiling the Source Codes
### Prerequisite for Compilation
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) libray installed
- [Boost](https://www.boost.org/) library installed

The library [ALGLIB](https://www.alglib.net/) is also used in the model, but it is already in the source codes so you need not to download it in addition.

Once compilation is complete, the input csv files, source files, object files and binary file will take about 120 MB of space.

### Compiler Choice and Settings
I use [TDM-GCC 9.2.0](https://jmeubank.github.io/tdm-gcc/articles/2020-03/9.2.0-release) to compile the codes on windows, but other compilers should also work.

I use [Code::Blocks](https://www.codeblocks.org/docs/main_codeblocks_en.html) for project management, and [a codeblock project file](https://github.com/TonyYenTWN/trEnD/blob/main/trEnD.cbp) is provided in the repository. Again it should be possible to use other IDEs for project management and code editing.

Remember to add the project's top-level directory (the directory which stores this repository) to compiler search directories when compiling. This can be done quite easily in codeblocks.

## Output Files

# Documentation of Source Codes
Below is a summary of what the source codes do, categorized by the subfolders in [src](https://github.com/TonyYenTWN/trEnD/tree/main/src) folder. The full documentation can be found at [the github-pages of this repository](https://tonyyentwn.github.io/trEnD).

- [agent](https://github.com/TonyYenTWN/trEnD/tree/main/src/agent): source codes for the operation and investment strategies of relevant agents in the model.
- [alglib](https://github.com/TonyYenTWN/trEnD/tree/main/src/alglib): source codes from the external library [ALGLIB](https://www.alglib.net/).
- [basic](https://github.com/TonyYenTWN/trEnD/tree/main/src/basic): source codes for generic functions and headers.
- [power_market](https://github.com/TonyYenTWN/trEnD/tree/main/src/power_market): source codes for the market clearing at different levels.
- [power_network](https://github.com/TonyYenTWN/trEnD/tree/main/src/power_network): source codes for the processing of power network data.
- [spatial_field](https://github.com/TonyYenTWN/trEnD/tree/main/src/spatial_field): source codes for the processing of spatial fields.
- [main.cpp](https://github.com/TonyYenTWN/trEnD/blob/main/src/main.cpp): the cpp file where the main function is located.
