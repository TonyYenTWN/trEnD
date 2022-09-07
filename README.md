![alt text](https://github.com/TonyYenTWN/trEnD/blob/main/pic/main_logo.jpg?raw=true)
# trEnD: Model for Transition of End-users and the Distribution Power Network
Welcome to trEnD, a c++ based model of the energy transition, focusing on the end-users and the distribution network. The source codes are currently undergoing revision and being tested by real data of Norway in my phD project. Feel free to comment on any issues.

# User Guide
## Preparing the Input .csv Files
An example of the input csv files can be found [here](https://github.com/TonyYenTWN/trEnD/tree/main/csv/input). In case the model is applied for other time periods of the same region, the format of the csv files should be kept the same. For regions where conventional thermal power plants still play a role in the power system, additional input data might be needed, and corresponding parts of the source codes might therefore have to be modified in the future. 

Below is a brief summary of the input .csv files currently used in the model.

### .csv Files for Power Market
The following files can be found in the folder [power_market](https://github.com/TonyYenTWN/trEnD/tree/main/csv/input/power_market). If not specifically mentioned, the time series were extracted from [ENTSO-E Transparency Platform](https://transparency.entsoe.eu/).

- [cbt_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/cbt_forecast_2021.csv) stores the time series of scheduled cross border transmission flows of the modeled bidding zones and their directly connected neighbors. Since there can be two directions of transmission flow between 2 connected bidding zones, each connection contributes to 2 columns of time series.
- [control_reserve_activated_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/control_reserve_activated_2021.csv) stores the time series of activated control reserve in the modeled bidding zones. Since the activated control reserve can be both positive and negative, each bidding zone contributes to 2 columns of time series.
- [da_price_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/da_price_2021.csv) stores the time series of day-ahead prices in the modeled bidding zones and their directly connected neighbors. The prices in Great Britian was extracted from [Elexon Portal](https://www.elexonportal.co.uk/).
- [generation_solar_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_solar_forecast_2021.csv) stores the time series of scheduled solar power generation in the modeled bidding zones and their directly connected neighbors.
- [generation_total_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_total_forecast_2021.csv) stores the time series of  scheduled total generation in the modeled bidding zones and their directly connected neighbors.
- [generation_wind_offshore_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_wind_offshore_forecast_2021.csv) stores the time series of  scheduled offshore wind power generation in the modeled bidding zones and their directly connected neighbors.
- [generation_wind_onshore_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/generation_wind_onshore_forecast_2021.csv) stores the time series of  scheduled onshore wind power generation in the modeled bidding zones and their directly connected neighbors. The data of Great Britian were extracted from [Elexon Portal](https://www.elexonportal.co.uk/).
- [merit_order_curve_q_assimilated_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/merit_order_curve_q_assimilated_2021.csv) stores the merit order curves of the modeled bidding zones and their directly connected neighbors. This file might not be needed in the future once merit order curves can be modeled endogeneously in the model.
- [solar_radiation_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/solar_radiation_2021.csv) stores the time series of solar radiation data from meteorology stations in the modeled bidding zones.

### .csv Files for Power Network
The files can be found in the folder [power_network](https://github.com/TonyYenTWN/trEnD/tree/main/csv/input/power_network). If not specifically mentioned, the time series were extracted from the [NVE dataset](https://nedlasting.nve.no/gis/).

- [DSO_Bidding_Zone.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/DSO_Bidding_Zone.csv) stores the bidding zones where original DSO areas from the NVE dataset are located.
- [cbt_constraint.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/cbt_constraint.csv) stores the cross border transmission flow constraint between the modeled bidding zones and their directly connected neighbors. "-1" means that there are no connections between the 2 bidding zones. This dataset was inferred from [cbt_forecast_2021.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_market/cbt_forecast_2021.csv).
- [cbt_entry_nodes.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/cbt_entry_nodes.csv) stores the transmission nodes where power flows from / to neighboring bidding zones injects / leaves. These transmission nodes were inferred from the transmission network graph from [PyPSA-Eur](https://pypsa-eur.readthedocs.io/).
- [hydro_plants.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/hydro_plants.csv) stores the information of hydroelectric power plants in the modeled bidding zones.
- [point_info.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/point_info.csv) stores the information of spatial points considered in the model.
- [point_matrix.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/point_matrix.csv) stores the geometrical relations of the spatial points in [point_info.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/point_info.csv) with a matrix. It makes finding the neighbors of a spatial point faster.
- [solar_radiation_stations.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/solar_radiation_stations.csv) stores the information of the meteorology stations in the modeled bidding zones.
- [transmission_edges.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/transmission_edges.csv) stores the information of the edges of the transmission network in the modeled bidding zones.
- [transmission_nodes.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/transmission_nodes.csv) stores the information of the nodes of the transmission network in the modeled bidding zones.
- [wind_plants.csv](https://github.com/TonyYenTWN/trEnD/blob/main/csv/input/power_network/wind_plants.csv) stores the information of wind power plants in the modeled bidding zones.

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
