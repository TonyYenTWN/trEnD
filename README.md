# Project Introduction
Welcome to trEnD, a c++ based energy system model for the transition of end-users and the distribution network.
The codes are for my PhD project and it is undergoing revision. Feel free to comment on any issues.

The library [ALGLIB](https://www.alglib.net/), [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), and [Boost](https://www.boost.org/) are used to compile the codes.

# Categorization Rationale of Files
- Input files are csv files imported into the project.
- Processed files are csv files written by the source code and will be used later in the workflow of the project.
- Output files are csv files written by the source code as the result of the project.

# Comprehensive Documentation of Folders and Files
## [agent_operational](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/agent_operational)
Operational strategies of agents relevant to the distribution power network.

## [area_demand](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/area_demand)
Inference of the mean electricity demand per area (aka density field).

### [Input files](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/area_demand/input)

### Program files
- [area_demand.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/area_demand/area_demand.cpp): main source file for the inference of mean electricity demand density field using Bayesian maximum entropy method.
- [Geostat.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/area_demand/Geostat.cpp): source file for geostatisitc related functions, including the calculation of geodestic on an ellipsoid.
- [Geostat.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/area_demand/Geostat.h): header file for [Geostat.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/area_demand/Geostat.cpp).

### Output files

## [basic](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/basic)
Basic functions and header files for generic usage, including

- [Basic_Definitions.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/Basic_Definitions.h): header file for basic definitions.
- [Eigen_Sparse.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/Eigen_Sparse.h): header file for using sparse matrices from the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library.
- [LP_gpa.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/LP_gpa.cpp): source file of a linear programming solver using gradient projection algorithm.
- [LP_gpa.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/LP_gpa.h): header file for [LP_gpa.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/LP_gpa.cpp).
- [rw_csv.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/rw_csv.cpp): source file for reading and writing csv files.
- [rw_csv.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/rw_csv.h): header file for [rw_csv.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/rw_csv.cpp).

## [power_market](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/power_market)
Operation of the power market in Norway and between its neighbors.

### [Input files](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/power_market/input)

### Program files

### Output files

## [power_network](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/power_network/)
Physical model of the transmission and distribution power network in Norway.

### [Input files](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/tree/main/power_network/input)

### Program files

### Output files
