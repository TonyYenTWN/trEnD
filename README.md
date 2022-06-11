# distribution_grid_transition_Norway
The codes are for my PhD prgram and it is undergoing revision. Feel free to comment on any issues.

The library [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](https://www.boost.org/) are required to compile the codes.

## agent_operational
Operational strategies of agents relevant to the distribution power network.

## area_demand
Inference of the mean electricity demand per area.

### Input files

### Program files

### Output files

## basic
Basic functions and header files for generic usage, including

- [Eigen_Sparse.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/Eigen_Sparse.h): header file for using sparse matrices from the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library.
- [LP_gpa.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/LP_gpa.cpp): source file of a linear programming solver using gradient projection algorithm.
- [LP_gpa.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/LP_gpa.h): header file for [LP_gpa.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/LP_gpa.cpp).
- [rw_csv.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/rw_csv.cpp): source file for reading and writing csv files.
- [rw_csv.h](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/rw_csv.h): header file for [rw_csv.cpp](https://github.com/TonyYenTWN/distribution_grid_transition_Norway/blob/main/basic/rw_csv.cpp).

## power_market
Operation of the power market in Norway and between its neighbors.

### Input files

### Program files

### Output files

## power_network
Physical model of the transmission and distribution power network in Norway.

### Input files

### Program files

### Output files
