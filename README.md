# trEnD: Model for Transition of End-users and the Distribution Power Network
Welcome to trEnD, a c++ based model of the energy transition, focusing on the end-users and the distribution network. The source codes are currently undergoing revision and being tested by real data of Norway in my phD project. Feel free to comment on any issues.

# User Guide
## Preparing the Input Files

## Compiling the Source Codes
Currently 3 external libraries are used in the source codes of the model. The library [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [Boost](https://www.boost.org/) are required to compile the source codes. The library [ALGLIB](https://www.alglib.net/) is already in the source codes so you need not to download it in addition.

## Output Files

# Documentation of Source Codes
Below is a summary of what the source codes do, categorized by the subfolders in [src](https://github.com/TonyYenTWN/trEnD/tree/main/src) folder. For full documentation, see [the Github page of this project](https://tonyyentwn.github.io/trEnD/).

- [agent_operational](https://github.com/TonyYenTWN/trEnD/tree/main/src/agent_operational): source codes for the operation strategies of relevant agents in the model.
- [alglib](https://github.com/TonyYenTWN/trEnD/tree/main/src/alglib): source codes from the external library [ALGLIB](https://www.alglib.net/).
- [basic](https://github.com/TonyYenTWN/trEnD/tree/main/src/basic): source codes for generic functions and headers.
- [power_market](https://github.com/TonyYenTWN/trEnD/tree/main/src/power_market): source codes for the market clearing at different levels.
- [power_network](https://github.com/TonyYenTWN/trEnD/tree/main/src/power_network): source codes for the processing of power network data.
- [spatial_field](https://github.com/TonyYenTWN/trEnD/tree/main/src/spatial_field): source codes for the processing of spatial fields.
- [main.cpp](https://github.com/TonyYenTWN/trEnD/blob/main/src/main.cpp): the cpp file where the main function is located.
