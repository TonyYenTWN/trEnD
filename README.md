![alt text](https://github.com/TonyYenTWN/trEnD/blob/main/logos/main_logo.jpg?raw=true)
# trEnD: Model for Transition of End-users and the Distribution Power Network
Welcome to trEnD, a c++ based model of the energy transition, focusing on the end-users and the distribution network. The source codes are currently undergoing revision and being tested by real data of Norway in my phD project. Feel free to comment on any issues.

# User Guide
## Preparing the Input Files

## Compiling the Source Codes
### Prerequisite for Compilation
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) libray installed
- [Boost](https://www.boost.org/) library installed

The library [ALGLIB](https://www.alglib.net/) is also used in the model, but it is already in the source codes so you need not to download it in addition.

### Compiler Choice and Settings
I use compiler [TDM-GCC 9.2.0](https://jmeubank.github.io/tdm-gcc/articles/2020-03/9.2.0-release) to compile the codes on windows, but other compilers should also work.

I use [Code::Blocks](https://www.codeblocks.org/docs/main_codeblocks_en.html) for project management, and [a codeblock project file](https://github.com/TonyYenTWN/trEnD/blob/main/trEnD.cbp) is provided in the repository. Again it should be possible to use other IDEs for project management and code editing.

Remember to add the project's top-level directory (the directory which stores this repository) to compiler search directories when compiling. This can be done quite easily in codeblocks.

## Output Files

# Documentation of Source Codes
Below is a summary of what the source codes do, categorized by the subfolders in [src](https://github.com/TonyYenTWN/trEnD/tree/main/src) folder. The full documentation can be found at [the github-pages of this repository](https://tonyyentwn.github.io/trEnD).

- [agent_operational](https://github.com/TonyYenTWN/trEnD/tree/main/src/agent_operational): source codes for the operation strategies of relevant agents in the model.
- [alglib](https://github.com/TonyYenTWN/trEnD/tree/main/src/alglib): source codes from the external library [ALGLIB](https://www.alglib.net/).
- [basic](https://github.com/TonyYenTWN/trEnD/tree/main/src/basic): source codes for generic functions and headers.
- [power_market](https://github.com/TonyYenTWN/trEnD/tree/main/src/power_market): source codes for the market clearing at different levels.
- [power_network](https://github.com/TonyYenTWN/trEnD/tree/main/src/power_network): source codes for the processing of power network data.
- [spatial_field](https://github.com/TonyYenTWN/trEnD/tree/main/src/spatial_field): source codes for the processing of spatial fields.
- [main.cpp](https://github.com/TonyYenTWN/trEnD/blob/main/src/main.cpp): the cpp file where the main function is located.
