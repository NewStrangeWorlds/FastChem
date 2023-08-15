This folder contains important input data for FastChem.

The following data can be found here:
element_abundances/*         Folder with element abundances compilations

logK/*                       Folder with the equilibrium constants

example_p_t_structures/*     Folder with some temperature-pressure example structures

config.input                 Config file for the C++ stand-alone version

parameters.config            An example for a parameter file the C++ version can use
                             Note: this file is *not* explicitly used in the stand-alone version.
                             To use this file, the C++ object of FastChem needs to be created with the
                             corresponding constructor (see the commented-out example in 
                             model_src/model_main.cpp or the manual)

parameters_py.config         An example for a parameter file the Python version can use
                             Note: the C++ parameter file is almost identical to this one. Only 
                             the relative paths for the input files are different because the
                             Python files are called from within a separate directory. To use this
                             file, the Python object of FastChem needs to be created with the
                             corresponding constructor (see the example Python scripts or the manual)

