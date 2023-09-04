# FastChem Cond (FastChem 3.0) #
#### Authors: Daniel Kitzmann, Joachim Stock ####

# Overview #

FastChem is an equilibrium chemistry code that calculates the chemical composition of the gas and condensed phase for given temperatures and pressures. The calculation of the gas phase is based on a semi-analytic approach, described in detail in Stock et al. (2018) and Stock et al. (2022). The new version 3.0 version of FastChem, called FastChem Cond, adds condensation to the code. It can now compute the chemical composition using equilibrium condensation or the rainout approximation that is commonly used in the field of exoplanets or brown dwarfs. FastChem Cond is described in detail in Kitzmann, Stock & Patzer (2023).

The code is written in object-oriented C++, including template programming that allows the model to run with either double or long double precision. The exact computational precision of long double depends on your compiler and operating system. Long double precision usually allows the model to properly converge for very low temperatures. FastChem has been tested for temperatures as low as 100 K. For many cases, we were also able to obtain converged results for temperatures well below 100 K. Overall, the model has been successfully tested for temperatures from 100 K to 6000 K and pressures from 1e-13 bar to 1000 bar for solar element abundances.

Besides the actual FastChem model, we provide a C++ stand-alone version in the *model_src* folder that allows to calculate the equilibrium chemistry for a given temperature-pressure structure. This stand-alone version can be adapted further to perform also more sophisticated calculations. In addition to the C++ stand-alone version, we also provide the Python interface PyFastChem.


# PyFastChem #

FastChem includes the Python interface, PyFastChem, that allows to run the C++ code as a normal Python module. PyFastChem is also available as a PyPy package and can easily be installed via pip. We provide several examples that show how to call FastChem from within a Python script, including the possibility to, for example, iterate over different metallicity values or C/O ratios or how to compute chemical compositions with condensation. The Python examples can be found in the *python* directory.  


# User guide #

FastChem comes with user guide in pdf format that can be found in the *manual* folder. It describes the installation and usage of FastChem and covers both the C++ stand-alone version as well as the Python interface. The manual also provides detailed information on the internal interface functions that the FastChem object class and its Python interface provide.


# Licence #

This project is Copyright (c) Daniel Kitzmann and Joachim Stock.

FastChem is released under the GNU Public Licence (GPL) 3.0. That means, it can be freely copied, edited, and re-distributed. If the code is re-distributed it has to be released under at least a GPL 3.0 licence as well. The full licence of FastChem can be found in the repository *licence.md* file or under https://www.gnu.org/licenses/gpl-3.0.html.

The user guide is released under the Creative Commons Licence (CC BY SA). Licensees may copy and distribute the work and make derivative works based on it only if they give the authors (Daniel Kitzmann & Joachim Stock) the credits by providing a reference to the original guide and this GitHub repository. Licensees may also distribute derivative works only under a license identical to ("not more restrictive than") the license that governs the original work.

