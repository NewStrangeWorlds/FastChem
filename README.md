# FastChem #
#### Authors: Daniel Kitzmann, Joachim Stock ####

This readme is currently incomplete and will be expanded in the near future. We will also provide a guide on using the internal methods provided by FastChem and  more details on the actual implementation of the methods discussed in Stock et al. (2018).

# Overview #

FastChem is an equilibrium chemistry code that calculates the chemical composition of the gas phase for given temperatures and pressures. It is based on a semi-analytic approach, described in detail in Stock et al. (2018). The code is optimised for extremely fast and accurate calculations. 

The code is written in object-oriented C++, including template programming that allows the model to run with either *double* or *long double* precision. The exact computational precision of *long double* depends on your compiler and operating system. *Long double* precision usually allows the model to properly converge for very low temperatures. FastChem has been tested for temperatures as low as 100 K. For many cases, we were also able to obtain converged results for temperatures well below 100 K. As shown by Stock et al. (2018), the model has been successfully tested for temperatures from 100 K to 2500 K and pressures from 1e-13 bar to 1000 bar.

The FastChem code itself has been designed to be coupled to other models. Because FastChem is written as a C++ object class, several instances of FastChem can be created and run simultaneously within your model code. For example, it is possible to have a *double* and *long double* version of FastChem working in parallel. The object class contains several methods that allow you to query the chemistry for certain molecules and their corresponding chemical symbols or names.

Besides the FastChem code itself, this release also provides several demo wrappers that show how FastChem is called from inside other programs. 
These demos can be used as stand-alone versions of FastChem, but are quite limited in terms of input and output options. You might want to modify them if you intend to use FastChem as a stand-alone code.



## Demo stand-alone versions ##

This release provides five different demos that can also act as stand-alone versions to calculate the gas phase chemical equilibrium composition for a given temperature-pressure structure.

The aim of these demos is, to give an overview of the different ways to call FastChem within another code. It also demonstrates how to query FastChem for indices, the names, and symbols of chemical species. Each demo shows a certain way to run FastChem. The following demos are provided in this release:

* **demo1**: standard version with long double and OpenMP support
* **demo2**: long double version, calls FastChem with a temperature/pressure array
* **demo3**: long double version, calls FastChem with single p-T points
* **demo4**: double version, calls FastChem with single p-T points
* **demo5**: same as above but with no diagnostic output


# Compiling and running FastChem #

Before you can run FastChem on your computer, the code needs to be compiled and linked to create the executable file. Additionally, you may need to edit the config files. 
If you intend to move FastChem from one computer to another, you should always recompile the code, even if both computers use the same operating system. The standard options in FastChem's makefile result in optimisations that are tied to the specific CPU type the code was compiled on. Running it on a different CPU type might result in degraded computational performance.

### Compiling the demo codes ###

FastChem includes a makefile that should make the compilation process easy. The makefile is currently configured to use the g++ compiler from the GNU Compiler Collection and requires no additional libraries. In case you want to use a different compiler, you have to edit the *make\_global.options* file with your compiler's executable file and the corresponding compiler options.
Note that in case you intend to use the demo1 version, that includes OpenMP support, the OpenMP extensions of the compiler need to be installed on your computer. If those are not present, you have to us another of the provided demo codes instead (e.g. demo2). To compile and link any of the demos, just type

> make demo\#   

in your console, where \# is the number of the demo version you want to compile. Just using 

> make  
 
without any further parameter will automatically compile the demo1 version.
After successful compilation and linking, the FastChem executable *fastchem* can be found in the root directory. You can change the name of the executable in the *make\_global.options* file, if necessary. Compiled object files are placed in the *obj* folder. 

Compilation has, so far, only been tested with the gcc compiler version 5.x. Since FastChem is written in standard C++, other compilers should work as well. Note, however, that the compiler needs to provide at least the C++11 standard to compile FastChem.

### Removing compiled files ###

Calling the makefile with

> make clean
 
will delete all compiled object files and the executable. Use this command, if you want to recompile all files from scratch (e.g. when moving FastChem to a different computer).


## Configuration files ##

The demo stand-alone versions of FastChem require two config files: one for the input and output of the chemistry calculations and one configuration file for the FastChem code itself. If you intend to couple FastChem to your own code, then you only need the latter one.

The input file for the demos has the following structure:

> \#FastChem parameter file  
> input/parameters.dat
>
> \#Atmospheric profile input file  
> input/AGB_stellar_wind.dat
>
> \#Abundances output file  
> output/chem_output.dat
>
> \#Monitor output file  
> output/monitor_output.dat
>
> \#FastChem console verbose level (1 - 4); 1 = almost silent, 4 = detailed console output  
>1
>
>\#Output mixing ratios (MR) or particle number densities (ND, default)?  
>MR

The first parameter is the location of the FastChem config file (see below), the second parameter the temperature-pressure input file, the third the location of the chemistry output file, the fourth is the location of the monitor file for the FastChem calculation, the sixth the FastChem verbose level (i.e. the amount of console output), and the seventh the output format (either mixing ratios MR or number densities ND). You can find an example of this file in the *input* folder.

The second, required config file is the one for FastChem itself. It contains the location of the element abundances and configuration parameters. If you intend to couple FastChem to your own code, the location of this file has to be supplied within the constructor of the FastChem class. Note that most of the parameters within this file can also be later adjusted using the corresponding methods provided by the FastChem class. The config file has the following structure:

>\#element abundance file  
> input/element\_abundances\_solar.dat
>  
>\#elements data file         
>chemistry/chem_input/chemical\_elements.dat
>
>\#species data file    
>input/logK.dat
>
>\#accuracy of chemistry iteration  
>1.0e-4
>
>\#accuracy of pressure iteration  
>1.0e-4
>
>\#max error of Newton's method  
>1.0e-4
>
>\#max number of chemistry iterations   
>300  
>
>\#max number of pressure iterations                     
>100
>
>\#max number of Nelder-Mead iterations  
>20000

The first parameter is the location of the file with the element abundances. You can find an example of this file with the abundances from Asplund et al. (2009) in the *input* directory of FastChem. The second parameter is the location of the input file with information about the chemical elements (such as their atomic masses). The standard version of this file can be found within the *fastchem\_src/chem\_input* directory. 
The third parameter is the location of the input file with information about the stoichiometry and mass action constants. The mass action constants are fitted according to Stock et al. (2018) based on thermochemical data from Barin (1995), Burcat and Ruscic (2005), Chase (1998), Goos et al. (2016) and Tsuji (1973). Molecules can be removed from or added to the list by deleting the corresponding line or adding a line in the same format with new user provided thermochemical data. The user should feel encouraged to optmise this file according to his or her needs. There are two different files inlcuded in the input folder, *logK.dat* for the full set of chemical species listed in Stock et al. (2018) and *logK_wo_ions* for the same set but without ions. The next three parameters set the accuracy of the chemistry and pressure iterations, as well as the convergence criterion for Newton's method. The last three parameters are the maximum number of iterations for the chemistry, the pressure iteration, and the Nelder-Mead method. 

## Temperature-pressure input file ##

The stand-alone versions of FastChem require an input file with temperature and pressure data. As explained above, the location of this file has to be listed in the first config file. The structure of the temperature-pressure input file is simple: the first line contains a header with information on the input data. When reading the input file, the stand-alone versions of FastChem will **always** skip the first line of the the file! Followed by this header, there should be two columns present in the file: the first column is the temperature in K, the second column the pressure in bar. Note that internally FastChem uses dyn cm-2 as units of pressure. The stand-alone versions will convert the pressure from the input file accordingly.
There are also some examples for temperature-pressure profiles of various atmospheres in the *input* folder that can be used as test cases. Their chemistry output can be found in the *output_benchmarks* folder (see below).

If you want to change the format of this input file, you have to edit the corresponding code in the *demos_src* folder.


## Starting a calculation ##

Starting the compiled FastChem demos is done via the console command

>./fastchem input/config.input

The executable must be called with the location of the main config file as a command line parameter. FastChem will terminate with an error message in case this file parameter is not present or the file could not be found. Depending on the verbose level set in the config file, FastChem will provide very detailed console outputs or be almost silent. With the lowest verbose level 1, FastChem will only report back critical errors, such as non-converging chemistry or pressure iterations.


## Output files ##

The demo versions of FastChem will create two output files, placed at the location indicated in the config file. The first output file contains the chemistry output. It will list either the mixing ratios or number densities of all species in the network, depending on what option is specified in the config file. Additionally, it will also repeat the temperatures and pressures that were used in the calculation, the total hydrogen density, and the mean molecular weight. The content of each column (including the unit) is marked in the first line of the output file. If you chose number densities as output option in the config file, cm-3 will be used as the corresponding unit for each species.

In a second file, the FastChem demos will provide a diagnostic output - except in the case of demo5 that runs FastChem without a such a diagnostic. The location of this monitor file is also supplied in the corresponding config file (see above). The monitor file will list the number of chemistry and pressure iterations for each temperature-pressure point. Additionally, it will also give information on the state of element conservation for each element in the network, as well as the charge conservation. If either of those is violated, *fail* will we printed for the corresponding element, otherwise *ok* if everything worked fine.

## Benchmark outputs ##

In the folder *output_benchmarks* we also provide the chemistry output for the sample of temperature-pressure profiles that can be found in the input folder. The files contain the full network provided with FastChem and were calculated by using the *demo1* version of FastChem for solar metallicities. The abundances of the chemical species are listed as mixing ratios. You can use those benchmarks to validate if FastChem is working correctly on your computer.


# Troubleshooting #

In case FastChem does not produce the desired output, yields strange results, or reports back errors, you can do the following to identify the root of the problem and solve it:

* Check the console output for error messages
* Increase the verbose level for more detailed output
* Increase the number of chemistry iterations (in case FastChem reports a non-converged chemistry)
* Increase the number of pressure iterations (in case the pressure iterations do not converge)
* Decrease/increase the accuracies for the iterations
* Check the monitor file for non-conserved elements - in this case you might need to increase the number of chemistry iterations


## Known issues ##

Including electrons and ions in the network can considerably decrease the computational performance, especially under conditions, where their number densities are expected to be extremely small (low temperature and high pressure). In this case, FastChem uses a Nelder-Mead downhill simplex method to calculate the electron density, which can require a lot of computational time. Thus, under conditions, where ions are not expected to form in large abundances, it is better to run FastChem without ions. Note that there are two files with equilibrium constants present in the *input* folder, one with the full network (*logK.dat*) and a second one without ions (*logK\_wo\_ions.dat*). For a calculation without ions you, therefore, just have to direct FastChem to use the second file in the parameter config file.

Under certain conditions, FastChem can require a rather large number of iterations to find the equilibrium solution. As explained in Stock et al, (2018), this often occurs when certain molecules compete for the same elements. This will be the case, for example, when the C/O ratio is close to unity. In this case, you have to increase the number of the chemistry iterations in FastChem's config file. Note, that this can sometimes require a number of iterations in excess of 5000.

Currently, FastChem is designed for hydrogen-dominated cases. We tested the code also for cases, where hydrogen is only a minor species. In principle, FastChem also works for such a case. However, the pressure iteration - that is tied to n\_\<H\> - will be **very** slow. We intend to improve this behaviour in the future. Additionally, FastChem will currently **not work at all** when hydrogen is missing in the chemistry network.


# References #

Asplund M., Grevesse N., Sauval A. J., Scott P., 2009, ARA&A, 47, 481  
Barin, I. (1995), Thermochemical Data of Pure Substances  
Burcat, A., Ruscic, B. (2005), Third Millenium Ideal Gas and Condensed Phase Thermochemical Database  
Chase, M., (1998), NIST-JANAF Thermochemical Tables  
Goos, E., Burcat, A., Ruscic, B. (2016), Extended Third Millenium Ideal Gas Thermochemical Database  
Stock, J., Kitzmann, D., Patzer, A.B.C., Sedlmayr, S. (2018), MNRAS submitted  
Tsuji, T. (1973), A&A, 23, 411
