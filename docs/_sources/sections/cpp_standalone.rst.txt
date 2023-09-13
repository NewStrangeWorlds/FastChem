C++ stand-alone executable
==========================

The ``FastChem`` object class is designed to be easily coupled to other
models. In addition to the object class itself, we also provide a
stand-alone executable that can call the module with some simple input
scripts. This stand-alone version, however, only provides a very basic
functionality, such as reading in a specific temperature-pressure
profile that ``FastChem`` will be run for. The stand-alone version does,
for example, not provide more advanced capabilities, such as looping
over different metallicity values or C/O ratios. If you intend to use
``FastChem`` for such purposes, you need to adapt the code that calls
``FastChem``.

The source code that is responsible for calling the actual ``FastChem``
chemistry is located in the folder ``model_src/``. It is split across
three different files: ``model_src/model_main.cpp``, the actual main
program, ``model_src/read_config.h`` for reading in the config file, and
``model_src/save_output.h`` for managing the output. Thus, if you want
to add another parameter to the config file, you would need to edit
``model_src/read_config.h``, while changes to the format of the output
files can be made in ``model_src/save_output.h``. Changing the contents
of these files obviously require a re-compilation of the code.

Running the FastChem executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following a successful configuration and compilation via ``CMake``, the
``FastChem`` executable ``fastchem`` should be present in the root
directory. The executable is started via

.. code:: bash

   ./fastchem input/config.input

where the second argument is the location of the config file that is
explained in the next section. ``FastChem`` will read in a pre-defined
pressure-temperature structures, the location of which is also specified
in the config file. After a successful calculation, ``FastChem`` will
produce two output files with a detailed chemistry output and one with
diagnostic output. The location of these files is also contained in the
config file and its contents are discussed :ref:`here<sec:fc_cpp_output>`.


Config file
~~~~~~~~~~~

| The config file that ``FastChem`` will read in at the beginning
  contains all important parameters and file locations necessary to
  initialise the chemistry and to perform the calculations. The
  numerical methods that these parameters refer to are described in
  Paper II. An example of such an input file is located in the input
  folder: ``input/config.input``. While this config file allows to set
  the most important ``FastChem`` parameters, some more advanced ones
  are not contained in this file and can only be set by invoking special
  ``FastChem`` functions during runtime. This, in particular, refers to
  the use of the optional scaling factors as described in the appendix
  of Paper II. More information on activating these scaling factors can
  be found in the :ref:`description<sec:fastchem_methods>` of the object 
  class and the set of :ref:`external parameters<sect:fc_param>`.
| The config file used for the ``C++`` stand-alone executable has the
  following structure:

.. code:: 

   #Atmospheric profile input file
   input/example_p_t_structures/Late_M-dwarf.dat

   #Chemistry calculation type
   ce

   #Chemistry output file
   output/chemistry.dat output/condensates.dat

   #Monitor output file
   output/monitor.dat

   #FastChem console verbose level (1 - 4)
   1

   #Output mixing ratios (MR) or particle number densities (ND, default)
   MR

   #Element abundance file   
   input/element_abundances/asplund_2009.dat

   #Species data files    
   input/logK/logK.dat input/logK/logK_condensates.dat

   #Accuracy of chemistry iteration
   1.0e-4

   #Accuracy of element conservation
   1.0e-4

   #Max number of chemistry iterations 
   80000  

   #Max number internal solver iterations
   20000

It contains the required parameters in the following order:

-  Location of the file with the pressure-temperature structure the
   chemistry should be calculated for. The file should contain two
   columns, where the first one is the pressure in units of bar and the
   second one the temperature in K. Header lines will be ignored.

-  the type of chemistry calculation:

   -  ``g`` calculate only the gas phase without condensation

   -  ``ce`` calculate the gas phase with equilibrium condensation

   -  ``cr`` calculate the gas phase with rainout condensation
      approximation

-  Desired location and file name for the output of the gas phase
   species and condensates

-  Desired location and file name for the diagnostic output

-  Verbose level, where a level of ``1`` is almost silent and ``4``
   produces a lot of diagnostic output on the terminal. Increase this
   level if you encounter issues to identify the source of the problems.

-  The output format for the gas-phase species’ abundances. By default,
   ``FastChem`` will use number densities in units of cm\ :math:`^{-3}`.
   If you use the keyword ``MR``, mixing ratios will be used instead.
   Any keyword other than ``MR`` will result in the default option of
   using number densities.

-  Location of the file with the element abundances, see 
   :ref:`here<sec:fastchem_input_files>` for details

-  Location of the file with the thermochemical and stoichiometric data
   for all species, see :ref:`here<sec:fastchem_input_files>` for
   details. The first entry refers to the gas-phase species, while the
   second one is for the condensate data.

-  Relative accuracy of the chemistry iterations, used as convergence
   criterion (also for Newton’s method)

-  Relative accuracy for the checks of the element conservation

-  Maximum number of chemistry iterations

-  Maximum number of internal solver method iterations (Newton,
   Nelder-Mead & bisection methods)

In the input file, the number of iterations for the Newton, Nelder-Mead,
and bisection methods are assumed to be the same. This number can be
adjusted individually for each of these internal solvers by using the
``FastChem.setParameter`` method (see :ref:`here<sec:fastchem_methods>`). 
The convergence criterion for Newton’s
method is also set to the accuracy of the chemistry iteration by
default. This convergence criterion can also be changed by the
``FastChem.setParameter`` method.

