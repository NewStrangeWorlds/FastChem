.. _sec:pyfc_details:

Detailed pyFastChem module description
=================================================

| By using the library ``PyBind11``, ``FastChem`` can be called directly
  within Python. This requires the compilation of ``FastChem``\ ’s
  Python wrapper that is located in the file
  ``python/fastchem_python_wrapper.cpp``. ``pyFastChem`` currently only
  links to the ``long double`` ``C++`` version of ``FastChem``.
| As described :ref:`here<sec:install_config>`, when using the
  ``CMake`` approach, ``pyFastChem`` will be automatically compiled when
  ``cmake`` is configured with the corresponding option. If the
  configuration and compilation is successful, a module file should be
  present in the ``python/`` folder that contains the Python module
  which acts as a wrapper between Python and the ``C++`` version of
  ``FastChem``. The file should be named ``pyfastchem.cpython-xxxx``,
  where ``xxxx`` will be a combination of your Python version and
  operating system.
| If ``pyFastChem`` has been built using the ``setup.py`` script or
  installed via ``PyPI``, then the module will be located in your normal
  Python module library path. It, thus, can be accessed from everywhere
  on your system like any other standard Python package. The location
  and additional module information can be obtained via

.. code:: bash

   pip show pyfastchem

| Depending on your Python installation, ``pip`` might need to be
  replaced by ``pip3``.
| A description of the module is given :ref:`here<sec:pyfc_module>`. Besides the ``pyFastChem`` module, we also
  provide several example Python scripts that show how to call
  ``FastChem`` from within Python for several different scenarios. We
  discuss the examples :ref:`here<sec:pyfc_examples>`.

.. _sec:pyfc_module:


The ``pyFastChem`` module provides access to the ``FastChem`` object
class as well as additional constants used within ``FastChem``. They are
essentially identical to their ``C++`` counterparts discussed :ref:`here<sec:fc_class>`.

To include the ``pyFastChem`` module in your Python project, just import
it using

.. code:: python

   import pyfastchem

This provides access to the ``FastChem`` object class as well as the
input and output structures and additional pre-defined constants used by
``FastChem``.

.. _sec:pfc_constants:

pyFastChem constants
~~~~~~~~~~~~~~~~~~~~

| The ``pyFastChem`` module contains a number of pre-defined constants.
  This includes the constant ``pyfastchem.FASTCHEM_UNKNOWN_SPECIES`` of
  type ``int`` that is returned by some ``pyFastChem`` methods when a
  chemical species is not found.
| The chemistry calculation can also return several output flags of type
  ``int`` defined as constants in the ``pyFastChem`` module:

``pyfastchem.FASTCHEM_SUCCESS``
  Indicates that the calculation has been successful, i.e. that the chemistry iterations converged.

``pyfastchem.FASTCHEM_NO_CONVERGENCE``
  Indicates that the calculation was not successful, i.e. that the chemistry did not converge within the allowed maximum number of iterations steps given in the config file or set manually via the ``setParameter`` method together with the ``nbIterationsChem`` parameter (see :ref:`here<sec:pfc_methods>`). One way to solve such a problem is to increase the maximum number of iteration steps.

``pyfastchem.FASTCHEM_INITIALIZATION_FAILED``
  Indicates that something went wrong during reading one of the input files. To find the source of the problem, one can set the verbose level in the config file or manually via ``setVerboseLevel`` (see :ref:`here<sec:pfc_methods>`) to a higher value and look at the terminal output.

``pyfastchem.FASTCHEM_IS_BUSY``
  The chemistry calculations of ``FastChem`` can only be called once for each object class instance. Attempting to start a new calculation while another is still running will result in ``FastChem`` returning this flag.

``pyfastchem.FASTCHEM_WRONG_INPUT_VALUES``
  ``FastChem`` returns this flag if some input values are wrong. Currently, this refers to the temperature and pressure vectors in the input structure not having the same size (see :ref:`here<sec:fastchem_input_struct>` for details on the input structure).

``pyfastchem.FASTCHEM_PHASE_RULE_VIOLATION``
  ``FastChem`` returns this flag if condensation is used and the system violates the phase rule. This happens when the number of elements contained in condensates equals the total number of elements. In this case, the gas phase lacks a degree of freedom to yield the correct gas pressure. Such a system cannot be solved as there has always to be at least one incondensable element in the gas phase (see the section about the phase rule in Paper III).

In addition to these flags, the ``pyFastChem`` module also includes a
constant string array ``pyfastchem.FASTCHEM_MSG`` that contains string
expressions for each of these flags. Using this array with any of the
aforementioned flags ``pyfastchen.FASTCHEM_MSG[flag]`` returns a string
with a description of the corresponding flag’s meaning. For example,

.. code:: python

  pyfastchem.FASTCHEM_MSG[pyfastchem.FASTCHEM_NO_CONVERGENCE] 
  
will return the string ``"convergence failed"``.

pyFastChem constructor
~~~~~~~~~~~~~~~~~~~~~~~~~~

| Since ``FastChem`` is written as an object class, an instance of that
  class (i.e. an object) needs to be created before ``FastChem`` can be
  used. This is done by calling the constructor of the ``FastChem``
  class that is contained within the ``pyFastChem`` module. There are
  three main ways to call the constructor and create an object.
  
.. code:: python

  pyfastchem.FastChem(str element_abundance_file, 
                      str gas_species_data_file, 
                      int verbose_level)

..

  This constructor requires three parameters: the locations of the element abundance and gas phase species data files, as well as the verbose level. All other options and parameters within ``FastChem`` will be set to their default values but can be later changed by using the appropriate methods described :ref:`here<sec:pfc_methods>`. The default maximum number of chemistry iterations is 3000, the number of Newton, bisection and Nelder-Mead method iterations is 3000, and the default accuracy of the of Newton method and the chemistry iterations is set to :math:`10^{-4}`. This constructor will not read in any condensate data. Trying to use an object created via this method for a calculation using condensation will result in an error message.

.. code:: python

  pyfastchem.FastChem(str element_abundance_file, 
                      str gas_species_data_file, 
                      str condensate_species_data_file, 
                      int verbose_level)

..

  This constructor requires four parameters: the locations of the element abundance and gas phase species data files, the condensate data file, as well as the verbose level. All other options and parameters within ``FastChem`` will be set to their default values but can be later changed by using the appropriate methods described :ref:`here<sec:pfc_methods>`. The default maximum number of chemistry iterations is 3000, the number of Newton, bisection and Nelder-Mead method iterations is 3000, and the default accuracy of the of Newton method and the chemistry iterations is set to :math:`10^{-4}`. Note that instead of a location for the condensate data, a ``str`` containing ``'none'`` can be used here as well. In that case, no condensate data will be read in and trying to use the object for a calculation using condensation will result in an error message.

.. code:: python

  pyfastchem.FastChem(str parameter_file, 
                      int intial_verbose_level)

..

  The constructor requires two different arguments: the location of the parameter file and the initial verbose level. The latter one will be replaced by the corresponding value read in from the parameter file. The structure of this parameter file is discussed :ref:`here<sec:fc_parameter_file>`. All of parameter values read in from the file can also be adjusted during runtime by using the methods listed :ref:`here<sec:pfc_methods>`.

.. _sec:pfc_input_ouput_struct:

pyFastChem input and output structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running a FastChem chemistry calculation requires input and output data
structures, resembling those of the ``C++`` :ref:`version<sec:fastchem_input_struct>`. 
In Python they are represented as classes rather than a ``C++``\ ``struct``.

Input structure
---------------
               

The original ``C++`` struct translated by ``PyBind11`` has the following
structure in Python:

.. code:: python

   class FastChemInput:
       temperature: list[float] = []
       pressure: list[float] = []

       equilibrium_condensation = False
       rainout_condensation = False

The input class contains the following variables:

``temperature``
  An array of ``float`` numbers that describe the temperature in K.

``pressure``
  An array of ``float`` numbers that describe the pressure in bar.

``equilibrium_condensation``
  A ``bool`` parameter that enables the calculation of equilibrium condensation. Its default value is ``False``.

``rainout_condensation``
  A ``bool`` parameter that enables the calculation of condensation via the rainout approximation. Its default value is ``False``. Note that when the flag ``rainout_condensation`` is set to ``True``, the value of the parameter ``equilibrium_condensation`` is ignored.

An input structure, in the example here called ``input_data``, can be
defined from the ``pyFastChem`` module in the following way:

.. code:: python

   input_data = pyfastchem.FastChemInput()

The two input arrays for temperature and pressure need to have the same
length. The ``PyBind11`` library allows normal Python lists or NumPy
arrays to be used here. For example, a NumPy array for the pressure
could be defined using NumPy’s ``logspace`` function:

.. code:: python

   input_data.pressure = np.logspace(-6, 1, num=1000)

Both of these input variables need to be an array-type variable, even if
only a single temperature-pressure point is going to be calculated.

.. _output-structure-1:

Output structure
----------------
                

The original ``C++`` output struct translated by ``PyBind11`` has the
following structure in Python:

.. code:: python

   class FastChemOutput:
       number_densities: list[list[float]]
       total_element_density: list[float]
       mean_molecular_weight: list[float]
       
       number_densities_cond: list[list[float]]
       element_cond_degree: list[list[float]]
       
       element_conserved: list[list[int]]
       nb_chemistry_iterations: list[int]
       nb_cond_iterations: list[int]
       nb_iterations: list[int]
       fastchem_flag: list[int]

It has the following variables:

``number_densities``
  The two-dimensional array contains the number densities in of all gas phase species (elements, molecules, ions) as ``float`` numbers. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure arrays of the input structure. The second dimension refers to the number of species and has a length of ``getGasSpeciesNumber()`` (see :ref:`here<sec:pfc_methods>`).

``total_element_density``
  One-dimensional array of ``float`` numbers that contains the total number density of all atoms :math:`i`, i.e. :math:`n_\mathrm{tot} = \sum_i \left( n_i + \sum_j n_j \nu_{ij} + \sum_c n_c \nu_{ic} \right)`, summed over their atomic number densities, as well as the ones contained in all other molecules/ions :math:`j` and condensates :math:`c`. This quantity is usually only a diagnostic output and not relevant for other calculations. The dimension of the array is equal to that of the input temperature and pressure vectors.

``mean_molecular_weight``
  One-dimensional array of ``float`` numbers. Contains the mean molecular weight of the mixture in units of the unified atomic mass unit. For all practical purposes, this can also be converted into units of g/mol. The dimension of the array is equal to that of the input temperature and pressure vectors.

``number_densities_cond``
  The two-dimensional array contains the fictitious number densities in of all condensate species as ``float`` numbers. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure arrays of the input structure. The second dimension refers to the number of species and has a length of ``getCondSpeciesNumber()`` (see :ref:`here<sec:pfc_methods>`).

``element_cond_degree``
  The two-dimensional array contains the degree of condensation for all elements. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of elements and has a length of ``getElementNumber()`` (see :ref:`here<sec:pfc_methods>`).

``element_conserved``
  The two-dimensional array of ``int`` numbers contains information on the state of element conservation. A value of 0 indicates that element conservation is not fulfilled, whereas a value of 1 means that the element has been conserved. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of elements and has a length of ``getElementNumber()`` (see :ref:`here<sec:pfc_methods>`).

``nb_chemistry_iterations``
  One-dimensional array of ``int`` numbers. Contains the total number of chemistry iterations that were required to solve the system for each temperature-pressure point. The dimension of the array is equal to that of the input temperature and pressure vectors.

``nb_cond_iterations``
  One-dimensional array of ``int`` numbers. Contains the total number of condensate calculation iterations that were required for each temperature-pressure point. The dimension of the array is equal to that of the input temperature and pressure vectors.

``nb_iterations``
  One-dimensional array of ``int`` numbers. Contains the total number of coupled condensation-gas phase chemistry calculation iterations that were required to solve the system for each temperature-pressure point. The dimension of the vector is equal to that of the input temperature and pressure vectors.

``fastchem_flag``
  One-dimensional array of ``int`` numbers. Contains flags that give information on potential issues of the chemistry calculation for each temperature-pressure point. The set of potential values is stated :ref:`here<sec:pfc_constants>`. A string message for each corresponding flag can also be obtained from the constant ``pyfastchem.FASTCHEM_MSG`` vector of strings, via ``pyfastchem.FASTCHEM_MSG[flag]``. The dimension of the array is equal to that of the input temperature and pressure vectors.

The output structure from the ``pyFastChem`` module, in the example here
called ``output_data``, can be defined in the following way:

.. code:: python

   output_data = pyfastchem.FastChemOutput()

The arrays of the output structure don’t need to be pre-allocated. This
will be done internally within ``FastChem`` when running the chemistry
calculations. If the arrays already contain data, their contents will be
overwritten. The arrays from the output structure can also be easily
converted to more practical NumPy arrays by using, for example:

.. code:: python

   number_densities = np.array(output_data.number_densities)

.. _sec:pfc_methods:

pyFastChem functions
~~~~~~~~~~~~~~~~~~~~~~~~

| The ``pyFastChem`` object returned from ``pyfastchem.FastChem()`` has
  several methods that allow to interact with ``FastChem``. These
  methods are equivalent to those of the ``C++`` object class discussed
  :ref:`here<sec:fastchem_methods>`. 
  
| However, note that some methods listed below officially require an ``unsigned int`` 
  parameter variable in their original C++ version. Since this data type doesn’t exist 
  in Python, PyBind11 will convert the supplied ``int`` value to its unsigned integer 
  version for C++. Even though the parameter is defined as an ``int`` value for Python, 
  only positive numbers, including 0, are accepted as valid input. Using a negative 
  value will result in an error message from PyBind11.

.. code:: python

  int calcDensities(pyfastchem.FastChemInput() input, 
                    pyfastchem.FastChemOutput() output)

..

  Starts a chemistry calculation with the provided ``pyfastchem.FastChemInput()`` and ``pyfastchem.FastChemOutput()`` structures. Returns an ``int`` value that represents the highest value from the flag vector within the ``pyfastchem.FastChemOutput()`` structure.

.. code:: python

  setParameter(str param_name, 
               param_type param_value)

..

  Sets an internal ``FastChem`` parameter. Depending on the parameter, the variable type ``param_type`` can either be an ``int``, a ``bool``, or a ``float`` value. A list of parameters and their types can be found :ref:`here<sect:fc_param>`. The ``C++``\ ``double`` types listed there should be replaced by Python ``float`` values, while ``unsigned int`` are used as ``int``. ``PyBind11`` often converts the Python ``bool`` type to an integer value rather than a ``C++``\ ``bool`` type. Setting a boolean parameter will then result in an error message. In such a case, instead of using a simple ``True`` as parameter value, an explicit conversion has to be done instead, for example via ``np.bool_(True)``.   
  
.. code:: python

  int getGasSpeciesNumber()

..

  Returns the total number of gas phase species (atoms, ions, molecules) as ``int`` value.

.. code:: python

  int getElementNumber()

..

  Returns the total number of elements as ``int`` value.

.. code:: python

  int getMoleculeNumber()

.. 

  Returns the total number of molecules and ions (anything other than elements) as ``int`` value.

.. code:: python

  int getCondSpeciesNumber()

..

  Returns the total number of condensate species as ``int`` value.

.. code:: python

  str getGasSpeciesName(int species_index)

..

  Returns the name of a gas phase species with ``int`` index ``species_index`` as ``str``; returns empty string if species does not exist.

.. code:: python

  str getGasSpeciesSymbol(int species_index)

..

  Returns the symbol of an element or the formula of a molecule/ion with ``int`` index\ ``species_index`` as ``str``; returns empty string if species does not exist

.. code:: python

  int getGasSpeciesIndex(str symbol)

..

  Returns the index of a gas phase species (element/molecule/ion) with ``str`` symbol/formula ``symbol`` as ``int``; returns the constant ``pyfastchem.FASTCHEM_UNKOWN_SPECIES`` if species does not exist.

.. code:: python

  str getElementName(int species_index)

.. 

  Returns the name of an element with ``int`` index\ ``species_index`` as ``str``; returns empty string if species does not exist.

.. code:: python
  
  str getElementSymbol(int species_index)

..

  Returns the symbol of an element with ``int`` index\ ``species_index`` as ``str``; returns empty string if species does not exist

.. code:: python

  int getElementIndex(str symbol)

.. 

  Returns the index of an element with ``str`` symbol/formula ``symbol`` as ``int``; returns the constant ``pyfastchem.FASTCHEM_UNKOWN_SPECIES`` if species does not exist.

.. code:: python

  str getCondSpeciesName(int species_index)

..

  Returns the name of a condensate species with ``int`` index\ ``species_index`` as ``str``; returns empty string if species does not exist.

.. code:: python

  str getCondSpeciesSymbol(int species_index)

..
  Returns the formula of a condensate with ``int`` index\ ``species_index`` as ``str``; returns empty string if species does not exist

.. code:: python

  int getCondSpeciesIndex(str symbol)

..

  Returns the index of a condensate species with ``str`` symbol/formula ``symbol`` as ``int``; returns the constant ``pyfastchem.FASTCHEM_UNKOWN_SPECIES`` if species does not exist.

.. code:: python

  float getElementAbundance(int species_index)

..

Returns the abundance of an element with ``int`` index\ ``species_index`` as ``float``; returns 0 if the element does not exist

.. code:: python

  float [] getElementAbundance()

..

  Returns the abundances of all elements as an array of ``float`` values; array has a length of ``getElementNumber()``

.. code:: python

  setElementAbundances(float [] abundances)

..

  Sets the abundances of all elements; the abundances are supplied as an array of ``float`` values, where the array has to have a size of ``getElementNumber()``; if this is not the case, ``FastChem`` will print an error message and leave the element abundances unchanged
  
.. code:: python

  float FastChem.getGasSpeciesWeight(int species_index)

.. 

  Returns the weight of a gas phase species with ``int`` index\ ``species_index`` as ``float``; returns 0 if species does not exist; for an element this refers to the atomic weight

.. code:: python

  float FastChem.getElementWeight(int species_index)

..  

  Returns the atomic weight of an element with ``int`` index\ ``species_index`` as ``float``; returns 0 if species does not exist
  
.. code:: python

  float FastChem.getCondSpeciesWeight(int species_index)

.. 

  Returns the weight of a condensate species with ``int`` index\ ``species_index`` as ``float``; returns 0 if species does not exist

.. code:: python

  setVerboseLevel(int level)

..

  Sets the verbose level of ``FastChem``, i.e. the amount of text output in the terminal. A value of 0 will result in ``FastChem`` being almost silent, whereas a value of 4 would provide a lot of debug output. A value larger than 4 will be interpreted as 4. This value will overwrite the one from the ``FastChem`` config file.

