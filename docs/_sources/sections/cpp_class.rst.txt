Detailed C++ object class description
=========================================

.. _sec:fc_class:


``FastChem`` has been written in an object oriented way, split across
several different object classes. The entire source code of ``FastChem``
is contained in the folder ``fastchem_src/``. For including ``FastChem``
in another ``C++`` project, only adding the main fastchem header file
``fastchem.h`` is required. All ``FastChem`` code is encapsulated in its
own namespace called ``fastchem`` to avoid clashing with other
libraries.

Some comments on coding conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The entire ``FastChem`` code has been programmed using specific
conventions that make it easy to recognise and differentiate class,
method and variable names.

Class and structure names are always capitalised, for example

``class FastChem{...`` or ``struct Molecule{...``

| If the name is a compound noun, each noun is capitalised separately,
  e.g. ``struct ChemicalElementData{...`` No separators like \_ are used
  for class or structure names.
| Class methods and functions in general always start with a lowercase
  letter. If the name is a compound noun, the start of every other noun
  is capitalised and no separator is used, for example:
  ``FastChem.setElementAbundances(...)`` or
  ``Molecule.calcMassActionConstant(...)``

Variable names are always written in lowercase and compound nouns are
separated by a \_. For example: ``FastChem.element_calculation_order``
or ``Element.molecule_list`` The only exceptions are global constants.
They contain only capitalised letters, e.g.:
``constexpr unsigned int FASTCHEM_UNKNOWN_SPECIES`` or
``constexpr double CONST_AMU``

FastChem object class
~~~~~~~~~~~~~~~~~~~~~~~~~

The entire ``FastChem`` model is encapsulated in an object class called
``FastChem`` that is defined in the header file ``fastchem.h``. The
object class is programmed as a template that can be used in either
``double`` or ``long double`` precision. When creating an object of this
class by calling a corresponding constructor, one therefore has to
specify which of the two versions should be used:

.. code:: cpp

  FastChem<long double> fastchem("model_parameter_file.dat", 1);
  FastChem<double> fastchem("model_parameter_file.dat", 1); 

The ``long double`` version has a slightly higher computational overhead and
larger memory requirements than the ``double`` one. On the other hand,
it offers a higher numerical precision which is especially important
when dealing with chemical systems where the mass action constants and
number densities can vary by many orders of magnitudes. We strongly
suggest to always use the ``long double`` version despite the additional
computational overhead. In fact, in our experience the increased
numerical precision of ``long double`` can effectively lead to a smaller
number of iterations.

.. _sec:fastchem_constants:

FastChem constants
~~~~~~~~~~~~~~~~~~

| The ``FastChem`` namespace ``fastchem`` contains a number of constants that are all defined in the file ``fastchem_constants.h``. This includes the constant

.. code:: cpp

  constexpr unsigned int fastchem::FASTCHEM_UNKNOWN_SPECIES 
  
that is returned by some ``FastChem`` methods when a chemical species is not found. The chemistry calculation will return several output flags that are also defined in this file. This includes the following constants:

``constexpr unsigned int fastchem::FASTCHEM_SUCCESS``
  Indicates that the calculation has been successful, i.e. that the chemistry iterations converged.

``constexpr unsigned int fastchem::FASTCHEM_NO_CONVERGENCE``
  Indicates that the calculation was not successful, i.e. that the chemistry did not converge within the allowed maximum number of iterations steps given in the config file or set manually via ``FastChem.setMaxChemistryIter`` (see :ref:`here<sec:fastchem_methods>`). One way to solve such a problem is to increase the maximum number of iteration steps.

``constexpr unsigned int fastchem::FASTCHEM_INITIALIZATION_FAILED``
  Indicates that something went wrong during reading one of the input files. To find the source of the problem, one can set the verbose level in the config file or manually via ``FastChem.setVerboseLevel`` (see :ref:`here<sec:fastchem_methods>`) to a higher value and look at the terminal output.

``constexpr unsigned int fastchem::FASTCHEM_IS_BUSY``
  The chemistry calculations of ``FastChem`` can only be called once for each object class instance. Attempting to start a new calculation while another is still running, for example via OpenMP, will result in ``FastChem`` returning this flag.

``constexpr unsigned int fastchem::FASTCHEM_WRONG_INPUT_VALUES``
  ``FastChem`` returns this flag if some input values are wrong. Currently, this refers to the temperature and pressure vectors in the input structure not having the same size (see :ref:`here<sec:fastchem_input_struct>` for details on the input structure).

``constexpr unsigned int fastchem::FASTCHEM_PHASE_RULE_VIOLATION``
  ``FastChem`` returns this flag if condensation is used and the system violates the phase rule. This happens when the number of elements contained in condensates equals the total number of elements. In this case, the gas phase lacks a degree of freedom to yield the correct gas pressure. Such a system cannot be solved as there has always to be at least one incondensable element in the gas phase (see the section about the phase rule in Paper III).

In addition to these flags, ``fastchem_constants.h`` also includes a
constant string vector

.. code:: cpp

  const std::vector<std::string> fastchem::FASTCHEM_MSG
  
that contains string expressions for each of these flags. Using this vector with any
of the aforementioned flags ``fastchem::FASTCHEM_MSG[flag]`` returns a
string with a description of the corresponding flag’s meaning. For
example, 

.. code:: cpp

  fastchem::FASTCHEM_MSG[fastchem::FASTCHEM_NO_CONVERGENCE]

will return the string ``"convergence failed"``.

FastChem constructor
~~~~~~~~~~~~~~~~~~~~

| Since ``FastChem`` is written as an object class, an instance of that
  class (i.e. an object) needs to be created before ``FastChem`` can be
  used. This is done by calling the constructor of the ``FastChem``
  class. There are three primary ways to call the constructor and create
  an object.
  
.. code:: cpp

  FastChem(const std::string& model_parameter_file, 
           const unsigned int verbose_level_init)

..

  This constructor requires two parameters: the location of the parameter file, described :ref:`here<sec:fc_parameter_file>`, as well as the initial verbose, i.e. the amount of debug output in the terminal window. All main options and parameters will be read from the parameter file, but can be changed later by using the appropriate methods described :ref:`here<sec:fastchem_methods>`.

.. code:: cpp

   FastChem(const std::string& element_abundances_file, 
            const std::string& gas_species_data_file, 
            const unsigned int verbose_level)

..

  This constructor requires three parameters: the locations of the element abundance and gas phase species data files, as well as the verbose level. All other options and parameters within ``FastChem`` will be set to their default values but can be later changed by using the appropriate methods described :ref:`here<sec:fastchem_methods>`. The default maximum number of chemistry iterations is 3000, the number of Newton, bisection and Nelder-Mead method iterations is 3000, and the default accuracy of the of Newton method and the chemistry iterations is set to :math:`10^{-4}`. This constructor will not read in any condensate data. Trying to use an object created via this method for a calculation using condensation will result in an error message.

.. code:: cpp

   FastChem(const std::string& element_abundances_file, 
            const std::string& gas_species_data_file, 
            const std::string& cond_species_data_file,
            const unsigned int verbose_level)

..

  This constructor requires four parameters: the locations of the element abundance and gas phase species data files, the condensate data file, as well as the verbose level. All other options and parameters within ``FastChem`` will be set to their default values but can be later changed by using the appropriate methods described :ref:`here<sec:fastchem_methods>`. The default maximum number of chemistry iterations is 3000, the number of Newton, bisection and Nelder-Mead method iterations is 3000, and the default accuracy of the of Newton method and the chemistry iterations is set to :math:`10^{-4}`. Note that instead of a location for the condensate data, a string containing ``"none"`` can be used here as well. In that case, no condensate data will be read in and trying to use the object for a calculation using condensation will result in an error message.

A fourth way to create a ``FastChem`` object is to make a copy of an
existing one. ``FastChem`` contains an internal copy constructor that
manages the copy of all the object class’ data structures. Assuming that
``fastchem_a`` is a valid object instance of the ``FastChem`` class, a
second object, say ``fastchem_b``, can simply be created by using

.. code:: cpp

   fastchem::FastChem fastchem_b(fastchem_a);

In this example, ``fastchem_b`` is a direct copy of ``fastchem_a``, i.e.
all parameters, options, and species & element data structures are
identical. After the creation of ``fastchem_b``, both objects can be
used independently from each other and can even be run at the same time.

Input and output structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the chemistry calculation of FastChem,

.. code:: cpp

  FastChem.calcDensities(FastChemInput, FastChemOutput)

is called, input and output structures are required. Their definitions can be found
in the source file ``fastchem_src/input_output_struct.h``

.. _sec:fastchem_input_struct:

Input structure
^^^^^^^^^^^^^^^

The input structure is defined as follows:

.. code:: c++

   struct FastChemInput
     {
       std::vector<double> temperature; 
       std::vector<double> pressure;
       
       bool equilibrium_condensation = false;
       bool rainout_condensation = false;
     };

It contains vectors for the temperatures (in K) and pressures (in bar)
that the chemical composition should be calculated for. Both vectors
need to have the same length. Otherwise, ``FastChem.calcDensities`` will
return the constant ``fastchem::FASTCHEM_WRONG_INPUT_VALUES``. Note that
even if you want to run the chemistry for only a single temperature and
pressure point, you still need to provide the input in vectorial form.

The two boolean variables enable the calculation of either equilibrium
condensation or the rainout approximation. By default, both are set to
``false``. Note that when the flag ``rainout_condensation`` is set to
``true``, the value of the parameter ``equilibrium_condensation`` is
ignored.

Output structure
^^^^^^^^^^^^^^^^

The outout structure is defined as

.. code:: c++

   struct FastChemOutput
     {
       std::vector<std::vector<double>> number_densities;
       std::vector<double> total_element_density;
       std::vector<double> mean_molecular_weight;
       
       std::vector<std::vector<double>> number_densities_cond;
       std::vector<std::vector<double>> element_cond_degree;
       
       std::vector<std::vector<unsigned int>> element_conserved;
       std::vector<unsigned int> nb_chemistry_iterations;
       std::vector<unsigned int> nb_cond_iterations;
       std::vector<unsigned int> nb_iterations;
       std::vector<unsigned int> fastchem_flag;
     };

It has the following variables:

``std::vector<std::vector<double>> number_densities``
  The two-dimensional array contains the number densities in of all gas phase species (elements, molecules, ions). The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of species and has a length of ``FastChem.getGasSpeciesNumber()`` (see :ref:`here<sec:fastchem_methods>`).

``std::vector<double> total_element_density``
  Contains the total number density of all atoms :math:`i`, i.e. :math:`n_\mathrm{tot} = \sum_i \left( n_i + \sum_j n_j \nu_{ij} + \sum_c n_c \nu_{ic} \right)`, summed over their atomic number densities, as well as the ones contained in all other molecules/ions :math:`j` as well as condensate species :math:`c`. This quantity is usually only a diagnostic output and not relevant for other calculations. The dimension of the vector is equal to that of the input temperature and pressure vectors.

``std::vector<double> mean_molecular_weight``
  Contains the mean molecular weight of the mixture in units of the unified atomic mass unit. For all practical purposes, this can also be converted into units of g/mol. The dimension of the vector is equal to that of the input temperature and pressure vectors.

``std::vector<std::vector<double>> number_densities_cond``
  The two-dimensional array contains the fictitious number densities in of all condensate species. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of species and has a length of ``FastChem.getCondSpeciesNumber()`` (see :ref:`here<sec:fastchem_methods>`).

``std::vector<std::vector<double>> element_cond_degree``
  The two-dimensional array contains the degree of condensation for all elements. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of elements and has a length of ``FastChem.getElementNumber()`` (see :ref:`here<sec:fastchem_methods>`).

``std::vector<std::vector<unsigned int>> element_conserved``
  The two-dimensional array contains information on the state of element conservation. A value of 0 indicates that element conservation is violated, whereas a value of 1 means that the element has been conserved. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of elements and has a length of ``FastChem.getElementNumber()`` (see :ref:`here<sec:fastchem_methods>`).

``std::vector<unsigned int> nb_chemistry_iterations``
  Contains the total number of chemistry iterations that were required to solve the system for each temperature-pressure point. The dimension of the vector is equal to that of the input temperature and pressure vectors.

``std::vector<unsigned int> nb_iterations``
  Contains the total number of coupled condensation-gas phase chemistry calculation iterations that were required to solve the system for each temperature-pressure point. The dimension of the vector is equal to that of the input temperature and pressure vectors.

``std::vector<unsigned int> nb_chemistry_iterations``
  Contains the total number of chemistry iterations that were required to solve the system for each temperature-pressure point. The dimension of the vector is equal to that of the input temperature and pressure vectors.

``std::vector<unsigned int> fastchem_flag``
  Contains flags that give information on potential issues of the chemistry calculation for each temperature-pressure point. The set of potential values is stated :ref:`here<sec:fastchem_constants>`. A string message for each corresponding flag can also be obtained from the constant ``fastchem::FASTCHEM_MSG`` vector of strings, via ``fastchem::FASTCHEM_MSG[flag]``. The dimension of the vector is equal to that of the input temperature and pressure vectors.

The vectors of the output structure don’t need to be pre-allocated. This
will be done internally within ``FastChem`` when running the chemistry
calculations. If the vectors already contain data, their contents will
be overwritten.

.. _sec:fastchem_methods:

Public methods of the FastChem object class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: c++

  unsigned int FastChem.calcDensities(FastChemInput input, 
                                      FastChemOutput output)

..
  
  Starts a chemistry calculation with the provided ``FastChemInput`` and ``FastChemOutput`` structs. Returns an ``unsigned int`` that represents the highest value from the flag vector within the ``FastChemOutput`` struct.
  
.. code:: c++

  void FastChem.setParameter(std::string param_name, 
                             param_type param_value)

..

  Sets an internal ``FastChem`` parameter with the name ``param_name``. Depending on the parameter, the variable type ``param_type`` can either be an ``unsigned int``, a ``bool``, or a ``double`` value. A list of parameters and their types can be found in the next section.

.. code:: c++

  unsigned int FastChem.getGasSpeciesNumber()

..

  Returns the total number of gas phase species (atoms, ions, molecules) as ``unsigned int``

.. code:: c++

  unsigned int FastChem.getElementNumber()

..

  Returns the total number of elements as ``unsigned int``

.. code:: c++

  unsigned int FastChem.getMoleculeNumber()

..

  Returns the total number of molecules and ions (anything other than elements) as ``unsigned int``

.. code:: c++

  unsigned int FastChem.getCondSpeciesNumber()

..

  Returns the total number of condensate species as ``unsigned int``

.. code:: c++

  std::string FastChem.getGasSpeciesName(unsigned int species_index)

.. 

  Returns the name of a gas phase species with index ``species_index`` as ``std::string``; returns empty string if species does not exist

.. code:: c++

  std::string FastChem.getGasSpeciesSymbol(unsigned int species_index)

..

  Returns the symbol of an element or the formula of a molecule/ion with index ``species_index`` as ``std::string``; returns empty string if species does not exist

.. code:: c++

  unsigned int FastChem.getGasSpeciesIndex(std::string symbol)

.. 

  Returns the index of a species (element/molecule/ion) with symbol/formula ``symbol`` as ``unsigned int``; returns the constant ``fastchem::FASTCHEM_UNKOWN_SPECIES`` if species does not exist

  
.. code:: c++

  std::string FastChem.getElementName(unsigned int species_index)

..
  
  Returns the name of an element with index ``species_index`` as ``std::string``; returns empty string if species does not exist

.. code:: c++

  std::string FastChem.getElementSymbol(unsigned int species_index)

..
  
  Returns the symbol of an element with index ``species_index`` as ``std::string``; returns empty string if species does not exist

.. code:: c++

  unsigned int FastChem.getElementIndex(std::string symbol)

..

  Returns the index of an element with symbol ``symbol`` as ``unsigned int``; returns the constant ``fastchem::FASTCHEM_UNKOWN_SPECIES`` if species does not exist

.. code:: c++

  std::string FastChem.getCondSpeciesName(unsigned int species_index)

..

  Returns the name of a condensate species with index ``species_index`` as ``std::string``; returns empty string if species does not exist

.. code:: c++

  std::string FastChem.getCondSpeciesSymbol(unsigned int species_index)

..

  Returns the formula of a condensate species with index ``species_index`` as ``std::string``; returns empty string if species does not exist

.. code:: c++
  
  unsigned int FastChem.getCondSpeciesIndex(std::string symbol)

..

  Returns the index of a condensate species formula ``symbol`` as ``unsigned int``; returns the constant ``fastchem::FASTCHEM_UNKOWN_SPECIES`` if species does not exist

.. code:: c++

  double FastChem.getElementAbundance(unsigned int species_index)

..

  Returns the abundance of an element with index ``species_index`` as ``double``; returns 0 if element does not exist

.. code:: c++

  std::vector<double> FastChem.getElementAbundance()

..

  Returns the abundances of all elements as a vector of ``double``; vector has a length of ``FastChem.getElementNumber()``

.. code:: c++
  
  double FastChem.getGasSpeciesWeight(unsigned int species_index)

..

  Returns the weight of a gas phase species (element/molecule/ion) with index ``species_index`` as ``double``; returns 0 if species does not exist; for an element this refers to the atomic weight

.. code:: c++

  double FastChem.getElememtWeight(unsigned int species_index)

..

  Returns the atomic weight of an element with index ``species_index`` as ``double``; returns 0 if the element does not exist

.. code:: c++
  
  double FastChem.getCondSpeciesWeight(unsigned int species_index)

..

Returns the weight of a condensate species with index ``species_index`` as ``double``; returns 0 if species does not exist

.. code:: c++

  void FastChem.setElementAbundances(std::vector<double> abundances)

..

  Sets the abundances of all elements; the abundances are supplied as ``std::vector<double>``, where the vector has to have a size of ``FastChem.getElementNumber()``; if this is not the case, ``FastChem`` will print an error message and leave the element abundances unchanged

.. code:: c++

  void FastChem.setVerboseLevel(unsigned int level)

..

  Sets the verbose level of ``FastChem``, i.e. the amount of text output in the terminal. A value of 0 will result in ``FastChem`` being almost silent, whereas a value of 4 would provide a lot of debug output. A value larger than 4 will be interpreted as 4. This value will overwrite the one from the ``FastChem`` config file.

