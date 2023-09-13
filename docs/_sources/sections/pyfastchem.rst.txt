Running pyFastChem
==================

In addition to the ``C++`` executable, we provide several Python scripts
that can run the ``FastChem`` code through its Python interface
``pyFastChem``. The sample scripts can be found within the ``python/``
folder. These sample scripts show different use cases and can be used as
a basis for your own ``FastChem`` Python scripts.

.. _sec:pyfc_examples:

Provided Python examples
~~~~~~~~~~~~~~~~~~~~~~~~

| Currently, we provide the following examples:

``fastchem.py``
  Runs a simple ``FastChem`` calculation on a temperature-pressure structure defined within the script, writes output files, and creates a plot with selected species.

``fastchem_c_o.py``
  Runs a ``FastChem`` calculation on a temperature-pressure structure defined within the script and for a range of different C/O ratios. It will write output files, and create a plot with selected species.

``fastchem_metallicity.py``
  Runs a ``FastChem`` calculation on a temperature-pressure structure defined within the script and for a range of different metallicity factors. It will write output files, and create a plot with selected species.

``fastchem_cond.py``
  Runs a ``FastChem`` calculation on a read-in temperature-pressure structure of a brown dwarf with equilibrium condensation turned on. The rainout approximation can optionally be turned on as well. It will write output files, and create a plot with selected species.

``fastchem_cond_disk.py``
  Runs a ``FastChem`` calculation for the temperature-pressure structure of the midplane of a protoplanetary disk with equilibrium condensation. It will write output files, and create a plot with selected species. **Warning**: Due to the very low temperatures in the outer part of the disk, this calculation might take quite some time depending on your computer.

Note that the scripts should be executed from within the Python folder
since all file paths in the scripts are given relative to this
directory. These files can be used as templates to create your own
special Python scripts to run ``pyFastChem``. The following section
provides some details on the steps required to run ``FastChem`` from
within Python. A more detailed overview of all the methods and variables
available within ``pyFastChem`` can be found :ref:`here<sec:pyfc_details>`.

Detailed steps for running pyFastChem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| As a first step, we need to import the ``pyFastChem`` module:

.. code:: python

   import pyfastchem

| This will import the module compiled by ``PyBind11``. Next, we have to
  create a ``FastChem`` object (here named ``fastchem``) with a
  corresponding constructor ``pyfastchem.FastChem`` provided by
  ``pyFastChem``. There are three different possible versions of the
  constructor. F

If we’re only interested in a pure gas phase chemistry calculation, we
can use

.. code:: python

   fastchem = pyfastchem.FastChem('input/element_abundances/asplund_2009.dat', 
                                  'input/logK/logK.dat', 
                                  1)

| This constructor requires three different arguments: the location of
  the element abundance file, the location of the file with the
  equilibrium constants, and the verbose level.
| For a ``FastChem`` object that also includes condensed species, we
  have:

.. code:: python

   fastchem = pyfastchem.FastChem('input/element_abundances/asplund_2009.dat', 
                                  'input/logK/logK.dat', 
                                  'input/logK/logK_condensates.dat', 
                                  1)

where the additional argument is the location of the input file with the
equilibrium constants of the condensed species. It is also possible to
use ``'none'`` for this argument:

.. code:: python

   fastchem = pyfastchem.FastChem('input/element_abundances/asplund_2009.dat', 
                                  'input/logK/logK.dat', 
                                  'none', 
                                   1)

| in which case no condensate data will be read in and this constructor
  behaves like the first one by only allowing the calculation of the
  pure gas phase composition.
| Alternatively, a ``FastChem`` object can also be created via

.. code:: python

   fastchem = pyfastchem.FastChem('input/parameters_py.dat', 1)

| where the first argument is the location of the parameter file and the
  second one the initial verbose level. The latter one will later be
  replaced by the corresponding value read in from the parameter file.
  The structure of this parameter file is discussed :ref:`here<sec:fc_parameter_file>`.

| Creating a ``FastChem`` object with the first two methods will set
  internal parameters to their default values. The maximum number of
  chemistry iterations will be 3000, the number of Newton, bisection and
  Nelder-Mead method iterations is 3000, and the accuracy of the of
  Newton method and the chemistry iterations is set to :math:`10^{-4}`.
  All of these values can, however, be adjusted during runtime by using
  the methods listed :ref:`here<sec:pfc_methods>`.
  
| Next, we need to create the input and output structures used by
  ``pyFastChem``:

.. code:: python

   input_data = pyfastchem.FastChemInput()
   output_data = pyfastchem.FastChemOutput()

| Details on these structures can be found :ref:`here<sec:pfc_input_ouput_struct>`. 
  The input structure contains
  the temperature (in K) and pressure (in bar) arrays that the chemistry
  should be calculated for. They can be set, for example, by:

.. code:: python

   input_data.temperature = temperature
   input_data.pressure = pressure

| where ``temperature`` and ``pressure`` are standard Python lists or
  NumPy arrays. Both arrays need to have the same length. The input
  structure also contains two boolean flags that enable the calculation
  of the condensed phase:
  
| ``input_data.equilibrium_condensation``
| ``input_data.rainout_condensation``

| Setting the first flag to ``True`` will calculate the chemical
  composition assuming equilibrium condensation for each
  temperature-pressure point of the input structure separately. Setting
  the rainout condensation flag to ``True`` enables the calculation
  using the rainout approximation. Details on this can be found in
  `Kitzmann, Stock & Patzer (2023) <http://adsabs.harvard.edu/abs/2023arXiv230902337K>`_. 
  Note that if the rainout flag
  is set to ``True``, the value of the equilibrium condensation flag is
  ignored. By default, both flags are set to ``False``.
  
| With the input structure properly set up, we can now run the actual
  ``FastChem`` calculation by calling the ``calcDensities`` method:

.. code:: python

   fastchem_flag = fastchem.calcDensities(input_data, output_data)

| This method returns an integer flag that describes the overall outcome
  of the calculation. A description of the different flags can be found
  :ref:`here<sec:pfc_constants>`. After calling the
  ``calcDensity`` method, the output structure will be filled with the
  corresponding output data. For example,
| ``output_data.number_densities`` will contain the number densities of
  the chemical species. This is a 2D list, where the first dimension
  refers to the temperature and pressure input arrays and the second
  dimension refers to the different chemical species. The list can be
  easily converted into a NumPy array via:

.. code:: python

   number_densities = np.array(output_data.number_densities)

| The Python directory of the ``FastChem`` repository also contains
  functions that save the output into files, identical to those from the
  ``C++`` version. They can be called by:

.. code:: python

   saveChemistryOutput(output_dir + '/chemistry.dat', 
                       temperature,
                       pressure, 
                       output_data.total_element_density, 
                       output_data.mean_molecular_weight,  
                       output_data.number_densities, 
                       fastchem)
.. code:: python
   
   saveCondOutput(output_dir + '/condensates.dat', 
                  temperature, 
                  pressure, 
                  output_data.element_cond_degree, 
                  output_data.number_densities_cond, 
                  fastchem)

.. code:: python

   saveMonitorOutput(output_dir + '/monitor.dat', 
                     temperature, 
                     pressure, 
                     output_data.element_conserved, 
                     output_data.fastchem_flag, 
                     output_data.nb_chemistry_iterations, 
                     output_data.total_element_density, 
                     output_data.mean_molecular_weight, 
                     fastchem)

| A more detailed description of the output functions can be found in
  the next section.

Output functions of pyFastChem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| The Python directory includes several scripts that can save the
  ``FastChem`` chemistry and monitor output in either text or binary
  data files. All these functions are located within the file
  ``save_output.py``. Examples of their usage can be found in the three
  Python scripts discussed above.

Chemistry output scripts
^^^^^^^^^^^^^^^^^^^^^^^^

``save_output.py`` contains two functions for the general chemistry
output. The first, ``saveChemistryOutput``, saves the results in a text
file that is identical to the one of the ``C++`` version. If the
chemistry is calculated for a larger number of pressure and temperature
points, the output can become quite large. Saving these results into a
simple text file can, therefore, take a very long time - in extreme
cases even longer than the calculation itself.

Therefore, we provide an alternative function
``saveChemistryOutputPandas`` that saves the output in a pandas
``DataFrame`` format into a pickle file. Since this is a binary format,
saving a large output is substantially faster than the corresponding
ASCII text version.

| The function for saving the output as a normal text file is

.. code:: python

   saveChemistryOutput(file_path, 
                       temperature, 
                       pressure, 
                       total_element_density, 
                       mean_molecular_weight, 
                       number_densities, 
                       fastchem, 
                       output_species=None, 
                       additional_columns=None, 
                       additional_columns_desc=None)

with the following arguments:

``file_path``
  Contains the path to the output file as a ``str`` variable.

``temperature, pressure``
  Arrays of ``float`` values with the temperature and pressure structure the chemistry has been calculated for.

``total_element_density``
  ``float`` array of the total number density of all atoms :math:`i`, i.e. :math:`n_\mathrm{tot} = \sum_i \left( n_i + \sum_j n_j \nu_{ij} + \sum_c n_c \nu_{ic}\right)`, summed over their atomic number densities, as well as the ones contained in all other molecules/ions :math:`j` and condensate species :math:`c`. This quantity is usually only a diagnostic output and not relevant for other calculations. The dimension of the array is equal to that of the temperature and pressure vectors.

``mean_molecular_weight``
  ``float`` array of the computed mean molecular weight. The dimension of the array is equal to that of the temperature and pressure vectors.

``number_densities``
  Two-dimensional ``float`` array of the number densities. The first dimension of the array refers to the temperature and pressure input arrays, while the second dimension describes the different chemical species.

``fastchem``
  An object of the ``pyFastChem`` class that has been used to calculate the chemistry.

``output_species=None``
  Optional parameter. Is an array of ``str`` values that contains the chemical symbols of species the chemistry output file should be saved for. Without this optional parameter, the output function will by default save all species. The symbols have to match the ones used in the ``FastChem`` input file for the equilibrium constants. For the standard files supplied with ``FastChem``, the Hill notation, therefore, needs to be used here.

``additional_columns=None, additional_columns_desc=None``
  Optional parameters. Sometimes, ``FastChem`` calculations are not iterated only over temperature or pressure but also other variables, such as the metallicity or C/O ratio. The output function therefore contains these optional parameters that allow to print additional columns in the output file. The first parameter ``additional_columns`` is an :math:`N\times N_\mathrm{tp}`-dimensional array of ``float`` values, where the first dimension refers to the number of additional columns and the second dimension has to be equal to the dimensions of the temperature and pressure arrays.

  The second optional parameter ``additional_columns_desc`` contains an array of ``str`` values with the header descriptions of the additional columns. The dimension has to be equal to number of additional columns. If this is not the case, or if the parameter is missing entirely, the columns will be labelled ``unk`` instead.

| All of these function arguments, except for the optional parameters,
  are contained within the input and output structures of
  ``pyFastChem``, discussed :ref:`here<sec:pfc_input_ouput_struct>`.

| Saving the chemistry output with the panda ``DataFrame`` format in a
  pickle file is possible via the function:

.. code:: python

  saveChemistryOutputPandas(file_path,
                            temperature, 
                            pressure, 
                            total_element_density, 
                            mean_molecular_weight, 
                            number_densities, 
                            fastchem, 
                            output_species=None, 
                            additional_columns=None, 
                            additional_columns_desc=None)

All arguments are identical to those of the previous ASCII output
function. The saved panda DataFrame contains the same columns and
headers as the simple text output.

Condensate output script
^^^^^^^^^^^^^^^^^^^^^^^^

| The two condensate output scripts are almost identical to the gas
  phase chemistry one described above. The first script saves the output
  into a simple text file:

.. code:: python

   saveCondOutput(file_path, 
                  temperature, 
                  pressure, 
                  element_cond_degree, 
                  number_densities, 
                  fastchem, 
                  output_species=None, 
                  additional_columns=None, 
                  additional_columns_desc=None)

| It has the following arguments:

``file_path``
  Contains the path to the output file as a ``str`` variable.

``temperature, pressure``
  Arrays of ``float`` values with the temperature and pressure structure the chemistry has been calculated for.

``element_cond_degree``
  Two-dimensional ``float`` array of the degrees of condensation for all elements. The first dimension of the array refers to the temperature and pressure input arrays, while the second dimension describes the different elements.

``number_densities``
  Two-dimensional ``float`` array of the (fictitious) condensate number densities. The first dimension of the array refers to the temperature and pressure input arrays, while the second dimension describes the different condensate species.

``fastchem``
  An object of the ``pyFastChem`` class that has been used to calculate the chemistry.

``output_species=None``
  Optional parameter. Is an array of ``str`` values that contains the chemical symbols of condensates the chemistry output file should be saved for. Without this optional parameter, the output function will by default save all species. The symbols have to match the ones used in the ``FastChem`` input file for the equilibrium constants.

``additional_columns=None, additional_columns_desc=None``
  Optional parameters. Sometimes, ``FastChem`` calculations are not iterated only over temperature or pressure but also other variables, such as the metallicity or C/O ratio. The output function therefore contains these optional parameters that allow to print additional columns in the output file. The first parameter ``additional_columns`` is an :math:`N\times N_\mathrm{tp}`-dimensional array of ``float`` values, where the first dimension refers to the number of additional columns and the second dimension has to be equal to the dimensions of the temperature and pressure arrays.

  The second optional parameter ``additional_columns_desc`` contains an array of ``str`` values with the header descriptions of the additional columns. The dimension has to be equal to number of additional columns. If this is not the case, or if the parameter is missing entirely, the columns will be labelled ``unk`` instead.

| All of these function arguments, except for the optional parameters,
  are contained within the input and output structures of
  ``pyFastChem``, discussed :ref:`here<sec:pfc_input_ouput_struct>`.
| Saving the condensate output with the panda ``DataFrame`` format in a
  pickle file is possible via the function:

.. code:: python
  
  saveCondOutputPandas(file_path, 
                       temperature, 
                       pressure, 
                       element_cond_degree, 
                       number_densities_cond,
                       fastchem, 
                       output_species=None, 
                       additional_columns=None, 
                       additional_columns_desc=None)

All arguments are identical to those of the previous ASCII output
function. The saved panda DataFrame contains the same columns and
headers as the simple text output.

Monitor output scripts
^^^^^^^^^^^^^^^^^^^^^^

``save_output.py`` also contains two functions for the ``FastChem``
monitor output. The first, ``savMonitorOutput``, saves the debug output
in a text file that is identical to the one of the ``C++`` version. Just
like for the chemistry output, saving the results for a large number of
calculations can be quite slow. Therefore, we also provide an
alternative function ``saveMonitorOutputPandas`` that saves the output
as a pandas ``DataFrame`` format into a pickle file.

| The function for saving the output as a normal text file is

.. code:: python

  saveMonitorOutput(file_path, 
                    temperature, 
                    pressure, 
                    element_conserved, 
                    fastchem_flags, 
                    nb_iterations, 
                    nb_chemistry_iterations, 
                    nb_condensation_iterations, 
                    total_element_density, 
                    mean_molecular_weight, 
                    fastchem, 
                    additional_columns=None, 
                    additional_columns_desc=None)

with the following arguments:

``file_path``
  Contains the path to the output file as a ``str`` variable.

``temperature, pressure``
  Arrays of ``float`` values with the temperature and pressure structure the chemistry has been calculated for.

``element_conserved``
  The two-dimensional array of ``int`` numbers contains information on the state of element conservation. A value of 0 indicates that element conservation is not fulfilled, whereas a value of 1 means that the element has been conserved. The first dimension refers to the temperature-pressure grid and has the same size as the temperature and pressure vectors of the input structure. The second dimension refers to the number of elements and has a length of ``getElementNumber()`` (see :ref:`here<sec:pfc_methods>`).

``fastchem_flags``
  One-dimensional array of ``int`` numbers. Contains flags that give information on potential issues of the chemistry calculation for each temperature-pressure point. The set of potential values is stated :ref:`here<sec:pfc_constants>`. A string message for each corresponding flag can also be obtained from the constant ``pyfastchem.FASTCHEM_MSG`` vector of strings, via ``pyfastchem.FASTCHEM_MSG[flag]``. The dimension of the array is equal to that of the input temperature and pressure vectors.

``nb_iterations``
  One-dimensional array of ``int`` numbers. Contains the number of coupled chemistry-condensation iterations that were required to solve the system for each temperature-pressure point. The dimension of the array is equal to that of the input temperature and pressure vectors.

``nb_chemistry_iterations``
  One-dimensional array of ``int`` numbers. Contains the total number of chemistry iterations that were required to solve the system for each temperature-pressure point. The dimension of the array is equal to that of the input temperature and pressure vectors.

``nb_condensation_iterations``
  One-dimensional array of ``int`` numbers. Contains the total number of condensation calculation iterations that were required to solve the system for each temperature-pressure point. The dimension of the array is equal to that of the input temperature and pressure vectors.

``total_element_density``
  One-dimensional array of ``float`` numbers that contains the total number density of all atoms :math:`i`, i.e. :math:`n_\mathrm{tot} = \sum_i \left( n_i + \sum_j n_j \nu_{ij} + \sum_c n_c \nu_{ic}\right)`, summed over their atomic number densities, as well as the ones contained in all other molecules/ions :math:`j` and condensate species :math:`c`. This quantity is usually only a diagnostic output and not relevant for other calculations. The dimension of the array is equal to that of the input temperature and pressure vectors.

``mean_molecular_weight``
  One-dimensional array of ``float`` numbers. Contains the mean molecular weight of the mixture in units of the unified atomic mass unit. For all practical purposes, this can also be converted into units of g/mol. The dimension of the array is equal to that of the input temperature and pressure vectors.

``fastchem``
  An object of the ``pyFastChem`` class that has been used to calculate the chemistry.

``additional_columns=None, additional_columns_desc=None``
  Optional parameters. Sometimes, ``FastChem`` calculations are not iterated only over temperature or pressure but also other variables, such as the metallicity or C/O ratio. The output function therefore contains these optional parameters that allow to print additional columns in the output file. The first parameter ``additional_columns`` is an :math:`N\times N_\mathrm{tp}`-dimensional array of ``float`` values, where the first dimension refers to the number of additional columns and the second dimension has to be equal to the dimensions of the temperature and pressure arrays.

  The second optional parameter ``additional_columns_desc`` contains an array of ``str`` values with the header descriptions of the additional columns. The dimension has to be equal to number of additional columns. If this is not the case, or if the parameter is missing entirely, the columns will be labelled ``unk`` instead.

The monitor output file has the same format as the one produced by the
``C++`` version discussed :ref:`here<sec:fc_cpp_output>`. Saving
the chemistry output with the panda ``DataFrame`` format in a pickle
file is possible via the function:

.. code:: python

  saveMonitorOutputPandas(file_path, 
                          temperature, 
                          pressure, 
                          element_conserved, 
                          fastchem_flags, 
                          nb_iterations, 
                          nb_chemistry_iterations, 
                          nb_condensation_iterations, 
                          total_element_density, 
                          mean_molecular_weight, 
                          fastchem, 
                          additional_columns=None, 
                          additional_columns_desc=None)

All arguments are identical to those of the previous ASCII output
function. The saved panda ``DataFrame`` contains the same columns and
headers as the simple text output. The only difference between the
outputs is that for the ``DataFrame`` format, the element conservation
and ``FastChem`` flags are not converted to strings (i.e. to ``fail`` or
``ok``) but rather have their original integer values that are returned
by ``FastChem``. Their values are discussed :ref:`here<sec:pfc_input_ouput_struct>`
& :ref:`here<sec:pfc_constants>`. 
