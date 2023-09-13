
.. _sec:fastchem_input_files:

Standard input and output files
===============================

Besides optional other parameter files that are used within the ``C++``
stand-alone version or the Python version, ``FastChem`` requires two
special input files, one for the element abundances and a second
describing the mass action law constant parametrisations. Both are
described in the following.

Element abundance file
~~~~~~~~~~~~~~~~~~~~~~

This file should contain the element abundances for all chemical
elements that are used in ``FastChem``. The location of this file is
usually supplied either within a separate parameter file or directly in
the constructor of the ``FastChem`` object class.

A note on element abundances
----------------------------
                            

It is important to note that there are two different ways to define an
element abundance. Both variants, denoted by :math:`x_j` and
:math:`\epsilon_j`, are related via:

.. math:: \epsilon_j = 10^{x_j - 12}

or

.. math:: x_j = \log\left(\epsilon_j\right) + 12 \ .

In the :math:`x_i` version, widely used in the astronomical literature,
hydrogen has a value of 12 for solar element abundances, such that its
:math:`\epsilon_i` is unity.

| In its input file, ``FastChem`` uses the :math:`x_j` notation, also
  employed in the usual standard abundance compilations (e.g. `Asplund et al. 2009 <http://adsabs.harvard.edu/abs/2009ARA&A..47..481A>`_). 
  For example, in the
  :math:`x_j` notation, the solar element abundance for oxygen is
  :math:`x_\mathrm{O}` = 8.69, whereas its value for
  :math:`\epsilon_\mathrm{O}` would be :math:`0.00048978`.
| Internally, ``FastChem`` converts the :math:`x_j` from the input file
  to the computationally more appropriate :math:`\epsilon_j`. This also
  refers to all methods of the ``FastChem`` object class that are used
  to interact with the element abundances: these will **always** refer
  to :math:`\epsilon_j`.

Thus, if one wants to change the oxygen element abundance in the input
file (which refers to :math:`x_j`) to ten times its solar value, one
would need to use a value of :math:`x_\mathrm{O}` = 1 + 8.69 = 9.69. If
one, on the other hand, uses one of the internal ``FastChem`` methods to
change element abundances on the fly, one would need to set it to a
value of :math:`\epsilon_\mathrm{O} = 10 \cdot 0.00048978 =  0.0048978`.

File structure
--------------
              

The element abundance file should have the following structure to be
readable by ``FastChem``:

.. code:: 

   #Solar element abundances based on Asplund et al. (2009)
   e-  0.00
   Al  6.45
   Ar  6.40
   C   8.43
   Ca  6.34
   Cl  5.50
   Co  4.99
   Cr  5.64
   Cu  4.19
   F   4.56
   Fe  7.50
   Ge  3.65
   H   12.00          
   He  10.93
   K   5.03
   Mg  7.60
   Mn  5.43
   N   7.83
   Na  6.24
   Ne  7.93
   Ni  6.22
   O   8.69
   P   5.41
   S   7.12
   Si  7.51
   Ti  4.95
   V   3.93
   Zn  4.56

The first line is always a header line that provides important
information for the user and is ignored by ``FastChem``. All subsequent
lines contain each the symbol for an element and its element abundance.
Molecules that contain elements not present in this file are ignored.
The element abundance for the electron has an arbitrary value. It is
only present in the file to inform ``FastChem`` that the electrons (and
thus ions) should be included in the chemistry calculations. Its element
abundance :math:`\epsilon_e` will internally be set to 0 because its
number density is determined by charge balance. The elements are not
required to be in any particular order.

Standard files
--------------
              

Together with ``FastChem``, we provide several different element
abundance files, located in the ``input/element_abundances/`` folder.
The folder includes three different sets of element abundances:
`Asplund et al. (2009) <http://adsabs.harvard.edu/abs/2009ARA&A..47..481A>`_,
`Asplund et al. (2021) <http://adsabs.harvard.edu/abs/2021A&A...653A.141A>`_, and
`Lodders (2003) <http://adsabs.harvard.edu/abs/2003ApJ...591.1220L>`_ .

The standard files provide the solar element abundances for species that
are at least as abundant as germanium. As alternative versions, we also
include additional files ``*_extended.dat`` for each compilation, that
includes more elements, up to uranium. These files can be used for the
extended set of ion species described in the next section.

Species data files
~~~~~~~~~~~~~~~~~~

| Another important input is the thermochemical data for all gas phase
  species (molecules and ions) as well as condensates. This includes in
  particular their stoichiometric information as well as a
  parametrisation for their mass action constants. As described in the
  first ``FastChem`` publication
  `Stock et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.479..865S>`_, we use the natural
  logarithm of the dimensionless mass action constant of species
  :math:`i`

  .. math::

     \ln\bar{K}_i(T) = - \frac{\Delta_\mathrm{r} G_i(T)}{R\,T} \ ,
     \label{eq:lnK}

| where :math:`G_i(T)` is the Gibbs free energy of dissociation.
  For ``FastChem``, these mass action constants are fitted with the
  expression

  .. math::

     \ln\bar{K}_i(T) = \frac{a_0}{T} + a_1\,\ln T + b_0 + b_1\,T + b_2\,T^2 \ ,
     \label{eq:fit}

| where :math:`a_0`, :math:`a_1`, :math:`b_0`, :math:`b_1`, and
  :math:`b_2` are the fit coefficients.
| It is in principle possible to use your own parametrisation. For that,
  you need to edit the source code that performs the calculation of the
  mass action constants. For the gas phase species, the corresponding
  code can be found in the source file
  ``gas_phase/molecule_struct.cpp``, while the one for condensates is
  located in ``condensed_phase/condensate_struct.cpp``.

File structure for gas phase species
------------------------------------
                                    

For ``FastChem``, the species information file should have the following
structure:

.. code:: 

   #logK = a1/T + a2 ln T + a3 + a4 T + a5 T^2 for FastChem:
   #includes elements with eps >= eps_Ge
   #fit coefficients calculated from indicated data source.
   Al1Cl1 Aluminum_Chloride : Al 1 Cl 1 # JANAF tables
      6.01726e+04  -9.82181e-01  -5.80778e+00   1.65774e-04  -6.11197e-09

   Al1Cl1F1 Aluminum_Chloride_Fluoride : Al 1 Cl 1 F 1 # JANAF tables
      1.22295e+05  -1.60844e+00  -1.43675e+01   3.72486e-04  -1.98493e-08

   Al1Cl1F2 Aluminum_Chloride_Fluoride : Al 1 Cl 1 F 2 # JANAF tables
      1.93126e+05  -1.90100e+00  -3.00531e+01   6.68640e-04  -3.72957e-08

| The first three lines of the file are treated as header lines and
  discarded when reading in the file.
| The data for each species consists of two lines, while different
  species are separated by a blank line. The first line starts with the
  species’ sum formula. In the standard ``FastChem`` files, we use the
  modified Hill notation for the formulas. Isomeric species would in
  principle have the same formula in the Hill notation. For example, the
  two species and would both be referred to as ``C1H1N1``. To
  distinguish the two in the standard set of ``FastChem``, underscores
  are used, such that ``C1H1N1_1`` refers to , while ``C1H1N1_2``
  represents . The use of the Hill notation is not a requirement. In a
  custom version of the species file, a different chemical notation
  could be used.
| The sum formula is followed by an optional name for the species, for
  example ``Aluminum Chloride``. The name is read until the seperator
  ``:`` is encountered. Note that there has to be a white space between
  the last part of the species’ name and the ``:``.

| After the seperator ``:``, ``FastChem`` expects the stoichiometric
  information of the species, i.e. the elements and their stoichiometric
  coefficients. The elements need to be present in the element abundance
  file as well, otherwise the species will be discarded. They don’t need
  to be in any specific order.
| The stoichiometric information is followed by an optional reference
  for the data. If a reference is used, a separator ``#`` is required
  between the stoichiometry and the reference.
| The second line contains the fit coefficients for the mass action
  constants. ``FastChem`` will read in as many coefficients as it can
  find in that line but for its own parametrisation it will only use the first five.

File structure for condensate species
-------------------------------------
                                     

For ``FastChem``, the species information file should have the following
structure:

.. code:: 

   #logK = a1/T + a2 ln T + a3 + a4 T + a5 T^2 for FastChem:
   #includes elements with eps >= eps_Ge
   #fit coefficients calculated from indicated data source.
   Al2O3(s,l) Aluminum Oxide, Corundum : Al 2 O 3 # JANAF tables
      sl
      2327.0 4000.0
      3.68482e+05 -1.04194e+01 -2.79243e+01 1.17904e-02 -1.58786e-06 
      3.71349e+05 5.273685e+00 -1.35345e+02 1.68173e-03 -8.80259e-08 

   Al2SiO5(s) Aluminum Silicate, Kyanite : Al 2 O 5 Si 1 # JANAF tables
      s
      3000.0
      5.91829e+05 -1.46795e+01 -5.53926e+01 1.428257e-02 -1.5322e-06
      
   CO(l) Carbon Monoxide : C 1 O 1 # Goodwin, R. D., 1985.
     l
     132.85
     1.298145e+05 2.919440e+00 -3.27363e+01 -1.57201e-02 -5.51055e-05

| The first three lines of the file are again header lines and discarded
  when reading in the file.
| For condensates, the data set for each species consists of at least
  five lines. The first one is identical to the one for gas phase
  species described above. One difference with respect to the species’
  formulas compared to those for molecules and ions is that we don’t use
  the Hill notation here. To be consistent with general mineralogy, we
  rather use the normal sum formulas here, together with the phase
  state.

The second line contains the phase state information. ``FastChem``
currently recognises three different states: ``s`` (solid), ``l``
(liquid), and ``sl`` (solid and liquid).

Condensates can have multiple data fits attached to them. Phase changes
can result in non-monotonic slopes within the Gibb’s free energy of
formation. To allow for a better fit to the resulting mass action
coefficients, it is therefore usually beneficial to fit each phase
separately. Additionally, thermochemical data for condensates is often
only tabulated over a restricted temperature range, which makes
extrapolation to higher temperatures problematic.

The third line, thus, lists the temperature limits of the different
fitted data sets in the order they appear below. For solids, the limit
is usually the melting point, whereas for liquids this can typically be
the critical point or the last tabulated data point.

``FastChem`` will use these temperature limits of the phase states when
it calculates the activities of the condensates. By default,
``FastChem`` will also *not* extrapolate the date fits beyond the last
listed temperature. The user can override this behaviour by setting a
special ``FastChem`` parameter. More information on this can be found :ref:`here<sect:fc_param>`.

The following lines finally contain the fit coefficients for the
different data fits. Their format is identical to the corresponding one
for the gas phase species discussed above. The number of data fits has
to be equal to the number of temperature limits discussed above.

.. _standard-files-1:

Standard files
--------------
              

Together with ``FastChem``, we provide two different files for gas phase
species, located in the ``input/logK/`` folder. The file ``logK.dat``
provides the standard set, discussed in
`Stock et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.479..865S>`_ and
`Stock, Kitzmann & Patzer (2022) <http://adsabs.harvard.edu/abs/2022MNRAS.517.4070S>`_. This includes species for
all elements at least as abundant as germanium.

As an alternative version, we also provide an additional file
``logK_extended.dat`` that includes more ions for elements up to
uranium. The data for this file is discussed in
`Hoeijmakers et al. (2019) <http://adsabs.harvard.edu/abs/2019A%26A...627A.165H>`_.

For condensate species we provide the file ``logK_condensates.dat`` in
the ``input/logK/`` folder that contains all species discussed in
`Kitzmann, Stock & Patzer (2023) <http://adsabs.harvard.edu/abs/2023arXiv230902337K>`_. We note that these condensates
should only be used in combination with the standard ``logK.dat`` file
for the gas phase since it lacks condensate species for the additional
elements contained in the species listed in ``logK_extended.dat``.

Basic element data file (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| In addition to the element abundances, ``FastChem`` also needs to have
  additional basic data for the elements, such as their atomic weight to
  calculate the molecular weights of molecules, for example. For most
  elements up to uranium, this data is hard-coded in a standard set
  located in the header file ``chemical_element_data.h``. If you want to
  change this standard set by removing or adding elements or add
  isotopes, you can change it directly in the header file and re-compile
  ``FastChem``.
| Alternatively, ``FastChem`` also has the option to read an external
  file with the required information.

.. _file-structure-1:

File structure
--------------
              

The optional file has the following, simple structure, starting with a
header line that is ignored when reading in the file:

.. code::

   #Basic element data based on Meija et al. (2016)
   e-	Electron	5.4857990907e-4
   H	Hydrogen	1.008
   He	Helium		4.002602
   Li	Lithium		6.94
   Be	Beryllium	9.0121831
   B	Boron		10.81
   C	Carbon		12.011
   N	Nitrogen	14.007
   O	Oxygen		15.999
   F	Fluorine	18.998403163
   Ne	Neon		20.1797
   Na	Sodium		22.98976928
   Mg	Magnesium	24.305
   Al	Aluminium	26.9815385
   Si	Silicon		28.085
   P	Phosphorus	30.973761998
   S	Sulfur		32.06
   Cl	Chlorine	35.45
   Ar	Argon		39.948
   K	Potassium	39.0983
   Ca	Calcium		40.078
   Sc	Scandium	44.955908
   Ti	Titanium	47.867
   Mn	Manganese	54.938044
   Fe	Iron		55.845
   Co	Cobalt		58.933194
   Ni	Nickel		58.6934
   Cu	Copper		63.546
   Zn	Zinc		65.38
   Ga	Gallium		69.723
   Ge	Germanium	72.630
   As	Arsenic		74.921595
   Se	Selenium	78.971
   Br	Bromine		79.904

It contains three columns, where the first one lists the elements’
symbols, the second their names, and the third their atomic weights. An
example of this file can be found in the folder
``fastchem_src/chem_input/``.

.. _sec:fc_parameter_file:

FastChem parameter file (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``FastChem`` is able to load a specific parameter file when one of its
instances is created through the object class constructor. This
parameter file includes the most important parameters and options used
within ``FastChem``. All of these quantities can also be changed during
runtime by using the appropriate methods listed :ref:`here<sec:fastchem_methods>`
for the ``C++`` object class and :ref:`here<sec:pfc_methods>`
for the Python module. Using the parameter
file in principle allows changing these options and parameters outside
of the code and, therefore, does not require the code to be recompiled.

The ``input/`` folder of the ``FastChem`` repository contains two example files:
``parameters.dat`` and ``parameters_py.dat``. The former is to be used for the C++ 
version, while the latter is intended for the Python version. The only difference 
between the two files are the relative locations of the other input files since the Python 
scripts are called from within a separate folder.

.. _file-structure-2:

File structure
--------------
              

The optional parameter file has the following structure:

.. code:: 

   #element abundance file   
   input/element_abundances/asplund_2009.dat

   #species data files    
   input/logK/logK.dat input/logK/logK_condensates.dat

   #accuracy of chemistry iteration
   1.0e-4

   #accuracy of element conservation
   1.0e-4

   #max number of chemistry iterations 
   80000  

   #max number internal solver iterations
   20000

   #element data file (optional)
   input/basic_element_data.dat

The first two entries are the locations of the element abundance and gas
phase & condensate species data files. The condensate species data file
is optional and can be left out or also replaced by ``none``. In that
case, ``FastChem`` won’t be able to perform calculations including
condensation, though.

The next parameter determines the convergence criterion of the chemistry
iteration. This value is also used for the internal Newton’s method and
the condensation calculation iteration. The latter ones can be adjusted
within in the code by the methods listed :ref:`here<sec:fastchem_methods>` 
C++ version and :ref:`here<sec:pfc_methods>` for Python. This is
followed by the (relative) accuracy of the element conservation checks.

The maximum number of chemistry iterations is listed next. The value
here is also used for the maximum number of condensation iterations and
for the number of coupled gas and condensed phase iterations.

The next parameter sets the maximum numbers of iterations for the
different internal numerical methods employed within ``FastChem``. This
includes the Newton, Nelder-Mead, and bisection methods. Using the
corresponding functions of the ``FastChem`` object class, this
number can be adjusted for each of these numerical methods individually.
The last parameter is optional and does not need to be present in the
file. It contains the path to the file for an alternative basic element
data file. If this parameter is not present, ``FastChem`` will use the
standard set that is directly located in the ``FastChem`` source code
(see previous section).

.. _sec:fc_cpp_output:

Output files
~~~~~~~~~~~~

The ``C++`` stand-alone version will produce two output files: a
detailed chemistry output and a monitor file with diagnostic
information. The file names of both files can be chosen in the config
file discussed in the previous section.

Gas phase chemistry output
--------------------------
                          

The chemistry output is organised in columns. The first line of the file
is a header that describes the content of each column.

The first and second column contain the pressure in bar and the
temperature in K, respectively. The third column lists the total number
density of all atoms :math:`i`, i.e.
:math:`n_\mathrm{tot} = \sum_i \left( n_i + \sum_j n_j \nu_{ij} + \sum_c n_c \nu_{ic} \right)`,
summed over their atomic number densities, as well as the ones contained
in all other molecules/ions :math:`j` and the fictitious number
densities of the condensate species :math:`n_c`. This is usually only a
diagnostic quantity and rarely used in other applications.

The fourth column is the number density of the gas in units of
cm\ :math:`^{-3}`, derived from the ideal gas law. This is followed by a
column of mean molecular weights of the mixture of species in units of
the unified atomic mass unit. For all practical purposes, this can also
be converted into units of g/mol.

All subsequent columns contain the number densities (in
cm\ :math:`^{-3}`) or the mixing ratios of all species, depending on the
choice of output made in the config file. By default, elements will be
placed in the beginning, followed by molecules and ions. Note that in
its species data files, ``FastChem`` employs the modified Hill notation
as used in the JANAF thermochemical tables (`Chase, 1986 <https://janaf.nist.gov/>`_)
for the formulas of all
non-element species. If, for example, you are looking for the abundance
of carbon dioxide, you need to locate the ``C1O2`` column rather than
``CO2``, whereas ``NH3`` would be listed as ``H3N1``.

Condensate species output
-------------------------
                         

The output for the condensate species has a similar structure than the
one for the gas phase discussed above. The first two columns again refer
to the pressures and temperatures. The subsequent columns for the
various elements contain the corresponding degrees of condensation.

The remaining columns for the condensate species are the fictitious
number densities in units of cm\ :math:`^{-3}`. As discussed in
`Kitzmann, Stock & Patzer (2023) <http://adsabs.harvard.edu/abs/2023arXiv230902337K>`_, for a given temperature and
pressure the total number of stable condensates that can be present is
limited by the number of elements. Thus, most condensates will usually
have a fictitious number density of 0 cm\ :math:`^{-3}`, which indicates that they are
not present.

Monitor file
------------
            

The monitor output file is a **very** important diagnostic output that
provides crucial details on the outcome of the chemistry calculations.
You should further investigate any chemistry calculations that shows
problems in this file. It is, therefore, advisable to check this file
after each calculation to verify that everything went fine. The first
line of the file is a header that describes the content of each column.

The monitor output is organised in columns, where the first column
contains a simple integer that refers to index of the input
temperature-pressure structure. The second column lists the number of
coupled gas-phase chemistry and condensed phase iterations that were
required to solve the system. If the number is zero, then no condensates
were stable and only a pure gas-phase chemistry calculation was
required.

The third columns contains the total number of gas-phase chemistry
iterations, while the fourth column lists the total number of iterations
for the condensed phase.

The next columns contain information on the convergence of the chemistry
and on the status of overall element conservation. If the chemistry did
converge properly ``ok`` will be listed as output, whereas ``fail`` is
used when the chemistry failed to converge in the maximum allowed number
of steps. The same keywords are used for the element conservation
status: ``ok`` if all elements were conserved, ``fail`` if any element
was not conserved.

The next four columns contain basic chemistry output, that is also found
in the chemistry output file: the pressure, temperature, total element
density, gas number density, and mean molecular weight.

All remaining columns list the status of the element conservation for
each element separately. The same keywords as for the overall element
conservation status are used again in these columns. For the electrons,
this status refers to the charge balance rather than element
conservation. 


Benchmark input and output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The input folder contains a selected sample of atmospheric structures of
various objects, from AGB stars to exoplanets.  The
pre-computed chemistry output of these benchmark structures can be found
in the folder ``output_benchmarks``. This chemistry output has been
generated with the standard ``FastChem`` options and the standard solar
element abundance and equilibrium constants files. These benchmarks can
be used to validate that the FastChem installation works correctly. 
