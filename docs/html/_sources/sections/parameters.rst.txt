.. _sect:fc_param:

Internal parameters
===================

| This section provides an overview over all the internal ``FastChem``
  parameters that can be adjusted via the aforementioned
  
.. code:: c++

  void FastChem.setParameter(std::string param_name, 
                             param_type param_value)
  
method of the C++ object class or with its Python equivalent 

.. code:: python

  pyfastchem.setParameter(str param_name,
                          param_type param_value)

More detailed explanations on what these parameters do can be found in `Stock et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.479..865S>`_, `Stock, Kitzmann & Patzer (2022) <http://adsabs.harvard.edu/abs/2022MNRAS.517.4070S>`_, and `Kitzmann, Stock & Patzer (2023) <http://adsabs.harvard.edu/abs/2023arXiv230902337K>`_

If the ``setParameter`` method is used within Python, the ``C++``\ ``double`` types listed in the following should be replaced by Python ``float`` values, while ``unsigned int`` are used as ``int``. ``PyBind11`` often converts the Python ``bool`` type to an integer value rather than a ``C++``\ ``bool`` type. Setting a boolean parameter will then result in an error message. In such a case, instead of using a simple ``True`` as parameter value, an explicit conversion has to be done instead, for example via ``np.bool_(True)``. 


Standard parameters
-------------------

| ``accuracyChem``
|  Relative accuracy requirement for chemistry iteration and internal solvers (except Newton’s method).
|  Type: ``double``
|  Default value: 1e-5

| ``accuracyElementConservation``
|  Relative accuracy for checking the element conservation.
|  Type: ``double``
|  Default value: 1e-4

| ``accuracyNewton``
|  Relative accuracy requirement for accuracy requirement for Newton’s method.
|  Type: ``double``
|  Default value: 1e-5

| ``accuracyCond``
|  Relative accuracy requirement for condensation calculations.
|  Type: ``double``
|  Default value: 1e-5


| ``nbIterationsChemCond``
|  Maximum number of combined gas phase chemistry and condensate calculations.
|  Type: ``unsigned int``
|  Default value: 3000

| ``nbIterationsChem``
|  Maximum number of chemistry iterations.
|  Type: ``unsigned int``
|  Default value: 3000

| ``nbIterationsCond``
|  Maximum number of condensate calculation iterations.
|  Type: ``unsigned int``
|  Default value: 3000

| ``nbSwitchToNewton``
|  Number of iterations above which the gas phase calculations switch to a multi-dimensional Newton’s method.
|  Type: ``unsigned int``
|  Default value: 400

| ``nbIterationsNewton``
|  Maximum number of iterations for Newton’s method.
|  Type: ``unsigned int``
|  Default value: 3000

| ``nbIterationsBisection``
|  Maximum number of iterations for the bisection method.
|  Type: ``unsigned int``
|  Default value: 3000

| ``nbIterationsNelderMead``
|  Maximum number of iterations for the Nelder-Mead method.
|  Type: ``unsigned int``
|  Default value: 3000


Advanced parameters
-------------------

| The following parameters can be considered potentially dangerous and
  should only be adjusted by users who know exactly what these
  parameters actually do and which side effects they may have.

| ``useScalingFactor``
|  Use the additional scaling factor for the :math:`\log K` in case they become too large. ``FastChem`` will estimate a corresponding scaling factor. An additional one can be added via the ``additionalScaling`` parameter.
|  Type: ``bool``
|  Default value: ``false``

| ``additionalScaling``
|  Additional scaling factor, added to the one ``FastChem`` estimates internally. Only used if ``useScalingFactor`` is set to ``true``.
|  Type: ``double``
|  Default value: 0

| ``minDensityExponentElement``
|  Exponent for the minimum number density of elements.
|  Type: ``double``
|  Default value: -155 (double) or -512 (long double)

| ``minDensityExponentMolecules``
|  Exponent for the minimum number density of molecules&ions.
|  Type: ``double``
|  Default value: -155 (double) or -512 (long double)

| ``maxLogK``
|  Maximum value for the :math:`\log K` of molecules and ions.
|  Type: ``double``
|  Default value: 1e10

| ``condTau``
|  Baseline :math:`\tau` value for condensates.
|  Type: ``double``
|  Default value: 1e-15

| ``condIterChangeLimit``
|  Maximum change of condensate densities per iteration step in log10 space.
|  Type: ``double``
|  Default value: 5

| ``condSolveFullSystem``
|  Solve the full condensate system of equations without eliminating the :math:`\lambda_c`. Iteration will be slower but potentially more stable.
|  Type: ``bool``
|  Default value: ``false``

| ``condReduceSystemSize``
|  Reduce the size of the condensate system of equations by removing unstable condensates.
|  Type: ``bool``
|  Default value: ``true``

| ``condUseFullPivot``
|  Solve the condensate system of equations using an LU decomposition with full pivoting. This enables the detection of singular Jacobians but is slower than the default partial pivoting. If a singular matrix is detected, the system is solved by a perturbed Hessian matrix or by using a singular value decomposition.
|  Type: ``bool``
|  Default value: ``false``

| ``condUseSVD``
|  Solve the condensate system of equations using a singular value decomposition (SVD) in case the Jacobian is singular. This requires ``condUseFullPivot`` to be set to ``true`` as well to detect a singular matrix.
|  Type: ``bool``
|  Default value: ``false``

| ``useCondDataValidityLimit``
|  Use data temperature limits from the condensate input file and avoid extrapolation outside the validity range of the fitted data. Outside the validity range, the activity of condensates is automatically set to a very small value.
|  Type: ``bool``
|  Default value: ``true`` 
 
