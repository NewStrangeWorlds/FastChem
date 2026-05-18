.. _sect:fc_param:

Internal parameters
===================

| This section provides an overview of all the internal ``FastChem``
  parameters that can be adjusted via the aforementioned
  
.. code:: c++

  void FastChem.setParameter(std::string param_name, 
                             param_type param_value)
  
method of the C++ object class or with its Python equivalent 

.. code:: python

  pyfastchem.setParameter(str param_name,
                          param_type param_value)

More detailed explanations of what these parameters do can be found in `Stock et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.479..865S>`_, `Stock, Kitzmann & Patzer (2022) <http://adsabs.harvard.edu/abs/2022MNRAS.517.4070S>`_, and `Kitzmann, Stock & Patzer (2024) <https://doi.org/10.1093/mnras/stad3515>`_

When using the ``setParameter`` method in Python, the ``C++``\ ``double`` types listed below should be replaced by Python ``float`` values, while ``unsigned int`` should be passed as ``int``. Python ``bool`` values (``True``/``False``) are correctly dispatched to the ``bool`` overload without requiring explicit conversion. 


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
|  Relative accuracy requirement for Newton’s method.
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
|  Enable the additional scaling factor for the :math:`\log K` values in case they become too large. ``FastChem`` will estimate a corresponding scaling factor. An additional one can be added via the ``additionalScaling`` parameter.
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
|  Exponent for the minimum number density of molecules and ions.
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
|  Solve the condensate system of equations using a singular value decomposition (SVD) when the Jacobian is singular. This requires ``condUseFullPivot`` to be set to ``true`` as well in order to detect a singular matrix.
|  Type: ``bool``
|  Default value: ``false``

| ``useCondDataValidityLimit``
|  Use data temperature limits from the condensate input file and avoid extrapolation outside the validity range of the fitted data. Outside the validity range, the activity of condensates is automatically set to a very small value.
|  Type: ``bool``
|  Default value: ``true``

| ``condUseLM``
|  Use a Levenberg–Marquardt regularisation for the condensate Newton solver. When the Jacobian is singular or near-singular (e.g. fully depleted trace elements after rainout), the solver switches from a plain LU decomposition to the LM normal equations to obtain a well-defined step.
|  Type: ``bool``
|  Default value: ``true``

| ``nbSwitchToJoint``
|  Maximum number of combined gas phase / condensate outer-loop iterations before the solver switches to a joint gas+condensate Newton step. The joint step is also activated adaptively when stagnation or divergence is detected. Setting this to a large value effectively disables the fixed trigger.
|  Type: ``unsigned int``
|  Default value: 3000

| ``condTraceDensityThreshold``
|  Threshold for the maximum possible number density of a condensate (in cm\ :sup:`-3`). Condensates whose maximum possible number density (computed from element abundances and stoichiometry) falls below this value are excluded from the condensate Newton solver and their number density is set to zero. This avoids Jacobian ill-conditioning caused by trace elements that have been nearly completely rained out.
|  Type: ``double``
|  Default value: 1e-100

