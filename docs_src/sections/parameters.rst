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

| ``electronMinTemperature``
|  Temperature (in K) below which the ionisation equilibrium is not solved. For temperatures below this value, all ions are switched off and the electron density is set to its minimum value, i.e. the gas phase is treated as if no ions were present. At low temperatures, the mass action constants of ions become extreme and the ionisation equilibrium is numerically meaningless, while still costing a large number of iterations. Setting this parameter to 0 disables the cap and always solves for the ionisation equilibrium.
|  Type: ``double``
|  Default value: 1000


Advanced parameters
-------------------

| The following parameters can be considered potentially dangerous and
  should only be adjusted by users who know exactly what these
  parameters actually do and which side effects they may have.

| ``minDensityExponentElement``
|  Exponent for the minimum number density of elements.
|  Type: ``double``
|  Default value: -155

| ``minDensityExponentMolecules``
|  Exponent for the minimum number density of molecules and ions.
|  Type: ``double``
|  Default value: -155

| ``condTau``
|  Baseline :math:`\tau` value for condensates.
|  Type: ``double``
|  Default value: 1e-15

| ``condIterChangeLimit``
|  Maximum change of condensate densities per iteration step in log10 space.
|  Type: ``double``
|  Default value: 5

| ``condReduceSystemSize``
|  Reduce the size of the condensate system of equations by removing unstable condensates.
|  Type: ``bool``
|  Default value: ``true``

| ``condUseSVD``
|  Solve the condensate system of equations using a singular value decomposition (SVD) when the Jacobian is singular.
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
|  Number of combined gas phase / condensate outer-loop iterations after which the solver switches to a joint gas+condensate Newton step that couples the gas-phase element densities and the condensate densities. With the default value, the joint step is never triggered during the outer loop. Independent of this parameter, the joint Newton solver is always applied as a post-processing step when the final gas-phase calculation leaves elements unconserved with respect to ``accuracyElementConservation``, which can happen for stiff, strongly-condensing mixtures.
|  Type: ``unsigned int``
|  Default value: 3000

| ``condTraceDensityThreshold``
|  Relative threshold for the maximum possible number density of a condensate. The maximum possible number density of a condensate is limited by its most depleted constituent element. If this value falls below ``condTraceDensityThreshold`` times the total element density, the limiting element has effectively rained out and the condensate is treated as a ghost condensate that cannot form. Ghost condensates are excluded from the condensate Newton solver and from the convergence checks, and their number density is set to zero. This avoids an ill-conditioned system of equations caused by trace elements that have been nearly completely rained out, which would otherwise lead to NaN values in rainout calculations.
|  Type: ``double``
|  Default value: 1e-20

| ``newtonRegularisation``
|  Base value of the Levenberg–Marquardt regularisation that is added to the diagonal of the Jacobian in the multi-dimensional Newton's method of the gas phase. The regularisation is adapted automatically during the iteration: it is increased when a Newton step fails to reduce the element conservation residual and decreased when it succeeds, but never falls below this base value. This keeps the Newton solver stable when the Jacobian becomes near-singular, for example for strongly non-solar compositions or after rainout of elements. A value of 0 corresponds to a pure Newton's method.
|  Type: ``double``
|  Default value: 1e-4

| ``minConvergenceMixingRatioExponent``
|  Exponent for the minimum mixing ratio (number density of a species divided by the total gas number density) a species needs to have in order to be included in the convergence test of the gas-phase iteration. Species below this mixing ratio carry no mass, but high-order molecules (e.g. C\ :sub:`5`) that are proportional to a power of a trace element's density amplify its numerical noise and can oscillate indefinitely, preventing convergence even though they are physically absent. The default value is far below any abundance that could affect the element conservation and therefore does not change the converged result.
|  Type: ``double``
|  Default value: -150

