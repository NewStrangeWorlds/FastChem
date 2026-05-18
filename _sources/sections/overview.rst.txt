Overview
~~~~~~~~

``FastChem`` is an open-source computer program that can calculate the gas phase chemical equilibrium composition as well as equilibrium condensation and rainout of general systems for a given temperature, pressure, and element abundances.

``FastChem`` uses a semi-analytical approach to solve the non-linear system of mass action law equations for the gas-phase calculations, which results in a massive increase in computational performance over other approaches like Gibbs minimisation.
The general concept and the original version 1.0 are described by `Stock et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.479..865S>`_ (Paper I). Version 1.0, however, is restricted to systems that are dominated by hydrogen and helium and required an additional iteration to account for the pressure of the system. Version 2.0 can be applied to arbitrary element compositions and is described by `Stock, Kitzmann & Patzer (2022) <http://adsabs.harvard.edu/abs/2022MNRAS.517.4070S>`_ (Paper II). Paper III (`Kitzmann, Stock & Patzer (2024) <https://doi.org/10.1093/mnras/stad3515>`_) introduced stoichiometric condensates for calculations of equilibrium condensation and rainout to FastChem.

The current version ``FastChem 4`` resolves convergence issues in both gas phase and condensate calculations, uses log-based element densities as primary variables to avoid numerical underflow, and adds a wide range of new elements and chemical species. This allows ``FastChem`` to converge at temperatures as low as 1 K (without ions).
