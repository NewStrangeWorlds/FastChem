Overview
~~~~~~~~

``FastChem`` is an open-source computer program that can calculate the gas phase chemical equilibrium composition as well as equilibrium condensation and rainout of general systems for a given temperature pressure and element abundances. Its latest version 3.x is also referred to as ``FastChem Cond``

``FastChem`` uses a semi-analytical approach to solve the non-linear system of mass action law equations for the gas-phase calculations, which results in a massive increase in computational performance over other approaches like Gibbs minimisation. 
The general concept and the original version 1.0 is described by `Stock et al. (2018) <http://adsabs.harvard.edu/abs/2018MNRAS.479..865S>`_ (Paper I). Version 1.0, however, is restricted to systems that are dominated by hydrogen and helium and required an additional iteration to account for the pressure of the system. The current version 2.0 can now be applied to arbitrary element compositions. This version will is described by `Stock, Kitzmann & Patzer (2022) <http://adsabs.harvard.edu/abs/2022MNRAS.517.4070S>`_, from here on referred to as Paper II. Paper III (`Kitzmann, Stock & Patzer (2023) <http://adsabs.harvard.edu/abs/2023arXiv230902337K>`_) introduced stoichiometric condensates for calculations of equilibrium condensation and rainout to FastChem. 

``FastChem`` has already been applied to numerous different systems,
from brown dwarfs (`Kitzmann et al. 2020 <http://adsabs.harvard.edu/abs/2020ApJ...890..174K>`_), to
mini-Neptunes, hot-Jupiters (`Bourrier et al. 2020 <http://adsabs.harvard.edu/abs/2020A%26A...637A..36B>`_), to ultra-hot Jupiters
(`Hoeijmakers et al. 2019 <http://adsabs.harvard.edu/abs/2019A%26A...627A.165H>`_). It is directly
coupled to the retrieval model ``Helios-r2`` (`Kitzmann et al. 2020 <http://adsabs.harvard.edu/abs/2020ApJ...890..174K>`_), to the general
atmospheric model ``HELIOS`` (`Malik et al. 2019 <http://adsabs.harvard.edu/abs/2019AJ....157..170M>`_), and the non-equilibrium
chemistry ``VULCAN`` (`Tsai et al. 2018 <http://adsabs.harvard.edu/abs/2018ApJ...862...31T>`_), all of which are available under https://github.com/exoclime.

