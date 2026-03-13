# FastChem 4

**An ultra-fast equilibrium chemistry code for gas and condensed phases**

*Authors: Daniel Kitzmann, Joachim Stock*

[![PyPI version](https://img.shields.io/pypi/v/pyfastchem)](https://pypi.org/project/pyfastchem/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Overview

FastChem is an equilibrium chemistry code that calculates the chemical composition of the gas and condensed phase for given temperatures and pressures. The calculation of the gas phase is based on a semi-analytic approach, described in detail in [Stock et al. (2018)](#references) and [Stock et al. (2022)](#references). FastChem Cond (version 3) adds condensation to the code, computing chemical compositions using equilibrium condensation or the rainout approximation commonly used in the field of exoplanets and brown dwarfs, as described in [Kitzmann, Stock & Patzer (2024)](#references). The current version FastChem 4 resolved convergence issues in both gas phase and condensate calculations and added a wide range of new elements and chemical species.

The code is written in object-oriented C++. The current version uses log-based element densities as primary variables, which avoids numerical underflow and allows FastChem to converge at temperatures as low as 1 K (without ions). Overall, the model has been successfully tested for temperatures from 1 K to 6000 K and pressures from 1e-13 bar to 1000 bar for solar element abundances.

Besides the actual FastChem model, we provide a C++ stand-alone version in the *model_src* folder that allows to calculate the equilibrium chemistry for a given temperature-pressure structure. In addition, we provide the Python interface PyFastChem.


## Installation

### PyFastChem (Python)

The easiest way to install FastChem is via pip:

```bash
pip install pyfastchem
```

On most platforms, this will download a pre-compiled package, so no compiler is needed.

### Building from source (C++)

```bash
mkdir build && cd build
cmake ..
make -j
```

This produces the `fastchem` executable and the `fastchem_lib` library. To also build the Python wrapper via CMake:

```bash
cmake -DUSE_PYTHON=ON ..
make -j
```


## PyFastChem

FastChem includes the Python interface PyFastChem that allows running the C++ code as a normal Python module. We provide several examples that show how to call FastChem from within a Python script, including iterating over different metallicity values or C/O ratios, and computing chemical compositions with condensation. The Python examples can be found in the *python* directory.


## Notes on Apple Silicon

Previous versions of FastChem relied on long double (quadruple) precision for convergence at low temperatures. Since Apple Silicon (Mx) processors have no hardware support for quadruple-precision arithmetic, this caused convergence issues on these machines. The current version of FastChem has been redesigned to work fully in double precision, resolving this limitation.

### OpenMP on macOS

FastChem automatically runs multiple calculations in parallel using OpenMP. Apple's default compiler (Apple Clang) does not ship with OpenMP support. To enable OpenMP parallelisation on macOS, install `libomp` via Homebrew:

```bash
brew install libomp
```

CMake and the Python build system (`setup.py`) will automatically detect and use it.


## User Guide

FastChem comes with a user guide available at: https://newstrangeworlds.github.io/FastChem/

It describes the installation and usage of FastChem and covers both the C++ stand-alone version and the Python interface. The manual also contains detailed information on the internal interface functions that the FastChem object class and its Python interface provide.


## References

If you use FastChem in your work, please cite the relevant papers:

- **FastChem (original method):**
  Stock, J. W., Kitzmann, D., Patzer, A. B. C., & Sedlmayr, E. 2018, *MNRAS*, 479, 865.
  [DOI: 10.1093/mnras/sty1531](https://doi.org/10.1093/mnras/sty1531)

- **FastChem 2 (improved gas-phase solver):**
  Stock, J. W., Kitzmann, D., & Patzer, A. B. C. 2022, *MNRAS*, 517, 4070.
  [DOI: 10.1093/mnras/stac2623](https://doi.org/10.1093/mnras/stac2623)

- **FastChem Cond (condensation & rainout):**
  Kitzmann, D., Stock, J. W., & Patzer, A. B. C. 2024, *MNRAS*, 527, 7263.
  [DOI: 10.1093/mnras/stad3515](https://doi.org/10.1093/mnras/stad3515)

<details>
<summary>BibTeX entries</summary>

```bibtex
@article{Stock2018,
  author  = {Stock, J. W. and Kitzmann, D. and Patzer, A. B. C. and Sedlmayr, E.},
  title   = {{FastChem: A computer program for efficient complex chemical equilibrium calculations in the neutral/ionized gas phase with applications to stellar and planetary atmospheres}},
  journal = {Monthly Notices of the Royal Astronomical Society},
  year    = {2018},
  volume  = {479},
  pages   = {865--874},
  doi     = {10.1093/mnras/sty1531}
}

@article{Stock2022,
  author  = {Stock, J. W. and Kitzmann, D. and Patzer, A. B. C.},
  title   = {{FastChem 2: an improved computer program to determine the gas-phase chemical equilibrium composition for arbitrary element distributions}},
  journal = {Monthly Notices of the Royal Astronomical Society},
  year    = {2022},
  volume  = {517},
  pages   = {4070--4080},
  doi     = {10.1093/mnras/stac2623}
}

@article{Kitzmann2024,
  author  = {Kitzmann, D. and Stock, J. W. and Patzer, A. B. C.},
  title   = {{FastChem Cond: equilibrium chemistry with condensation and rainout for cool planetary and stellar environments}},
  journal = {Monthly Notices of the Royal Astronomical Society},
  year    = {2024},
  volume  = {527},
  pages   = {7263--7283},
  doi     = {10.1093/mnras/stad3515}
}
```

</details>


## Licence

This project is Copyright (c) Daniel Kitzmann and Joachim Stock.

FastChem is released under the [GNU General Public Licence (GPL) 3.0](https://www.gnu.org/licenses/gpl-3.0.html). It can be freely copied, edited, and re-distributed. If the code is re-distributed it has to be released under at least a GPL 3.0 licence as well. The full licence can be found in the repository *licence.md* file.

The user guide is released under the Creative Commons Licence (CC BY SA). Licensees may copy and distribute the work and make derivative works based on it only if they give the authors (Daniel Kitzmann & Joachim Stock) the credits by providing a reference to the original guide and this GitHub repository. Licensees may also distribute derivative works only under a license identical to ("not more restrictive than") the license that governs the original work.
