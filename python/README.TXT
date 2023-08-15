This folder contains some example Python scripts that showcases some use cases for FastChem's Python module pyFastChem.
The following scripts are available:

fastchem.py               A basic example that shows how FastChem is initialised and run.

fastchem_metallicity.py   An example that shows how the element abundances are changed on the fly.

fastchem_c_o.py           Shows how the C/O ratio can be changed on the fly.

fastchem_cond.py          Showcases how a calculation with condensation is set up. Uses a temperature structure of a brown dwarf as example.

fastchem_cond_disk.py     An example for condensation in a protoplanetary disk from Kitzmann et al. (2023). 
                          Showcases how some special internal parameters are set. 
                          Warning: Running this script for the full p-T structure of the disk can take *quite* long, depending on your computer! 
                          On a computer without support for quadruple-precision floating point numbers, this script will fail because FastChem 
                          will be unable to resolve the small numbers.
