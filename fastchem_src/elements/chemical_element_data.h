
/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2021 Daniel Kitzmann, Joachim Stock
*
* FastChem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* FastChem is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* FastChem directory under <license.md>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _chemical_element_data_h
#define _chemical_element_data_h

#include "../species_struct.h"
#include <vector>


namespace fastchem {


const std::vector<ChemicalElementData> standard_chemical_element_data{ 
  {"e-", "Electron", 5.4857990907e-4, 0.0},
  {"H",  "Hydrogen", 1.008,  0.0},
  {"He", "Helium", 4.002602, 0.0},
  {"Li", "Lithium", 6.94, 0.0},
  {"Be", "Beryllium", 9.0121831, 0.0},
  {"B",  "Boron", 10.81, 0.0},
  {"C",  "Carbon", 12.011, 0.0},
  {"N",  "Nitrogen", 14.007, 0.0},
  {"O",  "Oxygen", 15.999, 0.0},
  {"F",  "Fluorine", 18.998403163, 0.0},
  {"Ne", "Neon", 20.1797, 0.0},
  {"Na", "Sodium", 22.98976928, 0.0},
  {"Mg", "Magnesium", 24.305, 0.0},
  {"Al", "Aluminium", 26.9815385, 0.0},
  {"Si", "Silicon", 28.085, 0.0},
  {"P",  "Phosphorus", 30.973761998, 0.0},
  {"S",  "Sulfur", 32.06, 0.0},
  {"Cl", "Chlorine", 35.45, 0.0},
  {"Ar", "Argon", 39.948, 0.0},
  {"K",  "Potassium", 39.0983, 0.0},
  {"Ca", "Calcium", 40.078, 0.0},
  {"Sc", "Scandium", 44.955908, 0.0},
  {"Ti", "Titanium", 47.867, 0.0},
  {"V",  "Vanadium", 50.9415, 0.0},
  {"Cr", "Chromium", 51.9961, 0.0},
  {"Mn", "Manganese", 54.938044, 0.0},
  {"Fe", "Iron", 55.845, 0.0},
  {"Co", "Cobalt", 58.933194, 0.0},
  {"Ni", "Nickel", 58.6934, 0.0},
  {"Cu", "Copper", 63.546, 0.0},
  {"Zn", "Zinc", 65.38, 0.0},
  {"Ga", "Gallium", 69.723, 0.0},
  {"Ge", "Germanium", 72.630, 0.0},
  {"As", "Arsenic", 74.921595, 0.0},
  {"Se", "Selenium", 78.971, 0.0},
  {"Br", "Bromine", 79.904, 0.0},
  {"Kr", "Krypton", 83.8, 0.0},
  {"Rb", "Rubidium", 85.4678, 0.0},
  {"Sr", "Strontium", 87.62, 0.0},
  {"Y",  "Yttrium", 88.9059, 0.0},
  {"Zr", "Zirconium", 91.224, 0.0},
  {"Nb", "Niobium", 92.9064, 0.0},
  {"Mo", "Molybdenum", 95.94, 0.0},
  {"Tc", "Technetium", 98, 0.0},
  {"Ru", "Ruthenium", 101.07, 0.0},
  {"Rh", "Rhodium", 102.9055, 0.0},
  {"Pd", "Palladium", 106.42, 0.0},
  {"Ag", "Silver", 107.8682, 0.0},
  {"Cd", "Cadmium", 112.411, 0.0},
  {"In", "Indium", 114.818, 0.0},
  {"Sn", "Tin", 118.71, 0.0},
  {"Sb", "Antimony", 121.76, 0.0},
  {"I",  "Iodine", 126.9045, 0.0},
  {"Te", "Tellurium", 127.6, 0.0},
  {"Xe", "Xenon", 131.293, 0.0},
  {"Cs", "Caesium", 132.9055, 0.0},
  {"Ba", "Barium", 137.327, 0.0},
  {"La", "Lanthanum", 138.9055, 0.0},
  {"Ce", "Cerium", 140.116, 0.0},
  {"Pr", "Praesodymium", 140.9077, 0.0},
  {"Nd", "Neodymium", 144.24, 0.0},
  {"Pm", "Promethium", 145, 0.0},
  {"Sm", "Samarium", 150.36, 0.0},
  {"Eu", "Europium", 151.964, 0.0},
  {"Gd", "Gadolinium", 157.25, 0.0},
  {"Tb", "Terbium", 158.9253, 0.0},
  {"Dy", "Dysprosium", 162.5, 0.0},
  {"Ho", "Holmium", 164.9303, 0.0},
  {"Er", "Erbium", 167.259, 0.0},
  {"Tm", "Thulium", 168.9342, 0.0},
  {"Yb", "Ytterbium", 173.04, 0.0},
  {"Lu", "Lutetium", 174.967, 0.0},
  {"Hf", "Hafnium", 178.49, 0.0},
  {"Ta", "Tantalum", 180.9479, 0.0},
  {"W",  "Tungsten", 183.84, 0.0},
  {"Re", "Rhenium", 186.207, 0.0},
  {"Os", "Osmium", 190.23, 0.0},
  {"Ir", "Iridium", 192.217, 0.0},
  {"Pt", "Platinum", 195.078, 0.0},
  {"Au", "Gold", 196.9665, 0.0},
  {"Hg", "Mercury", 200.59, 0.0},
  {"Tl", "Thallium", 204.3833, 0.0},
  {"Pb", "Lead", 207.2, 0.0},
  {"Bi", "Bismuth", 208.9804, 0.0},
  {"Th", "Thorium", 232.0381, 0.0},
  {"U",  "Uranium", 238.0289, 0.0}
  };

}

#endif