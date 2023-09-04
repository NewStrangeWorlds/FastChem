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


#ifndef _fastchem_constants_h
#define _fastchem_constants_h

#include <string>
#include <vector>

namespace fastchem {


//FastChem constants
constexpr unsigned int FASTCHEM_UNKNOWN_SPECIES = 9999999;

constexpr unsigned int FASTCHEM_SUCCESS = 0;
constexpr unsigned int FASTCHEM_NO_CONVERGENCE = 1;
constexpr unsigned int FASTCHEM_INITIALIZATION_FAILED = 2;
constexpr unsigned int FASTCHEM_IS_BUSY = 3;
constexpr unsigned int FASTCHEM_WRONG_INPUT_VALUES = 4;
constexpr unsigned int FASTCHEM_PHASE_RULE_VIOLATION = 5;

const std::vector<std::string> FASTCHEM_MSG = 
{"convergence ok", "convergence failed", "initialisation failed", "is busy", "wrong input value", "phase rule violation"};


//Physical constants
constexpr double CONST_K = 1.380649e-16;    //Boltzmann's constant in erg K-1
constexpr double CONST_AMU = 1.66055e-24;   //Atomic mass unit in g


}
#endif
