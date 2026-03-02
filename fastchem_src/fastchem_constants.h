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
#include <cmath>
#include <limits>

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


//Log-space utilities for numerical stability
constexpr double LOG_DENSITY_FLOOR = -1e8;


//Clamped exponential: returns 0 for very negative arguments, clamps at max representable value
template <class double_type>
inline double_type safeExp(double_type x)
{
  const double_type max_arg = std::log(std::numeric_limits<double_type>::max()) * 0.99;

  if (x < -max_arg) return 0.0;
  if (x > max_arg) x = max_arg;

  return std::exp(x);
}


//Stable computation of ln(sum_k coeffs[k] * exp(log_args[k]))
//Uses the log-sum-exp shift trick to avoid overflow/underflow
template <class double_type>
inline double_type logSumExp(
  const std::vector<double_type>& log_args,
  const std::vector<double_type>& coeffs)
{
  if (log_args.empty()) return static_cast<double_type>(LOG_DENSITY_FLOOR);

  double_type x_max = static_cast<double_type>(LOG_DENSITY_FLOOR);

  for (size_t k = 0; k < log_args.size(); ++k)
    if (log_args[k] > x_max) x_max = log_args[k];

  if (x_max <= static_cast<double_type>(LOG_DENSITY_FLOOR))
    return static_cast<double_type>(LOG_DENSITY_FLOOR);

  double_type sum = 0.0;

  for (size_t k = 0; k < log_args.size(); ++k)
  {
    double_type shifted = log_args[k] - x_max;

    if (shifted > -700)
      sum += coeffs[k] * std::exp(shifted);
  }

  if (sum <= 0.0) return static_cast<double_type>(LOG_DENSITY_FLOOR);

  return x_max + std::log(sum);
}


//Stable computation of log(exp(a) + exp(b)) without overflow/underflow
//The exp argument is always <= 0, so log1p argument is always in [0, 1]
template <class double_type>
inline double_type logAddExp(double_type a, double_type b)
{
  if (a > b) return a + std::log1p(std::exp(b - a));
  else       return b + std::log1p(std::exp(a - b));
}


}
#endif
