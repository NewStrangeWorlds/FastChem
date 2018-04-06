/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
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


#include "fastchem.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//Analytic solution for linear equation, see Sect. 2.4.2 and Eq. (2.32)
template <class double_type>
void FastChem<double_type>::linSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index)
{

  double_type scaling_factor = 0.0;

  //in case we use the scaling factor, see Appendix A for details
  if (use_scaling_factor)
    scaling_factor = solverScalingFactor(species, number_density_min, h_density, grid_index);



  if (scaling_factor > 700.0 && verbose_level >= 3)
    std::cout << "FastChem: WARNING: Underflow in LinSol for element " << species.symbol << "\n";



  //calculation of coefficient A_j1, see Eq. (2.28)
  //referred to as B_j in the following
  unsigned int index = species.index;

  double_type Bj = std::exp(-scaling_factor);


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[index] == 1 && molecules[i].abundance == species.abundance)
    {
      molecules[i].sum[grid_index] = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];

        if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
          molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density[grid_index]);
      }

      Bj += std::exp(molecules[i].mass_action_constant[grid_index] + molecules[i].sum[grid_index] - scaling_factor);
    }
  }


  //calculation of coefficient A_j0, see Eq. (2.27)
  //referred to as C_j in the following
  double_type Cj = std::exp(-scaling_factor) * (number_density_min - species.abundance * h_density);


  //calculation of n_j, Eq. (2.32)
  species.number_density[grid_index] = -Cj/Bj;
  checkN(species, h_density, grid_index);
}



//Analytic solution for quadratic equation, see Sect. 2.4.2 and Eq. (2.32)
template <class double_type>
void FastChem<double_type>::quadSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index)
{

  double_type scaling_factor = 0.0;


  //in case we use the scaling factor, see Appendix A for details
  if (use_scaling_factor)
    scaling_factor = solverScalingFactor(species,number_density_min, h_density, grid_index);



  if (scaling_factor > 700.0 && verbose_level >= 3)
    std::cout << "FastChem: WARNING: Underflow in QuadSol for element " << species.symbol << "\n";



  //calculation of coefficient A_j2, see Eq. (2.29)
  //referred to as A_j in the following
  unsigned int index = species.index;

  double_type Aj = 0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[index] == 2 && molecules[i].abundance == species.abundance)
    {
      molecules[i].sum[grid_index] = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];

        if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
          molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density[grid_index]);

      }

      Aj += std::exp(molecules[i].mass_action_constant[grid_index] + molecules[i].sum[grid_index] - scaling_factor);
    }
  }


  Aj *= 2.;


  if (Aj<1.e-4900L)
  {
    if (verbose_level >= 3) std::cout << "FastChem: Underflow in QuadSol for species " <<  species.symbol << " : switching to LinSol.\n";

    linSol(species, h_density, number_density_min, grid_index);
  }
  else
  {
    //calculation of coefficient A_j1, see Eq. (2.28)
    //referred to as B_j in the following
    double_type Bj = std::exp(-scaling_factor);


    for (size_t j=0; j<species.molecule_list.size(); ++j)
    {
      unsigned int i = species.molecule_list[j];


      if (molecules[i].stoichometric_vector[index] == 1 && molecules[i].abundance == species.abundance)
      {
        molecules[i].sum[grid_index] = 0;

        for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
        {
          unsigned int l = molecules[i].element_indices[k];

          if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
            molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density[grid_index]);
        }

        Bj += std::exp(molecules[i].mass_action_constant[grid_index] + molecules[i].sum[grid_index] - scaling_factor);
      }
    }


    //calculation of coefficient A_j0, see Eq. (2.27)
    //referred to as C_j in the following
    double_type Cj = std::exp(-scaling_factor) * (number_density_min - species.abundance*h_density);


    //calculation of n_j, Eq. (2.32)
    double_type Qj = -0.5 * (Bj + std::sqrt(Bj*Bj - 4.*Aj*Cj));

    species.number_density[grid_index] = Cj/Qj;
    checkN(species, h_density, grid_index);
  }


}





template class FastChem<double>;
template class FastChem<long double>;

}



