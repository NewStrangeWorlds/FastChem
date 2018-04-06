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


//Standard Newton's method so solve for element densities
//See Sect. 2.4.2
template <class double_type>
void FastChem<double_type>::newtSol(Element<double_type>& species, const double_type h_density, const double_type number_density_min, const unsigned int grid_index)
{

  double_type scaling_factor = 0.0;

  //in case we use the scaling factor, see Appendix A for details
  if (use_scaling_factor)
    scaling_factor = solverScalingFactor(species, number_density_min, h_density, grid_index); // + std::numeric_limits<double_type>::min_exponent/6.0;


  unsigned int order = species.solver_order;


  //for performance reasons, we first calculate the logarithms of the element densities once here
  //and use those values later in the calculation of the molecule densities
  std::vector<double_type> log_element_densities(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    log_element_densities[i] = std::log(elements[i].number_density[grid_index]);



  std::vector<double_type> Aj(order+1, 0.0);

  Aj[0] = std::exp(-scaling_factor) * (number_density_min - species.abundance*h_density);
  Aj[1] = std::exp(-scaling_factor);


  unsigned int index = species.index;


  for (size_t k=1; k<order+1; ++k)
  {

    for (size_t j=0; j<species.molecule_list.size(); ++j)
    {
      unsigned int i = species.molecule_list[j];


      if ( molecules[i].stoichometric_vector[index] == int(k) && molecules[i].abundance == species.abundance)
      {
        molecules[i].sum[grid_index] = 0;

        for (size_t l=0; l<molecules[i].element_indices.size(); ++l)
        {
          unsigned int n = molecules[i].element_indices[l];

          if (n != index)
            molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[n] * log_element_densities[n];
        }


        Aj[k] += std::exp(molecules[i].mass_action_constant[grid_index] + molecules[i].sum[grid_index] - scaling_factor);
      }
    }

    Aj[k] *= k;
  }




  //Newton's method
  bool converged = false;

  double_type x = species.abundance * h_density; //initial guess, ensures monotonous convergence.



  //one Newton step as lambda function
  auto newton_step = [&] (const double_type &x)
    {
      double_type P_j = Aj[order];        //Horner scheme
      double_type P_j_prime = order*Aj[order];

      for (int k = order-1; k >= 1; --k)
      {
        P_j = Aj[k] + x * P_j;
        P_j_prime = k * Aj[k] + x * P_j_prime;
      }

      P_j = Aj[0] + x * P_j;


      return x - P_j/P_j_prime; //Newton step
    };



  //Newton iteration
  unsigned int mu = 0;
  for (mu=0; mu<nb_max_newton_iter; ++mu)
  {
    double_type x_new = newton_step(x);


    if (std::fabs(x_new - x) < newton_err * std::fabs(x_new))  //root found?
    {
      x = x_new;
      converged = true;

      break;
    }


    //prevent x to become negative due to numerical underflow
    if (x_new < 1.e-8*x)
    {
      x_new = 1.e-8*x;
    }


    x = x_new;


    //if (std::isnan(x) || x < 0) break;
  }



  //test if root is in (max(0,x*(1-newton_err)),x*(1+newton_err))
  double_type x_lower = std::fmax(0., x * (1. - newton_err));
  double_type x_upper = x * (1. + newton_err);

  double_type P_j_lower = Aj[order];
  double_type P_j_upper = Aj[order];

  for(int k = order-1; k >=0 ; k--)
  {
    P_j_lower = Aj[k] + x_lower * P_j_lower;
    P_j_upper = Aj[k] + x_upper * P_j_upper;
  }


  species.number_density[grid_index] = x;


  //in case the normal Newton solver does not work, we switch to the alternative version
  if (x < 0 || !converged || P_j_lower*P_j_upper > 0.)
  {
    newtonSolveAlt(species, h_density, grid_index);


    if (verbose_level >= 3)
      std::cout << "FastChem: WARNING: NewtSol failed for species " << species.symbol << " switched to Backup " << x << "\t" << species.number_density[grid_index] << "\n";
  }



  checkN(species, h_density, grid_index);
}






//Alternative version of Newton's method that does not use the n_j_min values
//Takes the law of mass action fully into account; to be used in cases where n_j_min > element_abundance * h_density
//See Sect. 2.4.2
template <class double_type>
void FastChem<double_type>::newtonSolveAlt(Element<double_type>& species, const double_type h_density, const unsigned int grid_index)
{

  double_type scaling_factor = 0.0;

  //in case we use the scaling factor, see Appendix A for details
  if (use_scaling_factor)
    scaling_factor = solverScalingFactor(species, 0.0, h_density, grid_index);


  unsigned int order = 0;

  for (unsigned int i=0; i<species.molecule_list.size(); ++i)
    if (molecules[species.molecule_list[i]].stoichometric_vector[species.index] > int(order) )
      order = molecules[species.molecule_list[i]].stoichometric_vector[species.index];



  //for performance reasons, we first calculate the logarithms of the element densities once here
  //and use those values later in the calculation of the molecule densities
  std::vector<double_type> log_element_densities(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    log_element_densities[i] = std::log(elements[i].number_density[grid_index]);



  std::vector<double_type> Aj(order+1, 0.0);

  Aj[0] = - std::exp(-scaling_factor) * species.abundance*h_density;
  Aj[1] = std::exp(-scaling_factor);


  unsigned int index = species.index;


  for (size_t k=1; k<order+1; ++k)
  {

    for (size_t j=0; j<species.molecule_list.size(); ++j)
    {
      unsigned int i = species.molecule_list[j];


      if ( molecules[i].stoichometric_vector[index] == int(k))
      {
        molecules[i].sum[grid_index] = 0;

        for (size_t l=0; l<molecules[i].element_indices.size(); ++l)
        {
          unsigned int n = molecules[i].element_indices[l];

          if (n != index)
            molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[n] * log_element_densities[n];
        }


        Aj[k] += std::exp(molecules[i].mass_action_constant[grid_index] + molecules[i].sum[grid_index] - scaling_factor);
      }
    }

    Aj[k] *= k;
  }



  //Newton's method
  bool converged = false;

  double_type x = species.abundance * h_density; //Initial guess ensures monotonous convergence.


  //one Newton step as lambda function
  auto newton_step = [&] (const double_type &x)
    {
      double_type P_j = Aj[order];        //Horner scheme
      double_type P_j_prime = order*Aj[order];

      for (int k = order-1; k >= 1; --k)
      {
        P_j = Aj[k] + x * P_j;
        P_j_prime = k * Aj[k] + x * P_j_prime;
      }

      P_j = Aj[0] + x * P_j;


      return x - P_j/P_j_prime; //Newton step
    };


  //Newton iteration
  for (unsigned int mu=0; mu<nb_max_newton_iter; ++mu)
  {
    double_type x_new = newton_step(x);


    if (std::fabs(x_new - x) <= newton_err * std::fabs(x_new))  //root found?
    {
      x = x_new;
      converged = true;

      break;
    }


    //prevent x to become negative due to numerical underflow
    if (x_new < 1.e-8*x) x_new = 1.e-8*x;


    x = x_new;


    if (std::isnan(x)) break;
  }



  // Test if root is in (max(0,x*(1-newton_err)),x*(1+newton_err))
  double_type x_lower = std::fmax(0., x * (1. - newton_err));
  double_type x_upper = x * (1. + newton_err);

  double_type P_j_lower = Aj[order];
  double_type P_j_upper = Aj[order];

  for(int k = order-1; k >=0 ; k--)
  {
    P_j_lower = Aj[k] + x_lower * P_j_lower;
    P_j_upper = Aj[k] + x_upper * P_j_upper;
  }


  species.number_density[grid_index] = x;


  //in case something went wrong again, we try to use another backup
  if (x < 0 || !converged || P_j_lower*P_j_upper > 0.)
  {
    bisectionSolve(species, h_density, grid_index);


    if (verbose_level >= 3)
      std::cout << "FastChem: WARNING: NewtSol Alt failed for species " << species.symbol << " switched to 2nd Backup " << x << "\t" << species.number_density[grid_index] << "\n";
  }



  checkN(species, h_density, grid_index);
}






template class FastChem<double>;
template class FastChem<long double>;

}



