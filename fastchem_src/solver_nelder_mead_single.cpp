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

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>


namespace fastchem {


//Mass action law for a single element, used in the Nelder-Mead method
template <class double_type>
double_type FastChem<double_type>::nelderMeadFunction(Element<double_type>& species, const double_type x, const double_type h_density, const unsigned int grid_index)
{
  double_type f_i = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    unsigned int molecule_index = species.molecule_list[j];


    double_type f_j = 0.0;


    for (size_t l=0; l<molecules[molecule_index].element_indices.size(); ++l)
    {
      unsigned int element_index = molecules[molecule_index].element_indices[l];

      double_type species_density = std::log(elements[element_index].number_density[grid_index]);

      //just in case...
      if (std::isnan(species_density) || std::isinf(species_density)) species_density = 0.0;

      if (element_index != species.index)
        f_j += molecules[molecule_index].stoichometric_vector[element_index] * species_density;
      else
        f_j += molecules[molecule_index].stoichometric_vector[species.index] * x;
    }


    f_j += molecules[molecule_index].mass_action_constant[grid_index];

    f_j = std::exp(f_j) * molecules[molecule_index].stoichometric_vector[species.index];

    f_i += f_j;
  }


  if (species.symbol == "e-")
    f_i = (- f_i - std::exp(x));
  else
    f_i = (species.abundance * h_density - f_i - std::exp(x));


  return f_i;
}



//Nelder-Mead downhill simplex method in one dimension
template <class double_type>
bool FastChem<double_type>::nelderMeadSimplexSolve(Element<double_type>& species, const double_type initial_solution, const double h_density, const unsigned int grid_index)
{
  unsigned int N = 1; //dimension

  //construct initial simplex
  std::vector<double_type> x;

  double_type initial_distance = (1.0 + 0.05) * initial_solution;

  double_type simplex_point = initial_distance;
  x.push_back(simplex_point);

  x.push_back(initial_solution);



  std::vector<double_type> vf(N+1,0);   //values of function evaluated at simplex vertexes

  unsigned int x1 = 0;    //index of best solution
  unsigned int xn = 0;    //index of second worst solution
  unsigned int xnp1 = 0;  //index of worst solution



  double_type rho = 1.0, chi = 2.0, psi = 0.5, sigma = 0.5; //standard coefficients for Nelder-Mead method


  unsigned int nb_iterations = nb_max_neldermead_iter;
  bool converged = false;


  //downhill simplex method starts
  for (unsigned int iter_step = 0; iter_step < nb_iterations; ++iter_step)
  {

    for(size_t i = 0; i < N + 1; ++i)
	    vf[i] = std::fabs(nelderMeadFunction(species, x[i], h_density, grid_index));


    x1 = 0; xn = 0; xnp1 = 0;

    for (size_t i=0; i<vf.size(); ++i)
    {
      if(vf[i] < vf[x1]) x1 = i;

      if(vf[i] > vf[xnp1]) xnp1 = i;
    }

    xn = xnp1;  //in 1D they are equal


    double_type xg = x[x1]; //xg: centroid of the N best vertexes; in 1D corresponds to the solution x1


    //check if the function has a root in a delta region around xg
    double_type delta = xg * accuracy*1e-4;

    double_type vf_epsilon_plus = nelderMeadFunction(species, xg+delta, h_density, grid_index);
    double_type vf_epsilon_minus = nelderMeadFunction(species, xg-delta, h_density, grid_index);


    //if the function changes sign, we have the solution
    if ((vf_epsilon_minus < 0 && vf_epsilon_plus > 0) || (vf_epsilon_minus > 0 && vf_epsilon_plus < 0)  )
    {
      converged = true;

      break;
    }


    //reflection:
    double_type xr = (1.0 + rho) * xg - rho*x[xnp1];

    double_type fxr = std::fabs(nelderMeadFunction(species, xr, h_density, grid_index));


    if (fxr < vf[x1])
    {
      double_type xe;

      //expansion
      xe = (1.0 + rho*chi) * xg - rho*chi*x[xnp1];

      double_type fxe = std::fabs(nelderMeadFunction(species, xe, h_density, grid_index));

      if (fxe < fxr)
        x[xnp1] = xe;
      else
        x[xnp1] = xr;
    }
    else
    {
      if (fxr < vf[xn])
        x[xnp1] = xr;
      else
      {
        bool perform_shrink = false;

        if (fxr < vf[xnp1])
        {
          double_type xc;

          //outward contraction
          xc = (1.0 + psi*rho) * xg - psi*rho*x[xnp1];


          double_type fxc = std::fabs(nelderMeadFunction(species, xc, h_density, grid_index));

          if (fxc <= fxr)
            x[xnp1] = xc;
          else
            perform_shrink = true;

        }
        else
        {
          double_type xc;

          //inward contraction
          xc = (1.0 - psi) * xg + psi*rho*x[xnp1];


          double_type fxc = std::fabs(nelderMeadFunction(species, xc, h_density, grid_index));


          if (fxc < vf[xnp1])
            x[xnp1] = xc;
          else
            perform_shrink = true;
        }

        if (perform_shrink)
        {
          for(size_t i=0; i<x.size(); ++i )
          {

            //total contraction
            if (i != x1)
              x[i] = x[x1] + sigma*(x[i] - x[x1]);

          }

        }

      }

    }


  }//optimisation is finished


  species.number_density[grid_index] = std::exp(x[x1]);


  if (!converged && verbose_level >= 3)
    std::cout << "Nelder-Mead iteration limit reached, result may not be optimal." << "\t" << x[x1] << "\n";


  return converged;
}



template class FastChem<double>;
template class FastChem<long double>;

}
