/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2024 Daniel Kitzmann, Joachim Stock
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
#include <iostream>

#include "../_ext/Eigen/Dense"


namespace fastchem {


//Joint Newton step coupling gas-phase element densities and condensate densities.
//Solves for y_c = ln(n_c), y_lambda = ln(lambda_c), y_j = ln(n_j) simultaneously.
//Uses code convention: J * delta = rhs, where J = -dF/dx, rhs = F, update: x += delta.
//Equations match the condensed-phase formulation (Eq. 33 in Paper III):
//  Complementarity: ln(tau) = y_c + y_lambda  (i.e., tau = n_c * lambda)
//  Activity: -log_activity - lambda = 0
//  Element conservation: epsilon_j * N_total - A_j = 0  (with N_total constant in Jacobian)
//  Gas pressure: n_gas - P/(kT) = 0  (replaces one element conservation)
void FastChem::jointNewtonStep(
  const std::vector<Condensate*>& condensates_act_all,
  const double temperature,
  const double gas_density,
  double& total_element_density)
{
  // Filter to only include stable condensates (log_activity > LOG_ACTIVITY_THRESHOLD)
  std::vector<Condensate*> condensates_stable;
  condensates_stable.reserve(condensates_act_all.size());
  for (auto* c : condensates_act_all)
    if (c->log_activity > LOG_ACTIVITY_THRESHOLD)
      condensates_stable.push_back(c);

  const auto& condensates_act = condensates_stable;
  const size_t nc = condensates_act.size();
  if (nc == 0) return;

  // Build list of element indices with epsilon > 0 (excludes e-)
  std::vector<int> elem_idx;
  elem_idx.reserve(element_data.nb_elements);
  for (size_t k = 0; k < element_data.nb_elements; ++k)
    if (element_data.elements[k].epsilon > 0.0)
      elem_idx.push_back(static_cast<int>(k));

  const size_t ne = elem_idx.size();
  if (ne == 0) return;

  // Global element index -> local index map
  const size_t total_nb_elements = element_data.nb_elements;
  std::vector<int> global_to_local(total_nb_elements, -1);
  for (size_t j = 0; j < ne; ++j)
    global_to_local[elem_idx[j]] = static_cast<int>(j);

  // System size: 2*nc + ne
  const size_t N = 2 * nc + ne;

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(N, N);
  Eigen::VectorXd rhs(N);
  rhs.setZero();

  // Precompute condensate sigma_c = sum of stoichiometric coefficients
  std::vector<double> sigma_c(nc);
  for (size_t i = 0; i < nc; ++i)
  {
    double s = 0;
    for (auto idx : condensates_act[i]->element_indices)
      s += condensates_act[i]->stoichiometric_vector[idx];
    sigma_c[i] = s;
  }

  // ---- Block 1: Complementarity rows (0..nc-1) ----
  // F_i = ln(tau) - y_c - y_lambda,  J = -dF/dx
  // J(i,i) = 1 (= -d/dy_c[-1]), J(i,nc+i) = 1 (= -d/dy_lambda[-1])
  for (size_t i = 0; i < nc; ++i)
  {
    double log_n_c = safeLog(condensates_act[i]->number_density);
    double log_lambda = safeLog(condensates_act[i]->activity_correction);

    J(i, i) = 1.0;
    J(i, nc + i) = 1.0;

    rhs(i) = condensates_act[i]->log_tau - log_n_c - log_lambda;
  }

  // ---- Block 2: Activity rows (nc..2nc-1) ----
  // F_i = -log_activity - lambda,  where lambda = exp(y_lambda)
  // log_activity = mass_action_constant + sum_j(nu_cij * y_j)  [PLUS sign in FastChem]
  // dF/dy_lambda = -lambda  →  J = lambda
  // dF/dy_j = -nu_cij  →  J = nu_cij
  for (size_t i = 0; i < nc; ++i)
  {
    const size_t row = nc + i;
    double lambda_i = condensates_act[i]->activity_correction;

    J(row, nc + i) = lambda_i;

    for (auto idx : condensates_act[i]->element_indices)
    {
      int nu = condensates_act[i]->stoichiometric_vector[idx];
      int jl = global_to_local[idx];
      if (jl >= 0)
        J(row, 2 * nc + jl) = nu;
    }

    rhs(row) = -condensates_act[i]->log_activity - lambda_i;
  }

  // ---- Precompute C_k and D_k for element/pressure rows ----
  // C_k = n_k + sum_m sigma_m * nu_mk * n_m   (d(N_gas_elem)/d(y_k))
  // D_k = n_k + sum_m nu_mk * n_m              (d(n_gas)/d(y_k))
  std::vector<double> C_k(ne, 0.0);
  std::vector<double> D_k(ne, 0.0);

  for (size_t ki = 0; ki < ne; ++ki)
  {
    const int k = elem_idx[ki];
    const auto& ek = element_data.elements[k];
    C_k[ki] = ek.number_density;
    D_k[ki] = ek.number_density;

    for (auto m : ek.molecule_list)
    {
      const auto& mol = gas_phase.molecules[m];
      double nu_mk = mol.stoichiometric_vector[k];
      double n_m = mol.number_density;

      D_k[ki] += nu_mk * n_m;
      C_k[ki] += mol.sigma * nu_mk * n_m;
    }
  }

  // Choose j_ref: element with LARGEST D_k (most coupling to gas pressure).
  // Its conservation equation is most redundant with the gas pressure equation,
  // so removing it loses the least information.
  int j_ref_local = 0;
  {
    double max_Dk = D_k[0];
    for (size_t j = 1; j < ne; ++j)
    {
      if (D_k[j] > max_Dk)
      {
        max_Dk = D_k[j];
        j_ref_local = static_cast<int>(j);
      }
    }
  }

  // ---- Block 3: Element conservation rows (2nc .. 2nc+ne-2, skipping j_ref) ----
  // F_j = epsilon_j * N_total - A_j,  with N_total treated as CONSTANT for Jacobian
  // (avoids ill-conditioning from near-cancellation of eps*C and B terms)
  // J = -dF/dx (with N_total constant):
  //   condensate col i: J = nu_cij * n_ci
  //   element col k:    J = B_jk  (= delta_jk * n_j + sum_m nu_mj * nu_mk * n_m)
  size_t row_offset = 2 * nc;
  for (size_t ji = 0; ji < ne; ++ji)
  {
    if (static_cast<int>(ji) == j_ref_local) continue;

    const int j = elem_idx[ji];
    const auto& ej = element_data.elements[j];
    const size_t row = row_offset++;

    // Condensate columns: nu_cij * n_ci
    for (size_t i = 0; i < nc; ++i)
    {
      double nu_cij = condensates_act[i]->stoichiometric_vector[j];
      double n_ci = condensates_act[i]->number_density;
      J(row, i) = nu_cij * n_ci;
    }

    // Element columns: B_jk
    // Diagonal: n_j
    J(row, 2 * nc + ji) = ej.number_density;

    // Molecule coupling: sum_m nu_mj * nu_mk * n_m
    for (auto m : ej.molecule_list)
    {
      const auto& mol = gas_phase.molecules[m];
      double nu_mj = mol.stoichiometric_vector[j];
      double n_m = mol.number_density;
      if (n_m == 0.0) continue;

      for (auto k_global : mol.element_indices)
      {
        int kl = global_to_local[k_global];
        if (kl >= 0)
          J(row, 2 * nc + kl) += nu_mj * mol.stoichiometric_vector[k_global] * n_m;
      }
    }

    // RHS = F_j = epsilon_j * N_total - A_j
    // A_j sums over ALL condensates (stable and unstable)
    double Aj = ej.number_density;
    for (auto m : ej.molecule_list)
      Aj += gas_phase.molecules[m].stoichiometric_vector[j] * gas_phase.molecules[m].number_density;
    for (auto c : ej.condensate_list)
      Aj += condensed_phase.condensates[c].stoichiometric_vector[j] * condensed_phase.condensates[c].number_density;

    rhs(row) = ej.epsilon * total_element_density - Aj;
  }

  // ---- Block 4: Gas pressure row (last row) ----
  // F = n_gas - gas_density,  J = -dF/dx
  // dF/dy_k = D_k,  so J(row, 2nc+k) = -D_k
  {
    const size_t row = N - 1;

    for (size_t ki = 0; ki < ne; ++ki)
      J(row, 2 * nc + ki) = -D_k[ki];

    double n_gas = 0;
    for (size_t ki = 0; ki < ne; ++ki)
      n_gas += element_data.elements[elem_idx[ki]].number_density;
    for (size_t m = 0; m < gas_phase.molecules.size(); ++m)
      n_gas += gas_phase.molecules[m].number_density;
    if (element_data.e_ != FASTCHEM_UNKNOWN_SPECIES)
      n_gas += element_data.elements[element_data.e_].number_density;

    rhs(row) = n_gas - gas_density;
  }

  // ---- Diagnostics (verbose >= 4) ----
  if (options.verbose_level >= 4)
  {
    std::cout << "Joint Newton: N=" << N << " (nc=" << nc << ", ne=" << ne << ")\n";
    std::cout << "  j_ref: " << element_data.elements[elem_idx[j_ref_local]].symbol
              << " (D_k=" << D_k[j_ref_local] << ")\n";
    double rhs_norm = rhs.norm();
    std::cout << "  |rhs| = " << rhs_norm << "\n";
    int max_idx;
    rhs.cwiseAbs().maxCoeff(&max_idx);
    std::cout << "  max rhs[" << max_idx << "] = " << rhs(max_idx);
    if (max_idx < (int)nc) std::cout << " (compl " << condensates_act[max_idx]->symbol << ")";
    else if (max_idx < (int)(2*nc)) std::cout << " (activ " << condensates_act[max_idx-nc]->symbol << ")";
    else std::cout << " (elem/pres)";
    std::cout << "\n";
  }

  // ---- Row scaling ----
  Eigen::VectorXd scaling = J.cwiseAbs().rowwise().maxCoeff();
  for (int i = 0; i < scaling.rows(); ++i)
    if (scaling(i) == 0.0) scaling(i) = 1.0;
  for (int col = 0; col < J.cols(); ++col)
    J.col(col).array() /= scaling.array();
  rhs.array() /= scaling.array();

  // ---- Solve ----
  Eigen::VectorXd delta;
  {
    Eigen::PartialPivLU<Eigen::MatrixXd> solver;
    solver.compute(J);
    delta = solver.solve(rhs);

    if (!delta.allFinite())
      delta = J.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
  }

  // ---- Damp large steps ----
  const double max_step = 2.0;
  double max_delta = delta.cwiseAbs().maxCoeff();

  if (options.verbose_level >= 4)
  {
    int max_idx;
    delta.cwiseAbs().maxCoeff(&max_idx);
    std::cout << "  max |delta| = " << max_delta << " at idx " << max_idx;
    if (max_idx < (int)nc) std::cout << " (nc:" << condensates_act[max_idx]->symbol << ")";
    else if (max_idx < (int)(2*nc)) std::cout << " (lam:" << condensates_act[max_idx-nc]->symbol << ")";
    else std::cout << " (elem:" << element_data.elements[elem_idx[max_idx-2*nc]].symbol << ")";
    std::cout << "\n";
    for (size_t ji = 0; ji < ne; ++ji)
      std::cout << "  d_n[" << element_data.elements[elem_idx[ji]].symbol << "]=" << delta(2*nc+ji) << "\n";
  }

  if (max_delta > max_step)
    delta *= max_step / max_delta;

  // ---- Apply updates ----
  // Condensate densities (y_c += delta_c)
  for (size_t i = 0; i < nc; ++i)
  {
    double log_n_c = safeLog(condensates_act[i]->number_density);

    log_n_c += delta(i);
    condensates_act[i]->number_density = safeExp(log_n_c);
  }

  // Activity corrections (y_lambda += delta_lambda)
  for (size_t i = 0; i < nc; ++i)
  {
    double log_lambda = safeLog(condensates_act[i]->activity_correction);

    log_lambda += delta(nc + i);
    condensates_act[i]->activity_correction = safeExp(log_lambda);
  }

  // Element densities (y_j += delta_j)
  for (size_t ji = 0; ji < ne; ++ji)
  {
    const int j = elem_idx[ji];
    element_data.elements[j].log_number_density += delta(2 * nc + ji);
    element_data.elements[j].number_density =
      safeExp(element_data.elements[j].log_number_density);
  }

  // ---- Refresh derived quantities ----
  gas_phase.updateMoleculeDensities();

  total_element_density =
    gas_phase.totalElementDensity() + condensed_phase.totalElementDensity();

  updatePhi(total_element_density);
}


}
