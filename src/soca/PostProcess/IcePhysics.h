/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

namespace soca {
namespace icephysics {

struct Constants {
  static constexpr double puny       = 1.0e-11;
  static constexpr double ain_min    = 1.0e-8;
  static constexpr double rho_ice    = 917.0;
  static constexpr double rho_snow   = 330.0;
  static constexpr double rho_ocean  = 1025.0;
};

// ---------------------------------------------------------------------------
// adjustThicknessCategories
//
// Re-bin per-category area (aicen) and volume (vicen) so that for every
// category k either aicen[k] == 0 or hicat[k] <= vicen[k]/aicen[k] <= hicat[k+1].
// Conserves Sum(aicen) = aiceTotal (taken as sum of input aicen) and
// Sum(vicen) = viceTarget on output. dhiMin enforces a minimum thickness gap
// from the lower category bound to avoid degenerate bins.
//
// hicat must have size aicen.size()+1 and be monotonically increasing.
// Returns true on success, false if no consistent allocation is possible
// (e.g. viceTarget is incompatible with aiceTotal and the bounds).
// ---------------------------------------------------------------------------
bool adjustThicknessCategories(std::vector<double> & aicen,
                               std::vector<double> & vicen,
                               double viceTarget,
                               const std::vector<double> & hicat,
                               double dhiMin = 0.01);

// ---------------------------------------------------------------------------
// enforceFreeboard
//
// Enforce hydrostatic balance: rho_ice*hi + rho_snow*hs <= rho_ocean*(hi+hs)
// per category. First redistribute snow across categories (move snow off the
// thinnest, most-flooded cats onto thicker cats with headroom). If any
// category remains negative-freeboard, slightly grow its ice volume to lift
// the snow-ice interface back to sea level. Operates on a single column.
//
// All vectors must have the same length (= ncat).
// Returns true if a non-negative freeboard configuration was reached.
// ---------------------------------------------------------------------------
bool enforceFreeboard(std::vector<double> & aicen,
                      std::vector<double> & vicen,
                      std::vector<double> & vsnon,
                      double rhoIce   = Constants::rho_ice,
                      double rhoSnow  = Constants::rho_snow,
                      double rhoOcean = Constants::rho_ocean);

// Per-category freeboard (m). freeboard[k] = hi[k]*(rho_ocean-rho_ice)/rho_ocean
//                                          - hs[k]*rho_snow/rho_ocean
// where hi = vicen/aicen, hs = vsnon/aicen. Zero where aicen <= puny.
std::vector<double> snowIceFreeboard(const std::vector<double> & aicen,
                                     const std::vector<double> & vicen,
                                     const std::vector<double> & vsnon,
                                     double rhoIce   = Constants::rho_ice,
                                     double rhoSnow  = Constants::rho_snow,
                                     double rhoOcean = Constants::rho_ocean);

}  // namespace icephysics
}  // namespace soca
