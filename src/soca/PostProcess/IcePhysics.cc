/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/PostProcess/IcePhysics.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace soca {
namespace icephysics {

namespace {

double sum(const std::vector<double> & v) {
  return std::accumulate(v.begin(), v.end(), 0.0);
}

}  // namespace

// ---------------------------------------------------------------------------
// adjustThicknessCategories
//
// Strategy (port of mc6util.adjust_thkncats_aice):
//   1. Build target per-cat thickness from current aicen/vicen, fall back to
//      mid-bin when a category has aicen but degenerate volume.
//   2. Clamp each per-cat thickness to [hicat[k]+dhiMin, hicat[k+1]-dhiMin].
//   3. Recompute vicen[k] = aicen[k] * h_clamped[k]. The total vice will then
//      not match viceTarget; correct by spreading the residual onto categories
//      that still have headroom (proportional to remaining headroom), iterating
//      a few times. If the residual cannot be absorbed (target outside the
//      feasible envelope sum_k aicen[k]*[h_lo,h_hi]) we return false.
// ---------------------------------------------------------------------------
bool adjustThicknessCategories(std::vector<double> & aicen,
                               std::vector<double> & vicen,
                               double viceTarget,
                               const std::vector<double> & hicat,
                               double dhiMin) {
  const size_t ncat = aicen.size();
  if (ncat == 0 || vicen.size() != ncat || hicat.size() != ncat + 1) {
    return false;
  }

  std::vector<double> hLo(ncat), hHi(ncat);
  for (size_t k = 0; k < ncat; ++k) {
    hLo[k] = hicat[k] + dhiMin;
    hHi[k] = hicat[k + 1] - dhiMin;
    if (hHi[k] <= hLo[k]) {
      // Bin too narrow for the requested gap; fall back to bin centre.
      const double mid = 0.5 * (hicat[k] + hicat[k + 1]);
      hLo[k] = hHi[k] = mid;
    }
  }

  // Feasibility: viceTarget must lie within [sum aicen*hLo, sum aicen*hHi].
  double vLo = 0.0, vHi = 0.0;
  for (size_t k = 0; k < ncat; ++k) {
    if (aicen[k] > Constants::puny) {
      vLo += aicen[k] * hLo[k];
      vHi += aicen[k] * hHi[k];
    }
  }
  if (viceTarget < vLo - 1.0e-12 || viceTarget > vHi + 1.0e-12) {
    return false;
  }

  // Initial per-cat thickness, clamped to bounds.
  std::vector<double> h(ncat, 0.0);
  for (size_t k = 0; k < ncat; ++k) {
    if (aicen[k] > Constants::puny) {
      const double hCur = vicen[k] / aicen[k];
      h[k] = std::min(std::max(hCur, hLo[k]), hHi[k]);
    }
  }

  // Distribute residual proportionally to available headroom (or shortfall).
  for (int iter = 0; iter < 50; ++iter) {
    double v = 0.0;
    for (size_t k = 0; k < ncat; ++k) v += aicen[k] * h[k];
    const double resid = viceTarget - v;
    if (std::fabs(resid) < 1.0e-14) break;

    // dh[k] should satisfy sum_k aicen[k] * dh[k] = resid, while respecting
    // bounds. Distribute proportionally to per-cat thickness room: choose
    // dh[k] = resid * room[k] / sum_k(aicen[k] * room[k]). That way the
    // per-cat aicen-weighting cancels in the volume sum.
    double weight = 0.0;
    for (size_t k = 0; k < ncat; ++k) {
      if (aicen[k] <= Constants::puny) continue;
      const double room = (resid > 0) ? (hHi[k] - h[k]) : (h[k] - hLo[k]);
      weight += aicen[k] * std::max(room, 0.0);
    }
    if (weight <= 0.0) return false;

    for (size_t k = 0; k < ncat; ++k) {
      if (aicen[k] <= Constants::puny) continue;
      const double room = (resid > 0) ? (hHi[k] - h[k]) : (h[k] - hLo[k]);
      const double dh   = resid * (std::max(room, 0.0) / weight);
      h[k] = std::min(std::max(h[k] + dh, hLo[k]), hHi[k]);
    }
  }

  for (size_t k = 0; k < ncat; ++k) {
    vicen[k] = (aicen[k] > Constants::puny) ? aicen[k] * h[k] : 0.0;
  }
  return true;
}

// ---------------------------------------------------------------------------
// snowIceFreeboard
// ---------------------------------------------------------------------------
std::vector<double> snowIceFreeboard(const std::vector<double> & aicen,
                                     const std::vector<double> & vicen,
                                     const std::vector<double> & vsnon,
                                     double rhoIce,
                                     double rhoSnow,
                                     double rhoOcean) {
  const size_t ncat = aicen.size();
  std::vector<double> fb(ncat, 0.0);
  const double cIce  = (rhoOcean - rhoIce) / rhoOcean;
  const double cSnow = rhoSnow / rhoOcean;
  for (size_t k = 0; k < ncat; ++k) {
    if (aicen[k] <= Constants::puny) continue;
    const double hi = vicen[k] / aicen[k];
    const double hs = vsnon[k] / aicen[k];
    fb[k] = hi * cIce - hs * cSnow;
  }
  return fb;
}

// ---------------------------------------------------------------------------
// enforceFreeboard
//
// Pass 1: redistribute snow. For each donor cat with fb<0, compute the snow
// volume that would bring it to fb=0. Move that excess to receiver cats with
// fb>0 in proportion to their freeboard headroom (weighted by aicen).
// If donors > receivers' total headroom, all snow is shoved onto the receivers
// up to the limit and the rest stays put.
// Pass 2: any residual fb<0 is fixed by growing ice volume in that cat just
// enough to make fb >= 0. This nudges aice/vice slightly (the freeboard-vs-
// ITD-bounds tradeoff noted in the Python script).
// ---------------------------------------------------------------------------
bool enforceFreeboard(std::vector<double> & aicen,
                      std::vector<double> & vicen,
                      std::vector<double> & vsnon,
                      double rhoIce,
                      double rhoSnow,
                      double rhoOcean) {
  const size_t ncat = aicen.size();
  if (ncat == 0 || vicen.size() != ncat || vsnon.size() != ncat) return false;
  const double cIce  = (rhoOcean - rhoIce) / rhoOcean;
  const double cSnow = rhoSnow / rhoOcean;

  // ---- Pass 1: redistribute snow off flooded cats ----
  // Compute per-cat fb and the snow volume change needed to bring fb to 0
  // (positive = excess to remove, negative = headroom to absorb).
  std::vector<double> fb = snowIceFreeboard(aicen, vicen, vsnon,
                                            rhoIce, rhoSnow, rhoOcean);

  double excessTotal = 0.0;
  std::vector<double> headroom(ncat, 0.0);   // positive: how much snow can be added
  std::vector<double> excess(ncat, 0.0);     // positive: how much snow must leave
  for (size_t k = 0; k < ncat; ++k) {
    if (aicen[k] <= Constants::puny) continue;
    // dh_snow that drives fb to 0: dhs = fb / cSnow  (subtract this from hs).
    const double dhs = fb[k] / cSnow;
    if (dhs < 0.0) {
      // fb < 0: need to remove |dhs| snow thickness, capped by what's there.
      const double hs   = vsnon[k] / aicen[k];
      const double take = std::min(-dhs, hs);
      excess[k]    = take * aicen[k];
      excessTotal += excess[k];
    } else {
      // fb > 0: cat can absorb up to dhs more snow thickness.
      headroom[k] = dhs * aicen[k];
    }
  }

  if (excessTotal > 0.0) {
    double headroomTotal = 0.0;
    for (size_t k = 0; k < ncat; ++k) headroomTotal += headroom[k];

    if (headroomTotal > 0.0) {
      const double absorbable = std::min(excessTotal, headroomTotal);
      // Remove excess proportionally from donors.
      const double removeFrac = absorbable / excessTotal;
      for (size_t k = 0; k < ncat; ++k) {
        if (excess[k] > 0.0) vsnon[k] -= excess[k] * removeFrac;
      }
      // Add to receivers proportionally to headroom.
      for (size_t k = 0; k < ncat; ++k) {
        if (headroom[k] > 0.0) {
          vsnon[k] += absorbable * (headroom[k] / headroomTotal);
        }
      }
    }
  }

  // ---- Pass 2: grow ice volume on any cat still flooded ----
  fb = snowIceFreeboard(aicen, vicen, vsnon, rhoIce, rhoSnow, rhoOcean);
  for (size_t k = 0; k < ncat; ++k) {
    if (aicen[k] <= Constants::puny) continue;
    if (fb[k] < 0.0) {
      // Solve hi_new*cIce = hs*cSnow  -> hi_new = hs*cSnow/cIce.
      const double hs     = vsnon[k] / aicen[k];
      const double hiNew  = hs * cSnow / cIce;
      vicen[k] = hiNew * aicen[k];
    }
  }

  // Final check.
  fb = snowIceFreeboard(aicen, vicen, vsnon, rhoIce, rhoSnow, rhoOcean);
  for (size_t k = 0; k < ncat; ++k) {
    if (aicen[k] > Constants::puny && fb[k] < -1.0e-9) return false;
  }

  // Suppress "unused" warning if compiled without the helper.
  (void)sum;
  return true;
}

}  // namespace icephysics
}  // namespace soca
