/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

#include "soca/PostProcess/IcePhysics.h"

namespace {

int gFailures = 0;

#define EXPECT(cond)                                                       \
  do {                                                                     \
    if (!(cond)) {                                                         \
      std::fprintf(stderr, "FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); \
      ++gFailures;                                                         \
    }                                                                      \
  } while (0)

bool nearEq(double a, double b, double tol) {
  return std::fabs(a - b) <= tol;
}

double sumv(const std::vector<double> & v) {
  return std::accumulate(v.begin(), v.end(), 0.0);
}

// ---------------------------------------------------------------------------
// adjustThicknessCategories
// ---------------------------------------------------------------------------
void testAdjustThkncats_Identity() {
  // Already-consistent distribution: hin per cat inside [hicat[k], hicat[k+1]].
  // viceTarget = current sum -> nothing should change appreciably.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.10, 0.20, 0.15, 0.05, 0.00};
  std::vector<double> vicen{0.03, 0.20, 0.30, 0.18, 0.00};  // hins ~ 0.3,1.0,2.0,3.6
  const double vTarget = sumv(vicen);
  const double aTotal  = sumv(aicen);

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, vTarget, hicat, 0.01);
  EXPECT(ok);
  EXPECT(nearEq(sumv(vicen), vTarget, 1.0e-12));
  EXPECT(nearEq(sumv(aicen), aTotal,  1.0e-12));
  for (size_t k = 0; k < aicen.size(); ++k) {
    if (aicen[k] > 1e-11) {
      const double h = vicen[k] / aicen[k];
      EXPECT(h >= hicat[k]     - 1.0e-9);
      EXPECT(h <= hicat[k + 1] + 1.0e-9);
    }
  }
}

void testAdjustThkncats_OutOfBoundClamps() {
  // Cat 0 has hin=2.0 (way above its bin [0, 0.64]) -> must be reduced.
  // Use a feasible target inside the post-clamp envelope.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.20, 0.10, 0.10, 0.05, 0.00};
  std::vector<double> vicen{0.40, 0.10, 0.20, 0.18, 0.00};
  const double aTotal  = sumv(aicen);
  // Post-clamp envelope is [vLo,vHi] ~ [0.452, 0.738]. Pick something inside.
  const double vTarget = 0.60;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, vTarget, hicat, 0.01);
  EXPECT(ok);
  EXPECT(nearEq(sumv(vicen), vTarget, 1.0e-9));
  EXPECT(nearEq(sumv(aicen), aTotal,  1.0e-12));
  for (size_t k = 0; k < aicen.size(); ++k) {
    if (aicen[k] > 1e-11) {
      const double h = vicen[k] / aicen[k];
      EXPECT(h >= hicat[k]     - 1.0e-9);
      EXPECT(h <= hicat[k + 1] + 1.0e-9);
    }
  }
}

void testAdjustThkncats_HitTarget() {
  // Force the algorithm to redistribute volume away from initial layout.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.10, 0.10, 0.10, 0.10, 0.10};
  std::vector<double> vicen{0.05, 0.10, 0.20, 0.30, 0.50};  // sum 1.15
  // Per-cat envelope sums to vLo ~ 0.911, vHi ~ 5.901. Pick inside.
  const double vTarget = 1.30;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, vTarget, hicat, 0.01);
  EXPECT(ok);
  EXPECT(nearEq(sumv(vicen), vTarget, 1.0e-9));
}

void testAdjustThkncats_InfeasibleTarget() {
  // Demand a viceTarget below the lower envelope -> should report failure.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.0, 0.0, 0.0, 0.20, 0.0};   // only cat 3 has area
  std::vector<double> vicen{0.0, 0.0, 0.0, 0.50, 0.0};
  // cat 3 lower bound = (2.47+0.01)*0.20 = 0.496; target below that.
  const double vTarget = 0.10;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, vTarget, hicat, 0.01);
  EXPECT(!ok);
}

// ---------------------------------------------------------------------------
// enforceFreeboard
// ---------------------------------------------------------------------------
void testFreeboard_AlreadyOk() {
  // Healthy column: thick ice, modest snow -> all freeboards positive.
  std::vector<double> aicen{0.10, 0.20, 0.15};
  std::vector<double> vicen{0.05, 0.40, 0.50};   // hi = 0.5, 2.0, 3.33
  std::vector<double> vsnon{0.005, 0.02, 0.015}; // hs = 0.05, 0.1, 0.1
  auto a0 = aicen, v0 = vicen, s0 = vsnon;

  const bool ok = soca::icephysics::enforceFreeboard(aicen, vicen, vsnon);
  EXPECT(ok);
  // Should be untouched (no flooded cats).
  for (size_t k = 0; k < 3; ++k) {
    EXPECT(nearEq(aicen[k], a0[k], 1.0e-12));
    EXPECT(nearEq(vicen[k], v0[k], 1.0e-12));
    EXPECT(nearEq(vsnon[k], s0[k], 1.0e-12));
  }
}

void testFreeboard_RedistributeSnow() {
  // One thin/heavy-snow cat (flooded) and one thick cat with headroom.
  // rho_ice=917, rho_snow=330, rho_ocean=1025
  // cat 0: aicen=0.5, hi=0.10, hs=0.30  -> fb = 0.10*(108/1025) - 0.30*(330/1025)
  //        = 0.01054 - 0.0966 = -0.086 (flooded)
  // cat 1: aicen=0.5, hi=2.00, hs=0.05  -> fb = 2.0*0.1054 - 0.05*0.3220
  //        = 0.2107 - 0.01610 = 0.195 (lots of headroom)
  std::vector<double> aicen{0.5, 0.5};
  std::vector<double> vicen{0.05, 1.00};
  std::vector<double> vsnon{0.15, 0.025};
  const double snowTotal0 = sumv(vsnon);

  const bool ok = soca::icephysics::enforceFreeboard(aicen, vicen, vsnon);
  EXPECT(ok);
  // Total snow must be conserved (no ice modification needed in this config).
  EXPECT(nearEq(sumv(vsnon), snowTotal0, 1.0e-12));
  // All freeboards now non-negative.
  auto fb = soca::icephysics::snowIceFreeboard(aicen, vicen, vsnon);
  for (double f : fb) EXPECT(f >= -1.0e-9);
  // Ice was not modified in this regime.
  EXPECT(nearEq(vicen[0], 0.05, 1.0e-12));
  EXPECT(nearEq(vicen[1], 1.00, 1.0e-12));
}

void testFreeboard_GrowIce() {
  // Single-category column, deeply flooded -> can't redistribute snow,
  // so ice volume must grow to restore freeboard.
  std::vector<double> aicen{1.0};
  std::vector<double> vicen{0.05};   // hi = 0.05
  std::vector<double> vsnon{0.30};   // hs = 0.30 -> fb < 0

  const bool ok = soca::icephysics::enforceFreeboard(aicen, vicen, vsnon);
  EXPECT(ok);
  auto fb = soca::icephysics::snowIceFreeboard(aicen, vicen, vsnon);
  EXPECT(fb[0] >= -1.0e-9);
  // Specifically hi_new = hs * rho_snow / (rho_ocean - rho_ice)
  //                     = 0.30 * 330 / 108 = 0.9167
  EXPECT(nearEq(vicen[0], 0.30 * 330.0 / (1025.0 - 917.0), 1.0e-9));
}

}  // namespace

int main() {
  testAdjustThkncats_Identity();
  testAdjustThkncats_OutOfBoundClamps();
  testAdjustThkncats_HitTarget();
  testAdjustThkncats_InfeasibleTarget();
  testFreeboard_AlreadyOk();
  testFreeboard_RedistributeSnow();
  testFreeboard_GrowIce();

  if (gFailures == 0) {
    std::printf("TestIcePhysics: all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "TestIcePhysics: %d failure(s)\n", gFailures);
  return 1;
}
