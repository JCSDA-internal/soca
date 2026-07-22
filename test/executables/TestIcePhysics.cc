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
      aicen, vicen, aTotal, vTarget, hicat);
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
      aicen, vicen, aTotal, vTarget, hicat);
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
  const double aTotal  = sumv(aicen);
  // Per-cat envelope sums to vLo ~ 0.911, vHi ~ 5.901. Pick inside.
  const double vTarget = 1.30;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aTotal, vTarget, hicat);
  EXPECT(ok);
  EXPECT(nearEq(sumv(vicen), vTarget, 1.0e-9));
}

void testAdjustThkncats_InfeasibleTarget() {
  // Demand a viceTarget far below the feasibility envelope -> should fail
  // even after Phase 2 (cannot reshuffle aicen to thin enough).
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.0, 0.0, 0.0, 0.20, 0.0};   // only cat 3 has area
  std::vector<double> vicen{0.0, 0.0, 0.0, 0.50, 0.0};
  const double aTarget = sumv(aicen);
  // After reshuffling aicen to the thinnest cat, the lower envelope is still
  // ~ aTarget * hLo[0] = 0.20 * 0.01 = 0.002. Target = 0.001 is below that.
  const double vTarget = 1.0e-3;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aTarget, vTarget, hicat);
  EXPECT(!ok);
}

// ---------------------------------------------------------------------------
// Full-form tests: aiceTarget + viceTarget, early-out, Phase 2 (aicen mutation)
// ---------------------------------------------------------------------------
void testAdjustThkncats_EarlyOut() {
  // Distribution already satisfies bin bounds AND totals match the targets.
  // Result must be bit-equal to input; aicenMutated must remain 0.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.10, 0.20, 0.15, 0.05, 0.00};
  std::vector<double> vicen{0.03, 0.20, 0.30, 0.18, 0.00};
  const auto aIn = aicen, vIn = vicen;
  const double aT = sumv(aicen);
  const double vT = sumv(vicen);
  int mutated = -1;
  double dMax = -1.0;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aT, vT, hicat, 0.01, 1.0e-8, 1.0e-8, 1.0e-6,
      &mutated, &dMax);
  EXPECT(ok);
  EXPECT(mutated == 0);
  EXPECT(dMax >= 0.0 && dMax < 1.0e-12);
  for (size_t k = 0; k < aicen.size(); ++k) {
    EXPECT(nearEq(aicen[k], aIn[k], 1.0e-12));
    EXPECT(nearEq(vicen[k], vIn[k], 1.0e-12));
  }
}

void testAdjustThkncats_Phase1Only() {
  // Bins violated (cat 0 too thick) but viceTarget reachable while holding
  // aicen fixed. Expect Phase 1 to fix it; mutated flag stays 0.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.20, 0.10, 0.10, 0.05, 0.00};
  std::vector<double> vicen{0.40, 0.10, 0.20, 0.18, 0.00};  // cat 0 h=2.0 out
  const auto aIn = aicen;
  const double aT = sumv(aicen);
  const double vT = 0.60;  // inside the post-clamp envelope
  int mutated = -1;
  double dMax = -1.0;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aT, vT, hicat, 0.01, 1.0e-8, 1.0e-8, 1.0e-6,
      &mutated, &dMax);
  EXPECT(ok);
  EXPECT(mutated == 0);
  // aicen only mutated by the rescale-to-target (which is a no-op here since
  // Sum already == aT). Tolerance generous to absorb the cleanup rescale.
  EXPECT(dMax < 1.0e-9);
  EXPECT(nearEq(sumv(aicen), aT, 1.0e-12));
  EXPECT(nearEq(sumv(vicen), vT, 1.0e-9));
  for (size_t k = 0; k < aicen.size(); ++k) {
    if (aicen[k] > 1.0e-11) {
      const double h = vicen[k] / aicen[k];
      EXPECT(h >= hicat[k]     - 1.0e-9);
      EXPECT(h <= hicat[k + 1] + 1.0e-9);
    }
  }
  // Per-cat aicen identical to input (no Phase 2).
  for (size_t k = 0; k < aicen.size(); ++k) {
    EXPECT(nearEq(aicen[k], aIn[k], 1.0e-9));
  }
}

void testAdjustThkncats_Phase2_NoiceToIce() {
  // Python script: noice->ice cell with placeholder aicen=[a_an, 0,0,0,0]
  // and a thick donor mean of 2.0 m -> viceTarget = 0.30 * 2.0 = 0.60.
  // Cat-0 envelope alone is [0.001, ~0.192], so Phase 1 cannot reach 0.60.
  // Phase 2 must shift aicen into cat 3.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.30, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> vicen{0.0,  0.0, 0.0, 0.0, 0.0};
  const double aT = 0.30;
  const double vT = 0.30 * 2.0;
  int mutated = -1;
  double dMax = -1.0;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aT, vT, hicat, 0.01, 1.0e-8, 1.0e-8, 1.0e-6,
      &mutated, &dMax);
  EXPECT(ok);
  EXPECT(mutated == 1);
  EXPECT(dMax > 0.0);
  EXPECT(nearEq(sumv(aicen), aT, 1.0e-9));
  EXPECT(nearEq(sumv(vicen), vT, 1.0e-9));
  // Some mass must have left cat 0 (the placeholder) since cat-0's envelope
  // ceiling 0.30 * 0.63 = 0.189 is below the target 0.60.
  EXPECT(aicen[0] < 0.30 - 1.0e-6);
  // All occupied cats must respect bin bounds.
  for (size_t k = 0; k < aicen.size(); ++k) {
    if (aicen[k] > 1.0e-11) {
      const double h = vicen[k] / aicen[k];
      EXPECT(h >= hicat[k]     - 1.0e-9);
      EXPECT(h <= hicat[k + 1] + 1.0e-9);
    }
  }
}

void testAdjustThkncats_TinyAice_Diagnostic() {
  // Reproduces a marginal-ice case from the rebin test grid:
  //   ai_an = 1.22e-5, hi_an = 2.415 m  -> vtot_target = 2.95e-5
  // Seed: aicen = [ai_an, 0, 0, 0, 0]; vicen = [vtot_target, 0, 0, 0, 0].
  // The target lies in cat 3's bin ([2.47, 4.57]); rebin must reshuffle
  // aicen there.
  std::vector<double> hicat{0.0, 0.6445072, 1.391433, 2.470179, 4.567288, 9.333887};
  std::vector<double> aicen{1.22e-5, 0, 0, 0, 0};
  std::vector<double> vicen{2.95e-5, 0, 0, 0, 0};
  const double aT = 1.22e-5;
  const double vT = 2.95e-5;
  int mut = -1;
  double dmax = 0.0;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aT, vT, hicat, 0.01, 1.0e-8, 1.0e-8, 1.0e-6, &mut, &dmax);
  // Diagnostic dump regardless of pass/fail.
  double aSum = 0.0, vSum = 0.0;
  for (size_t k = 0; k < 5; ++k) { aSum += aicen[k]; vSum += vicen[k]; }
  std::fprintf(stderr,
      "TinyAice: ok=%d mutated=%d  Sum(aicen)=%.4e (target %.4e)  Sum(vicen)=%.4e (target %.4e)\n",
      ok, mut, aSum, aT, vSum, vT);
  for (size_t k = 0; k < 5; ++k) {
    if (aicen[k] > 0.0) {
      std::fprintf(stderr, "  cat %zu: a=%.4e v=%.4e h=%.4f (bin [%.4f, %.4f])\n",
                   k, aicen[k], vicen[k], vicen[k]/aicen[k], hicat[k], hicat[k+1]);
    }
  }
  EXPECT(ok);
  EXPECT(nearEq(aSum, aT, 1.0e-12));
  EXPECT(nearEq(vSum, vT, 1.0e-9));
  for (size_t k = 0; k < 5; ++k) {
    if (aicen[k] > 1.0e-11) {
      const double h = vicen[k] / aicen[k];
      EXPECT(h >= hicat[k]     - 1.0e-9);
      EXPECT(h <= hicat[k + 1] + 1.0e-9);
    }
  }
}

void testAdjustThkncats_RescaleToAiceTarget() {
  // Caller supplies aiceTarget != Sum(aicen_in): expect rescale-to-target
  // (NOT counted as Phase 2 mutation), then redistribute vicen.
  std::vector<double> hicat{0.0, 0.64, 1.39, 2.47, 4.57, 50.0};
  std::vector<double> aicen{0.10, 0.20, 0.15, 0.05, 0.00};
  std::vector<double> vicen{0.03, 0.20, 0.30, 0.18, 0.00};
  const double aIn = sumv(aicen);
  const double aT = 0.60;                // bump total aice from 0.50 to 0.60
  const double vT = sumv(vicen) * (aT / aIn);  // proportional volume scale
  int mutated = -1;
  double dMax = -1.0;

  const bool ok = soca::icephysics::adjustThicknessCategories(
      aicen, vicen, aT, vT, hicat, 0.01, 1.0e-8, 1.0e-8, 1.0e-6,
      &mutated, &dMax);
  EXPECT(ok);
  EXPECT(mutated == 0);  // pure rescale; Phase 1 sees bins OK after scaling
  EXPECT(nearEq(sumv(aicen), aT, 1.0e-12));
  EXPECT(nearEq(sumv(vicen), vT, 1.0e-9));
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

// ---------------------------------------------------------------------------
// thermo helpers
// ---------------------------------------------------------------------------
void testIceEnthalpyBL99_RoundTripSign() {
  // qin must be negative (energy required to melt ice is "missing energy"
  // by sign convention). Hold S fixed, sweep T.
  const double S = 3.0;
  for (double T = -20.0; T < -0.5; T += 1.0) {
    const double q = soca::icephysics::iceEnthalpyBL99(T, S);
    EXPECT(q < 0.0);
  }
}

void testIceEnthalpyBL99_KnownValue() {
  // Reference numerics (computed offline):
  //   T = -10 C, S = 3 ppt, Tmlt = -0.054*3 = -0.162 C
  //   q = -917 * ( 2106*(Tmlt - T) + Lfresh*(1 - Tmlt/T) - 4218*Tmlt )
  //     = -917 * ( 2106*9.838 + 0.334e6*(1 - 0.0162) - 4218*(-0.162) )
  //     ~ -3.224e8 J/m^3
  const double q = soca::icephysics::iceEnthalpyBL99(-10.0, 3.0);
  EXPECT(q < -3.20e8);
  EXPECT(q > -3.30e8);
}

void testSiceLayerCice4_BoundsAndShape() {
  // Profile must be monotone non-decreasing for the standard 7-layer case
  // (ice gets saltier with depth) and bounded above by saltmax.
  const int nlyr = 7;
  double prev = 0.0;
  for (int k = 1; k <= nlyr; ++k) {
    const double S = soca::icephysics::siceLayerCice4(k, nlyr);
    EXPECT(S >= 0.0);
    EXPECT(S <= soca::icephysics::Constants::saltmax + 1.0e-12);
    if (k > 1) EXPECT(S >= prev - 1.0e-12);
    prev = S;
  }
  // Out-of-range indices return 0 cleanly.
  EXPECT(soca::icephysics::siceLayerCice4(0, 7) == 0.0);
  EXPECT(soca::icephysics::siceLayerCice4(8, 7) == 0.0);
  EXPECT(soca::icephysics::siceLayerCice4(1, 0) == 0.0);
}

void testSnowEnthalpy_KnownValue() {
  // Tsfc = 0 gives qsn = -rho_snow * Lfresh exactly.
  const double q0 = soca::icephysics::snowEnthalpy(0.0);
  EXPECT(nearEq(q0,
                -soca::icephysics::Constants::rho_snow
                * soca::icephysics::Constants::Lfresh,
                1.0e-6));
  // Tsfc = -10 C gives qsn = -rho_snow*(Lfresh - cp_ice*(-10)) =
  //                     -rho_snow*Lfresh - rho_snow*cp_ice*10
  // (sign: colder snow has more-negative enthalpy).
  const double T = -10.0;
  const double qExpected = -soca::icephysics::Constants::rho_snow
      * (soca::icephysics::Constants::Lfresh
         - soca::icephysics::Constants::cp_ice * T);
  EXPECT(nearEq(soca::icephysics::snowEnthalpy(T), qExpected, 1.0e-6));
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
  testAdjustThkncats_EarlyOut();
  testAdjustThkncats_Phase1Only();
  testAdjustThkncats_Phase2_NoiceToIce();
  testAdjustThkncats_TinyAice_Diagnostic();
  testAdjustThkncats_RescaleToAiceTarget();
  testFreeboard_AlreadyOk();
  testFreeboard_RedistributeSnow();
  testFreeboard_GrowIce();
  testIceEnthalpyBL99_RoundTripSign();
  testIceEnthalpyBL99_KnownValue();
  testSiceLayerCice4_BoundsAndShape();
  testSnowEnthalpy_KnownValue();

  if (gFailures == 0) {
    std::printf("TestIcePhysics: all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "TestIcePhysics: %d failure(s)\n", gFailures);
  return 1;
}
