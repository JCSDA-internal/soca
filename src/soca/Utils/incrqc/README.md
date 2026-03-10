# `soca::incrqc` — Increment Quality Control

This module provides quality control (QC) functionality to ensure that ocean analysis increments remain within physically meaningful bounds.

The QC process ensures:
- The final analysis state stays within user-defined physical bounds.
- The increment does not add new unstable cells to the background.
- The steric height increment remains within specified limits.

## Features

- **Analysis Bounds Enforcement**: Adjusts increments so that the analysis stays within user-specified bounds for temperature and salinity.
- **Water Column Stability Check**: Ensures that the increment does not introduce new static instabilities in the water column by checking density profiles.
- **Steric Height Constraint**: Limits the steric height increment to prevent unrealistically large changes to sea surface height.
- **Coastal Increment Filter**: Smoothly tapers T/S increments to zero near coastlines using a cosine taper on the `distance_from_coast` field.
- **Iterative Refinement**: Applies stability checks iteratively with optional smoothing to refine corrections.

## Usage

The main entry point is `soca::incrqc::qcIncrement()`, which performs in-place QC on an increment field:

```cpp
#include "soca/Utils/incrqc/include/soca_incr_qc.h"

// Assuming xb (background), dx (increment), config, and geom are available
soca::incrqc::qcIncrement(xb, dx, config, geom);
```

## Configuration Example

```yaml
state bounds:
  sea_water_potential_temperature: [-2.5, 36.0]
  sea_water_salinity: [0.0, 44.0]

absolute steric increment max: 0.5  # [m]

steric variable change:
  linear variable changes:
  - linear variable change name: BalanceSOCA

stability check:
  increment stability iterations: 20
  min stable density gradient: 1.0e-4  # [kg/m³/m]

coastal increment filter:
  min distance: 0.0       # [m] zero increments at the coast
  max distance: 200000.0  # [m] full increments beyond 200 km from coast
```

---

## API Reference

### `double adjustAnalysisBounds(double xB, double dX, double minBound, double maxBound)`

Adjusts a scalar increment so that the analysis value does not exceed the defined bounds.

- **Parameters**:
  - `xB`: Background value.
  - `dX`: Proposed increment.
  - `minBound`: Minimum allowed analysis value.
  - `maxBound`: Maximum allowed analysis value.
- **Returns**: A possibly adjusted `dX` such that `xB + dX` is within bounds.

### `void qcIncrement(const soca::State& xb, soca::Increment& dx, const eckit::Configuration& config, const soca::Geometry& geom)`

Main routine that performs in-place quality control on the increment `dx`.

- **Parameters**:
  - `xb`: Background ocean state.
  - `dx`: Increment to be quality-controlled. Modified in place.
  - `config`: Configuration with QC parameters and bounds.
  - `geom`: Geometry providing spatial metadata.

---

## Detailed Algorithm Descriptions

### 1. Water Column Stability

This check ensures that the increment does not introduce **new static instabilities** into the water column.

#### Method

The check operates iteratively (with a user-defined number of iterations), and for each grid point (node) it:

1. **Computes density** using the UNESCO 1983 equation of state (Fofonoff & Millard, 1983) for both the background and the analysis (i.e., `background + increment`) at every vertical level.
2. **Calculates vertical density gradients** (∂ρ/∂z) for both background and analysis.
3. **Identifies instability conditions**:
   - The analysis introduces a static instability (∂ρ/∂z < 0) where the background was stable (∂ρ/∂z ≥ 0).
   - The analysis increases the level of instability already present in the background.
4. **Applies a local correction** by scaling down temperature and salinity increments at affected levels based on the ratio of the analysis gradient to a user-defined minimum stable density gradient (`min stable density gradient`). The unit-less correction factor is clamped between 0.1 and 1.0.
5. **Smooths corrected values** using neighboring points to reduce grid noise and introduce local consistency.
6. **Repeats from step 1** until the maximum number of iterations is reached.

#### Stability Correction Details

**Nomenclature:**
- $\rho^{\text{bkg}}_k = \rho(T^{\text{bkg}}_k, S^{\text{bkg}}_k)$: background density at level $k$
- $\rho^{\text{ana}}_k = \rho(T^{\text{bkg}}_k + \delta T_k, S^{\text{bkg}}_k + \delta S_k)$: analysis (background + increment) density at level $k$
- $z_k$: depth at level $k$ (positive downward)
- $\frac{\partial \rho^{\text{bkg}}}{\partial z} \big|_k = \frac{\rho^{\text{bkg}}_k - \rho^{\text{bkg}}_{k-1}}{z_k - z_{k-1}}$
- $\frac{\partial \rho^{\text{ana}}}{\partial z} \big|_k = \frac{\rho^{\text{ana}}_k - \rho^{\text{ana}}_{k-1}}{z_k - z_{k-1}}$
- $\rho_{z}^{\text{min}} = \frac{\rho_0 N^2}{g}$ where $N^2$ is the Brunt–Väisälä frequency for a weakly stratified ocean.

The increment is flagged as **potentially destabilizing** if either:

1. The background is stable:

$$\frac{\partial \rho^{\text{bkg}}}{\partial z} \bigg|_k \geq 0 \quad \text{and} \quad \frac{\partial \rho^{\text{ana}}}{\partial z} \bigg|_k < 0$$
2. The background is already unstable but the analysis makes it worse:

$$\frac{\partial \rho^{\text{bkg}}}{\partial z} \bigg|_k < 0 \quad \text{and} \quad \frac{\partial \rho^{\text{ana}}}{\partial z} \bigg|_k < \frac{\partial \rho^{\text{bkg}}}{\partial z} \bigg|_k$$

In these cases, a correction factor is applied to the temperature and salinity increments:

$$\delta T_k \leftarrow \delta T_k \cdot \left(1 - 0.5 \cdot \text{clamp}\left(\frac{\left|\frac{\partial \rho^{\text{ana}}}{\partial z}\right|_k}{\rho_{z}^{\text{min}}},\ 0.1,\ 1.0\right)\right)$$

$$\delta S_k \leftarrow \delta S_k \cdot \left(1 - 0.5 \cdot \text{clamp}\left(\frac{\left|\frac{\partial \rho^{\text{ana}}}{\partial z}\right|_k}{\rho_{z}^{\text{min}}},\ 0.1,\ 1.0\right)\right)$$

Finally, corrected values are optionally smoothed using neighbor averages:

$$\delta T_k \leftarrow \overline{\delta T}_k^{\text{neighbors}}$$

$$\delta S_k \leftarrow \overline{\delta S}_k^{\text{neighbors}}$$

#### Notes

- Depth increases positively downward, so stable stratification corresponds to ∂ρ/∂z>0.
- Corrections are only applied where bathymetry is positive and layer thickness is non-zero.
- Neighbor information is used to perform **localized smoothing** of the corrected increments for both temperature and salinity.
- This procedure is intended to maintain hydrostatic stability and avoid introducing artificial density inversions that can degrade the forecast.

### 2. Steric Height Limit

This check constrains the **sea surface height (SSH) increment** derived from temperature and salinity changes to remain within physically meaningful limits.

#### Method

For each node, the check:

1. **Evaluates the SSH increment** derived from temperature and salinity profiles.
2. **Checks if the absolute SSH increment exceeds** a configured maximum (`absolute steric increment max`).
3. If the threshold is exceeded:
   - **Computes a rescaling factor** based on the ratio of the maximum allowed SSH increment to the current value.
   - **Scales down the temperature and salinity increments** by this factor to reduce their effect on SSH.
4. **Recomputes the steric height increment** using the rescaled temperature and salinity profiles.
5. **Updates the SSH increment** to match the recomputed steric height.

#### Notes

- Layer thickness is used as a vertical integration weight in computing steric height.
- The routine supports optional debugging output, including:
  - Original and rescaled SSH increment
  - Computed steric height using 3 alternative formulations
- This constraint ensures the increment does not introduce unrealistically large SSH adjustments that could degrade ocean model balance or lead to unrealistic surface gravity waves.

### 3. Coastal Increment Filter

This filter smoothly tapers temperature and salinity increments to zero near coastlines, where the ocean model grid is coarser relative to the dynamics and where observation coverage is often sparse, leading to potentially unreliable increments.

The filter uses the precomputed `distance_from_coast` field from the background state.

#### Configuration

```yaml
coastal increment filter:
  min distance: 0.0       # [m] zero increments at the coast
  max distance: 100000.0  # [m] full increments beyond 100 km from coast
```

#### Method

For each ocean node the distance from coast *d* is compared against the two thresholds:

- **d ≤ d_min**: the weight is 0 — increments are completely removed.
- **d_min < d < d_max**: a smooth cosine taper is applied:

$$w = \frac{1}{2}\left(1 - \cos\!\left(\pi\,\frac{d - d_{\min}}{d_{\max} - d_{\min}}\right)\right)$$

  This provides a C¹-smooth transition (zero derivative at both endpoints).
- **d ≥ d_max**: the weight is 1 — increments are left unchanged.

The weight *w* is applied uniformly to all vertical levels at a given node.

#### Notes

- The filter is optional; it is only activated when the `coastal increment filter` key is present in the configuration.
- The `distance_from_coast` field must be present in the background state.
- Applied before the stability and steric height checks so that downstream QC operates on already-tapered increments.

### 4. Hard Bounds Enforcement

- Brute-force check to make sure analysis values of temperature and salinity remain within defined minimum and maximum values.
- Adjusts the increment field accordingly.

---


## References

- Fofonoff, N. P., & Millard, R. C. (1983). Algorithms for computation of fundamental properties of seawater. UNESCO technical papers in marine science, 44, 53.
- Pedlosky, J. (1987). *Geophysical Fluid Dynamics*. Springer.
- Lellouche, J.-M. et al. (2018). Recent updates to the Copernicus Marine Service global ocean monitoring and forecasting system. *Ocean Sci.*, 14, 1093–1126.
