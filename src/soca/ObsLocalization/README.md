# Rossby observation localization (`localization method: Rossby`)

This page is user documentation for configuring Rossby-based horizontal localization in SOCA.

## Where to configure it

Configure Rossby localization under `obs localizations`:

```yaml
obs localizations:
- localization method: Rossby
  base value: 500.0e3
  rossby mult: 2.0
  min grid mult: 1.0
  min value: 200.0e3
  max value: 3000.0e3
```

## Required and optional parameters

- `localization method` (required): must be `Rossby`.
- `base value` (optional, default `0.0`): additive baseline length scale.
- `rossby mult` (optional, default `1.0`): multiplier on local `rossby_radius`.
- `min grid mult` (optional, default `1.0`): grid-based floor factor.
- `min value` (optional): lower bound on final length scale.
- `max value` (optional): upper bound on final length scale.

The computation uses local `rossby_radius` and cell `area` from geometry fields.

## How the localization length scale is computed

At each model location, SOCA computes:

$$L = \text{base value} + \text{rossby mult} \times \text{rossby radius}$$

where `rossby radius` is the local `rossby_radius` geometry field.

Then applies the grid-size floor:

$$L \leftarrow \max\left(L,\; \text{min grid mult} \times \sqrt{\text{area}}\right)$$

Then applies optional bounds:

$$L \leftarrow \max(L,\; \text{min value})\;\text{if set}$$
$$L \leftarrow \min(L,\; \text{max value})\;\text{if set}$$

The scale $L$ at this stage is a **Gaussian-style horizontal length scale**
(in meters): it is the Rossby/local-grid blended scale used before applying
the GC99 localization operator.

This bounded scale is then converted internally for GC99 localization width:

$$L_{GC99} = L \times \frac{2}{\sqrt{0.3}}$$

## Examples

### 1) Mostly Rossby-driven localization

```yaml
obs localizations:
- localization method: Rossby
  base value: 100.0e3
  rossby mult: 1.0
  min grid mult: 1.0
```

How this behaves:
- Scale follows local `rossby_radius` closely.
- `base value` keeps a minimum additive baseline.
- Grid floor still prevents unrealistically small scales on fine cells.

### 2) Strong lower/upper bounds

```yaml
obs localizations:
- localization method: Rossby
  base value: 0.0
  rossby mult: 2.0
  min grid mult: 2.0
  min value: 200.0e3
  max value: 900.0e3
```

How this behaves:
- Raw Rossby-based scale is amplified by `rossby mult`.
- Scale cannot go below `200 km` or above `900 km`.
- Grid floor can still dominate in regions where `sqrt(area)` is large.

### 3) Nearly constant-scale behavior

```yaml
obs localizations:
- localization method: Rossby
  base value: 600.0e3
  rossby mult: 0.0
  min grid mult: 1.0
```

How this behaves:
- Removes Rossby-radius dependence (`rossby mult: 0.0`).
- Uses approximately constant `600 km`, except where the grid floor forces larger values.

## Practical tips

- Increase `rossby mult` to make localization more flow-dependent.
- Increase `min grid mult` to avoid over-localization on small cells.
- Use `min value`/`max value` to control extremes and improve robustness.

## Solver compatibility
Rossby localization currently works only with LETKF / GETKF. It is not yet supported with sequential EnKF solvers (e.g. EAKF), which require obs–obs localization. `rossby_radius` is only defined on the model grid, not at arbitrary obs pairs. Invoking it from a sequential solver aborts at runtime.

## Example files in SOCA tests

- `soca/test/testinput/obslocalization.yml`
- `soca/test/testinput/letkf.yml`
- `soca/test/testinput/letkf3d.yml`
- `soca/test/testinput/letkf_split_solver.yml`
