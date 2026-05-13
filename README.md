# CSMQNM.jl

Authors: Shoya Ogawa and Okuto Morikawa

`CSMQNM.jl` is a cleaned-up Julia package skeleton for computing black-hole quasinormal-mode candidates with the complex scaling method (CSM). The implementation is intentionally narrowed to the basis systems currently needed in this project.

- `ical = 2` only, for legacy compatibility.
- Polynomial x Gaussian basis only.
- `range = :real` or `range = :complex` switches between polynomial x real-range Gaussian and polynomial x complex-range Gaussian.
- For `range = :complex`, only the `1 + im*beta` branch is included. The `1 - im*beta` branch is deliberately excluded because it is not convergent in the intended complex-scaled sector.
- Real-range Gaussian alone, sin/cos x real-range Gaussian, and CLD calculation are not included.
- Physical parameters are passed with a dictionary, and `ipot` is specified by a symbolic label rather than by a hard-to-remember integer.

## Installation

From the parent directory of this package:

```julia
using Pkg
Pkg.activate("CSMQNM")
Pkg.instantiate()
```

If `AutoTortoise.jl` is not available from your package environment, add it in the same way you used it in the original scripts. It is required for the dS/(A)dS tortoise-coordinate inversion used by `:sds4_rw` and `:sds_dim_rw`.

## Minimal example

```julia
using CSMQNM

cfg = RunConfig(
    basis = BasisConfig(imax=30, nmax=3, r0=0.1, rmax=60.0, range=:real),
    csm = CSMConfig(the0=42.0),
    integration = IntegrationConfig(xmin=-60.0, xmax=80.0, dx=0.02),
    physics = Dict(
        :ipot => :sds_dim_rw,
        :M => 1.0,
        :ell => 2,
        :s => 0,
        :lam => 0.4,
        :dim => 5,
        :c => :tensor,
    ),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm(cfg)
println(result.omega)
```

## Dictionary-style input

For notebook/script use, the dictionary wrapper is usually the quickest interface:

```julia
result = solve_qnm_from_dict(
    basis = Dict(
        :imax=>30,
        :nmax=>3,
        :r0=>0.1,
        :rmax=>60.0,
        :range=>:complex,
        :beta=>π/2,
    ),
    csm = Dict(:the0=>42.0),
    integration = Dict(:xmin=>-60.0, :xmax=>80.0, :dx=>0.02),
    physics = Dict(:ipot=>:sds4_rw, :M=>1.0, :ell=>2, :s=>2, :lam=>0.08),
    output = Dict(:write_spectrum=>true),
)
```

## Complex range Gaussian branch

For `range = :complex`, the Gaussian range is defined by

```julia
alpha_i = (1 + im * beta) / r_i^2
```

Only this branch is used for all polynomial orders `n = 0, ..., nmax`. The conjugate branch

```julia
(1 - im * beta) / r_i^2
```

is not included.

The convergence check follows the damping condition for the factor

```text
exp[-alpha * (1 + i beta) * exp(-2 i theta) * r_*^2]
```

and requires

```text
cos(2 theta) + beta sin(2 theta) > 0.
```

The package throws an error when this condition is violated. In particular:

- `beta > 0` is required.
- `0 < theta < pi/2` is required.
- For `pi/4 < theta < pi/2`, `beta` must be large enough to satisfy `beta > -cot(2 theta)`.

The scaling angle input `the0` is in degrees, so the above check is applied after converting `theta = the0 * pi / 180`.

For `range = :complex`, the package keeps the Polynomial x complex-range Gaussian basis with dimension

```julia
imax * (nmax + 1)
```

and solves the generalized eigenvalue problem directly:

```julia
H * c = E * N * c
```

The real-range `orthogonalize_basis` step is deliberately skipped for this branch because the complex-range overlap matrix is complex symmetric rather than Hermitian.


## Potential labels

Use symbolic labels in `physics[:ipot]`:

| label | meaning |
| --- | --- |
| `:schwarzschild_rw` | 4d Schwarzschild Regge-Wheeler potential |
| `:ern_scalar` | extremal Reissner-Nordstrom scalar potential |
| `:ern_em` | extremal Reissner-Nordstrom odd electromagnetic potential |
| `:ern_grav` | extremal Reissner-Nordstrom odd gravitational potential |
| `:sds4_rw` | 4d Schwarzschild-(A)dS Regge-Wheeler-type potential |
| `:sds_dim_rw` | higher-dimensional Schwarzschild-dS tensor/vector potential |
| `:stringy_ads_ds` | reserved label; coupled-channel backend is not implemented yet |

Legacy integer `ipot` values are still accepted internally:

```text
1 => :schwarzschild_rw
2 => :ern_scalar
3 => :ern_em
4 => :ern_grav
5 => :sds4_rw
6 => :stringy_ads_ds
7 => :sds_dim_rw
```

## Main parameter groups

`BasisConfig` controls the polynomial x Gaussian basis:

```julia
BasisConfig(imax=30, nmax=3, r0=0.1, rmax=60.0, range=:real)
BasisConfig(imax=30, nmax=3, r0=0.1, rmax=60.0, range=:complex, beta=π/2)
```

`CSMConfig` controls the complex-scaling angle:

```julia
CSMConfig(the0=42.0)
```

`IntegrationConfig` controls potential-matrix quadrature:

```julia
IntegrationConfig(xmin=-60.0, xmax=80.0, dx=0.02)
```

`physics` is a dictionary:

```julia
Dict(:ipot=>:sds_dim_rw, :M=>1.0, :ell=>2, :s=>0, :lam=>0.4, :dim=>5, :c=>:tensor)
```

For `:sds_dim_rw`, `:c=>:vector` and `:c=>:tensor` are aliases for `1` and `2`.

## Output

`solve_qnm` returns a `QNMResult`:

```julia
result.energy              # eigenvalues E = omega^2
result.omega               # QNM branch selected with Im(omega) <= 0
result.hamiltonian         # orthogonalized H for :real, raw H for :complex
result.basis_size_before
result.basis_size_after
result.dropped_basis_vectors
```

Set `OutputConfig(write_spectrum=true)` to write a spectrum file, and `OutputConfig(write_potential=true)` to write the potential sampled on the integration grid.

## Notes on the current scope

This package skeleton keeps the single-channel QNM computation only. The old CLD code, hand-picked CLD pole indices, stand-alone real-range Gaussian branch, complex-range Gaussian conjugate branch, and trigonometric basis branch were deliberately removed. The complex-range option is Polynomial x complex-range Gaussian, not the old even/odd-only complex Gaussian branch. The `:stringy_ads_ds` label remains as a placeholder because the original implementation is coupled-channel and should be refactored separately rather than forced into the single-channel API.

## ECS finite-difference backend

In addition to the original Gaussian-basis CSM backend, the package now contains an experimental ECS backend in `src/ECS.jl`.  The existing CSM code is kept in `src/CSM.jl`, and `src/CSMQNM.jl` only loads the two backends.

The ECS backend is intended for exploratory calculations beyond the practical `theta < pi/4` limitation of global real-range Gaussian matrix elements.  It discretizes the Schrödinger operator directly on an exterior-complex-scaled contour,

```text
z = g(x),
```

where `x` is a real finite-difference grid and `z` is the complex tortoise coordinate.  The finite-difference kinetic operator is built from

```text
- J^{-1} d/dx [ J^{-1} d/dx ],   J = dg/dx.
```

In `src/ECS.jl` this is discretized in conservative three-point flux form on the interior grid points.  The flux coefficient `J^{-1}` is evaluated at half-grid midpoints, so the matrix is tridiagonal rather than a product of two first-derivative matrices.  This is intended to avoid even/odd grid decoupling artifacts in the ECS continuum.

The potential is evaluated as `V(z)` through the same `potential_value` interface used by the CSM backend.  For dS/(A)dS potentials this calls `AutoTortoise.inverse_tortoise`, so this backend assumes that `AutoTortoise.jl` correctly handles the required complex inverse tortoise map along the chosen ECS contour.

A minimal Schwarzschild ECS example is

```julia
using CSMQNM

cfg = ECSRunConfig(
    ecs = ECSConfig(
        theta_deg = 60.0,
        xmin = -120.0,
        xmax = 120.0,
        nx = 1201,
        x0_left = 40.0,
        x0_right = 40.0,
    ),
    physics = Dict(
        :ipot => :schwarzschild_rw,
        :M => 1.0,
        :ell => 2,
    ),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm_ecs(cfg)
println(result.omega)
```

`ECSConfig` currently supports Dirichlet boundary conditions at the finite-box endpoints.  This is natural for ECS because outgoing waves are damped in the complex-scaled exterior regions when the box is sufficiently large.

The default contour is piecewise linear:

```text
g(x) = -x0_left  + (x + x0_left) exp(i theta),   x < -x0_left,
g(x) = x,                                         -x0_left <= x <= x0_right,
g(x) =  x0_right + (x - x0_right) exp(i theta),  x >  x0_right.
```

For preliminary smoothing tests, set `smoothing > 0`.  The smoothed contour uses a tanh-based ramp.  The piecewise-linear contour is easier to interpret and is the recommended first diagnostic.

The ECS result type is `ECSResult`:

```julia
result.energy       # eigenvalues E = omega^2
result.omega        # branch selected with Im(omega) <= 0
result.hamiltonian  # finite-difference ECS Hamiltonian
result.x            # interior real grid points
result.z            # interior ECS contour points
result.jacobian     # dg/dx on the interior grid
result.potential    # V(z) on the interior grid
```

For actual higher-overtone work, check stability under changes of `theta_deg`, `x0_left`, `x0_right`, `xmin`, `xmax`, and `nx`.  The ECS backend is deliberately separated from the Gaussian CSM backend so that this experimental development does not disturb the existing `solve_qnm` workflow.

## FEDVR-ECS backend

Version `0.1.0` also includes an experimental finite-element DVR backend on the
same exterior-complex-scaled contour.  This is intended as a higher-order
alternative to the uniform-grid finite-difference ECS solver.

The real computational coordinate is denoted by `x`, and the ECS contour is

```julia
z = g(x)
J = dg/dx
```

The Schrödinger operator is discretized in weak form,

```text
H = -d^2/dz^2 + V(z)
```

using local Gauss-Lobatto-Legendre DVR functions on each finite element.  In
terms of the real parameter `x`, the kinetic and overlap matrices are assembled
as

```text
K_ab = ∫ dx J(x)^(-1) dB_a/dx dB_b/dx,
N_ab = ∫ dx J(x) B_a(x) B_b(x),
V_ab = ∫ dx J(x) B_a(x) V(g(x)) B_b(x).
```

The DVR quadrature makes the potential and overlap matrices diagonal up to the
finite-element endpoint sharing.  Internal element boundaries are represented by
explicit bridge functions: the right endpoint of one element and the left
endpoint of the next element are merged into a single continuous basis function.
This enforces C0 continuity at element interfaces and avoids independent
left/right endpoint degrees of freedom.  Dirichlet boundary conditions are
imposed by removing the two external endpoint degrees of freedom.

A minimal example is

```julia
using CSMQNM

cfg = FEDVRECSRunConfig(
    ecs = FEDVRECSConfig(
        theta_deg = 60.0,
        xmin = -240.0,
        xmax = 240.0,
        nelements = 120,
        order = 10,
        x0_left = 25.0,
        x0_right = 25.0,
        smoothing = 2.5,
    ),
    physics = Dict(
        :ipot => :schwarzschild_rw,
        :M => 0.5,
        :ell => 3,
        :s => 2,
    ),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm_ecs_fedvr(cfg)
```

For comparison with the finite-difference ECS backend, `nelements = 120` and
`order = 10` gives roughly `120 * (10 - 1) + 1 = 1081` global grid points before
Dirichlet endpoint removal.  Increasing `order` improves local spectral
accuracy; increasing `nelements` improves spatial resolution and allows a more
local representation of the ECS transition region.

## Takagi-ECS Gaussian backend

Version `0.1.0` also includes an experimental Gaussian-packet ECS backend in
`src/TakagiECS.jl`.  This backend keeps the CSM/Gaussian-expansion flavor but
places localized Gaussian packets on the real ECS parameter `x`, rather than
using global Gaussian functions of the scaled coordinate.

**Status note.** Takagi-ECS is experimental. For the present Gaussian packet
basis, Takagi orthogonalization may fail because of severe conditioning issues.


The wavefunction is expanded as

```text
psi(x) = sum_i c_i phi_i(x),
phi_i(x) = envelope(x) exp[-(x - X_i)^2 / (2 sigma_i^2)].
```

The optional envelope is enabled by default and vanishes at `xmin` and `xmax`,
which imposes Dirichlet behavior at the outer box boundaries.  The ECS contour
is still

```text
z = g(x),    J = dg/dx.
```

The complex-symmetric weak-form matrices are assembled as

```text
N_ij = int dx J(x) phi_i(x) phi_j(x),
K_ij = int dx J(x)^(-1) phi_i'(x) phi_j'(x),
V_ij = int dx J(x) V(g(x)) phi_i(x) phi_j(x).
```

Because `N` is complex symmetric rather than Hermitian, the default
`orthogonalization = :auto` now uses a two-stage conditioning procedure:

```text
1. Hermitian auxiliary metric G_ij = int dx |J| conj(phi_i) phi_j
   removes near-linear dependencies of the Gaussian packets.
2. Takagi factorization is applied to the reduced complex-symmetric
   C-product overlap N_red.
```

The final transform still satisfies the C-product condition

```text
S^T N S ~= I,
H_ortho = S^T H S.
```

This is more stable than applying Takagi directly to the full Gaussian overlap.
For diagnostics one may set `orthogonalization = :takagi` for direct Takagi,
`:two_step` for the two-stage procedure, `:hermitian` for Hermitian auxiliary
conditioning only, or `:none` for the raw generalized eigenvalue problem.  In
`:auto`, if the two-stage Takagi step fails, the solver falls back to the
Hermitian-reduced generalized C-product problem.

A minimal example is

```julia
using CSMQNM

cfg = TakagiECSRunConfig(
    ecs = TakagiECSConfig(
        theta_deg = 60.0,
        xmin = -240.0,
        xmax = 240.0,
        nbasis = 240,
        sigma_scale = 1.5,
        x0_left = 30.0,
        x0_right = 60.0,
        smoothing = 5.0,
        quadrature_order = 12,
        overlap_cutoff = 1e-10,
        orthogonalization = :auto,
    ),
    physics = Dict(:ipot => :schwarzschild_rw, :M => 0.5, :ell => 3, :s => 2),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm_ecs_takagi(cfg)
```

Useful tuning parameters are `sigma_scale`, `nbasis`, `quadrature_order`, and
`overlap_cutoff`.  If the orthogonalization residual `result.orthogonalization_error` is large, reduce
basis redundancy by decreasing `sigma_scale` or increasing `overlap_cutoff`.
