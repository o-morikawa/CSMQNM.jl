# CSMQNM.jl

`CSMQNM.jl` is a cleaned-up Julia package skeleton for computing black-hole quasinormal-mode candidates with the complex scaling method (CSM).  The implementation is intentionally narrowed to the basis that is useful in the current project:

- `ical = 2` only, for legacy compatibility.
- Polynomial x Gaussian basis only.
- `range = :real` or `range = :complex` switches between polynomial x real-range Gaussian and polynomial x complex-range Gaussian.
- Real-range Gaussian alone, sin/cos x real-range Gaussian, and CLD calculation are not included.
- Physical parameters are passed with a dictionary, and `ipot` is specified by a symbolic label rather than by a hard-to-remember integer.

## Installation

From the parent directory of this package:

```julia
using Pkg
Pkg.activate("CSMQNM")
Pkg.instantiate()
```

If `AutoTortoise.jl` is not available from your package environment, add it in the same way you used it in the original scripts.  It is required for the dS/(A)dS tortoise-coordinate inversion used by `:sds4_rw` and `:sds_dim_rw`.

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
    basis = Dict(:imax=>30, :nmax=>3, :r0=>0.1, :rmax=>60.0, :range=>:complex),
    csm = Dict(:the0=>42.0),
    integration = Dict(:xmin=>-60.0, :xmax=>80.0, :dx=>0.02),
    physics = Dict(:ipot=>:sds4_rw, :M=>1.0, :ell=>2, :s=>2, :lam=>0.08),
    output = Dict(:write_spectrum=>true),
)
```

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

The input is split into four groups.

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
result.hamiltonian         # orthogonalized Hamiltonian matrix
result.basis_size_before
result.basis_size_after
result.dropped_basis_vectors
```

Set `OutputConfig(write_spectrum=true)` to write a spectrum file, and `OutputConfig(write_potential=true)` to write the potential sampled on the integration grid.

## Notes on the current scope

This package skeleton keeps the single-channel QNM computation only.  The old CLD code, hand-picked CLD pole indices, stand-alone real-range Gaussian branch, stand-alone complex-range Gaussian branch, and trigonometric basis branch were deliberately removed.  The `:stringy_ads_ds` label remains as a placeholder because the original implementation is coupled-channel and should be refactored separately rather than forced into the single-channel API.
