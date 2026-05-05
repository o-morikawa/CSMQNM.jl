module CSMQNM

using LinearAlgebra
using Printf
using LambertW
using AutoTortoise

export BasisConfig, CSMConfig, IntegrationConfig, OutputConfig, RunConfig, QNMResult
export solve_qnm, solve_qnm_from_dict, write_potential, potential_value, potential_labels

const DEG = π / 180
const ORTHOGONALIZATION_CUTOFF = 1e-5

"Polynomial times Gaussian basis. `range = :real` or `:complex`."
Base.@kwdef struct BasisConfig
    imax::Int = 30
    nmax::Int = 3
    r0::Float64 = 0.1
    rmax::Float64 = 60.0
    range::Symbol = :real
    beta::Float64 = π / 2
end

"Complex scaling parameters. `the0` is in degrees. `ical` is kept fixed to 2 for legacy compatibility."
Base.@kwdef struct CSMConfig
    the0::Float64 = 40.0
    ical::Int = 2
end

"Grid used for numerical integration of potential matrix elements."
Base.@kwdef struct IntegrationConfig
    xmin::Float64 = -60.0
    xmax::Float64 = 80.0
    dx::Float64 = 0.02
end

"Output controls."
Base.@kwdef struct OutputConfig
    write_potential::Bool = false
    potential_path::Union{Nothing,String} = nothing
    write_spectrum::Bool = false
    spectrum_path::Union{Nothing,String} = nothing
end

"Full run configuration. `physics` accepts keys such as :ipot, :M, :ell, :s, :lam, :dim, :c."
Base.@kwdef struct RunConfig
    basis::BasisConfig = BasisConfig()
    csm::CSMConfig = CSMConfig()
    integration::IntegrationConfig = IntegrationConfig()
    physics::Dict{Symbol,Any} = Dict(:ipot => :sds_dim_rw, :M => 1.0, :ell => 2, :s => 0, :lam => 0.4, :dim => 5, :c => :tensor)
    output::OutputConfig = OutputConfig()
end

struct QNMResult
    energy::Vector{ComplexF64}
    omega::Vector{ComplexF64}
    hamiltonian::Matrix{ComplexF64}
    basis_size_before::Int
    basis_size_after::Int
    dropped_basis_vectors::Int
    config::RunConfig
end

potential_labels() = (
    :schwarzschild_rw,
    :ern_scalar,
    :ern_em,
    :ern_grav,
    :sds4_rw,
    :stringy_ads_ds,
    :sds_dim_rw,
)

function solve_qnm_from_dict(; basis=Dict{Symbol,Any}(), csm=Dict{Symbol,Any}(), integration=Dict{Symbol,Any}(), physics=Dict{Symbol,Any}(), output=Dict{Symbol,Any}())
    cfg = RunConfig(
        basis = BasisConfig(; symbol_keys(basis)...),
        csm = CSMConfig(; symbol_keys(csm)...),
        integration = IntegrationConfig(; symbol_keys(integration)...),
        physics = symbol_keys(physics),
        output = OutputConfig(; symbol_keys(output)...),
    )
    return solve_qnm(cfg)
end

function symbol_keys(d::AbstractDict)
    out = Dict{Symbol,Any}()
    for (k, v) in d
        out[Symbol(k)] = v
    end
    return out
end

function solve_qnm(cfg::RunConfig)
    validate_config(cfg)
    basis = cfg.basis
    csm = cfg.csm
    integ = cfg.integration
    theta = csm.the0 * DEG

    cn, alpha = normalization_polynomial(basis)
    nmat = norm_matrix_polynomial(basis, cn, alpha)
    kmat = kinetic_matrix_polynomial(basis, cn, alpha, theta)
    vmat = potential_matrix_polynomial(cfg, cn, alpha, theta)
    h0 = kmat + vmat

    h, kept, dropped = orthogonalize_basis(nmat, h0)
    ene = eigvals(h)
    omega = select_qnm_branch.(sqrt.(ene))

    result = QNMResult(ComplexF64.(ene), ComplexF64.(omega), h, size(h0, 1), kept, dropped, cfg)

    if cfg.output.write_potential
        path = something(cfg.output.potential_path, default_potential_filename(cfg))
        write_potential(cfg, path)
    end
    if cfg.output.write_spectrum
        path = something(cfg.output.spectrum_path, default_spectrum_filename(cfg))
        write_spectrum(result, path)
    end
    return result
end

function validate_config(cfg::RunConfig)
    cfg.csm.ical == 2 || error("Only ical = 2 is supported: Polynomial x Gaussian basis.")
    cfg.basis.range in (:real, :complex) || error("basis.range must be :real or :complex.")
    cfg.basis.imax >= 2 || error("basis.imax must be at least 2.")
    cfg.basis.nmax >= 0 || error("basis.nmax must be non-negative.")
    cfg.basis.r0 > 0 || error("basis.r0 must be positive.")
    cfg.basis.rmax > cfg.basis.r0 || error("basis.rmax must be larger than basis.r0.")
    cfg.integration.dx > 0 || error("integration.dx must be positive.")
    cfg.integration.xmax > cfg.integration.xmin || error("integration.xmax must be larger than integration.xmin.")
    canonical_potential_label(cfg.physics)
    return nothing
end

function canonical_potential_label(physics::AbstractDict)
    raw = get(physics, :ipot, get(physics, "ipot", :sds_dim_rw))
    if raw isa Integer
        table = Dict(1=>:schwarzschild_rw, 2=>:ern_scalar, 3=>:ern_em, 4=>:ern_grav, 5=>:sds4_rw, 6=>:stringy_ads_ds, 7=>:sds_dim_rw)
        haskey(table, raw) || error("Unknown legacy ipot=$raw")
        return table[raw]
    end
    label = Symbol(raw)
    label in potential_labels() || error("Unknown potential label: $label. Choose one of $(potential_labels()).")
    return label
end

function physics_get(physics::AbstractDict, key::Symbol, default)
    if haskey(physics, key)
        return physics[key]
    elseif haskey(physics, String(key))
        return physics[String(key)]
    else
        return default
    end
end

channel_c(c) = c === :vector ? 1 : c === :tensor ? 2 : Int(c)

function normalization_polynomial(basis::BasisConfig)
    imax, nmax = basis.imax, basis.nmax
    nrange = basis.range === :real ? 1 : 2
    base_ranges = [basis.r0 * (basis.rmax / basis.r0)^((i - 1) / (imax - 1)) for i in 1:imax]
    alpha_base = basis.range === :real ? ComplexF64.(1.0 ./ base_ranges.^2) :
        ComplexF64.(vcat((1.0 + im * basis.beta) ./ base_ranges.^2, (1.0 - im * basis.beta) ./ base_ranges.^2))

    nb = imax * nrange * (nmax + 1)
    alpha = Vector{ComplexF64}(undef, nb)
    cn = Vector{ComplexF64}(undef, nb)

    for n in 0:nmax
        for k in 1:(imax * nrange)
            idx = basis_index(k, n, imax * nrange)
            alpha[idx] = alpha_base[k]
            cn[idx] = sqrt(2.0^(3n + 0.5) * factorial(big(n)) * sqrt(alpha[idx]) / (factorial(big(2n)) * sqrt(π)))
        end
    end
    return cn, alpha
end

basis_index(k::Integer, n::Integer, ibase::Integer) = n * ibase + k
basis_dimension(basis::BasisConfig) = basis.imax * (basis.range === :real ? 1 : 2) * (basis.nmax + 1)

function gaussian_moment(power::Integer, a)
    isodd(power) && return zero(a)
    m = div(power, 2)
    return factorial(big(2m)) * sqrt(π) / (factorial(big(m)) * 2.0^(2m) * a^(m + 0.5))
end

function norm_matrix_polynomial(basis::BasisConfig, cn, alpha)
    ibase = basis.imax * (basis.range === :real ? 1 : 2)
    dim = basis_dimension(basis)
    nmat = zeros(ComplexF64, dim, dim)
    for np in 0:basis.nmax, n in 0:np
        for i0 in 1:ibase, j0 in 1:(np == n ? i0 : ibase)
            i = basis_index(i0, np, ibase)
            j = basis_index(j0, n, ibase)
            val = cn[i] * cn[j] * alpha[i]^(np/2) * alpha[j]^(n/2) * gaussian_moment(np + n, alpha[i] + alpha[j])
            nmat[i, j] = val
            nmat[j, i] = val
        end
    end
    return nmat
end

function kinetic_matrix_polynomial(basis::BasisConfig, cn, alpha, theta)
    ibase = basis.imax * (basis.range === :real ? 1 : 2)
    dim = basis_dimension(basis)
    kmat = zeros(ComplexF64, dim, dim)
    for np in 0:basis.nmax, n in 0:np
        for i0 in 1:ibase, j0 in 1:(np == n ? i0 : ibase)
            i = basis_index(i0, np, ibase)
            j = basis_index(j0, n, ibase)
            a = alpha[i] + alpha[j]
            val = zero(ComplexF64)
            if n >= 2
                val += n * (n - 1) * gaussian_moment(np + n - 2, a)
            end
            val += -2.0 * alpha[j] * (2n + 1) * gaussian_moment(np + n, a)
            val += 4.0 * alpha[j]^2 * gaussian_moment(np + n + 2, a)
            val *= -exp(-2im * theta) * cn[i] * cn[j] * alpha[i]^(np/2) * alpha[j]^(n/2)
            kmat[i, j] = val
            kmat[j, i] = val
        end
    end
    return kmat
end

function potential_matrix_polynomial(cfg::RunConfig, cn, alpha, theta)
    basis, integ = cfg.basis, cfg.integration
    ibase = basis.imax * (basis.range === :real ? 1 : 2)
    dim = basis_dimension(basis)
    vmat = zeros(ComplexF64, dim, dim)
    xs = collect(integ.xmin:integ.dx:integ.xmax)
    pot = [potential_value(cfg.physics, x) for x in xs]

    for np in 0:basis.nmax, n in 0:np
        phase = exp(-im * (np + n + 1) * theta)
        for i0 in 1:ibase, j0 in 1:(np == n ? i0 : ibase)
            i = basis_index(i0, np, ibase)
            j = basis_index(j0, n, ibase)
            a = (alpha[i] + alpha[j]) * exp(-2im * theta)
            fsum = zero(ComplexF64)
            for k in eachindex(xs)
                x = xs[k]
                fsum += pot[k] * x^(np+n) * exp(-a * x^2) * integ.dx
            end
            val = phase * cn[i] * cn[j] * alpha[i]^(np/2) * alpha[j]^(n/2) * fsum
            vmat[i, j] = val
            vmat[j, i] = val
        end
    end
    return vmat
end

function orthogonalize_basis(nmat, h0; cutoff=ORTHOGONALIZATION_CUTOFF)
    vals, vecs = eigen(Matrix(nmat))
    keep = findall(abs.(vals) .> cutoff)
    isempty(keep) && error("All overlap eigenvalues were below cutoff=$cutoff.")
    u = vecs[:, keep]
    khalf_inv = Diagonal(1 ./ sqrt.(vals[keep]))
    h = khalf_inv * transpose(u) * h0 * u * khalf_inv
    return Matrix{ComplexF64}(h), length(keep), length(vals) - length(keep)
end

select_qnm_branch(w::Complex) = imag(w) <= 0 ? ComplexF64(w) : ComplexF64(-w)

function potential_value(physics::AbstractDict, x::Real)
    label = canonical_potential_label(physics)
    M = Float64(physics_get(physics, :M, 1.0))
    ell = Int(physics_get(physics, :ell, 2))
    s = Int(physics_get(physics, :s, 0))
    lam = Float64(physics_get(physics, :lam, 0.0))
    dim = Int(physics_get(physics, :dim, 4))
    c = channel_c(physics_get(physics, :c, :tensor))

    if label === :schwarzschild_rw
        r = schwarzschild_r_of_rstar(x, M)
        return rw_potential(ell, M, r)
    elseif label === :ern_scalar
        r = ern_r_of_rstar(x, M)
        return ern_scalar_potential(ell, M, r)
    elseif label === :ern_em
        r = ern_r_of_rstar(x, M)
        return ern_em_potential(ell, M, r)
    elseif label === :ern_grav
        r = ern_r_of_rstar(x, M)
        return ern_grav_potential(ell, M, r)
    elseif label === :sds4_rw
        r = adsds_r_of_rstar(x, M, lam; dim=4)
        return rw_adsds_potential(s, ell, lam, M, r)
    elseif label === :sds_dim_rw
        r = adsds_r_of_rstar(x, M, lam; dim=dim)
        return rw_adsds_dim_potential(c, dim, ell, lam, M, r)
    elseif label === :stringy_ads_ds
        error("The single-channel solver does not support :stringy_ads_ds yet. Use one of the scalar potential labels, or implement a coupled-channel backend.")
    end
end

schwarzschild_r_of_rstar(x, M) = 2.0 * M * (1.0 + lambertw(exp(x / (2.0 * M) - 1.0)))

function adsds_r_of_rstar(x, M, lam; dim::Int=4)
    red_lam = dim == 4 ? lam / 3.0 : 2.0 * lam / ((dim - 1) * (dim - 2))
    f = MetricFunction(-(dim - 3) => -2.0 * M, 0 => 1.0, 2 => -red_lam)
    h = horizons(f)
    region = static_region(h; prefer=:finite)
    tm = tortoise_map(f, region)
    return inverse_tortoise(tm, x)
end

rw_potential(ell, M, r) = (1.0 - 2.0 * M / r) * (ell * (ell + 1) / r^2 - 6.0 * M / r^3)

function rw_adsds_potential(s, ell, lam, M, r)
    sigma = 1 - s^2
    f = 1.0 - 2.0 * M / r - lam * r^2 / 3.0
    return f * (ell * (ell + 1) / r^2 + 2.0 * sigma * M / r^3 - (s - 1.0) * (s - 2.0) * lam / 3.0)
end

function rw_adsds_dim_potential(c, dim, ell, lam, M, r)
    red_lam = 2.0 * lam / ((dim - 1) * (dim - 2))
    f = 1.0 - red_lam * r^2 - 2.0 * M / r^(dim - 3)
    fp = -2.0 * red_lam * r + 2.0 * M * (dim - 3) / r^(dim - 2)
    if c == 1
        return f / r^2 * (ell * (ell + dim - 3) + (dim - 4) + dim * (dim - 2) * f / 4.0 - (dim - 2) * r * fp / 2.0)
    elseif c == 2
        return f / r^2 * (ell * (ell + dim - 3) + (dim - 2) * r * fp / 2.0 + (dim - 2) * (dim - 4) * f / 4.0)
    else
        error("For higher-dimensional SdS choose c=:vector/1 or c=:tensor/2.")
    end
end

function ern_r_of_rstar(rstar, M; tol=1e-14, maxiter=100)
    y = rstar / M
    x = y == 0 ? 1.0 : y > 0 ? max(max(y, 1.0) - 2.0 * log(max(y, 1.0)), 1e-8) : max(-1.0 / y, 1e-8)
    for _ in 1:maxiter
        F = x - 1.0 / x + 2.0 * log(x) - y
        dF = (1.0 + 1.0 / x)^2
        xnew = x - F / dF
        if xnew <= 0 || !isfinite(xnew)
            xnew = x / 2.0
        end
        if abs(xnew - x) < tol * max(1.0, abs(x))
            return M * (1.0 + xnew)
        end
        x = xnew
    end
    error("Newton method did not converge while solving extremal RN tortoise map.")
end

ern_scalar_potential(ell, M, r) = (1.0 - M / r)^2 * (ell * (ell + 1) / r^2 + 2.0 * M * (r - M) / r^4)
ern_em_potential(ell, M, r) = (1.0 - M / r)^2 * (ell * (ell + 1) / r^2 - M * (3.0 - sqrt(4.0 * ell * (ell + 1.0) + 1.0)) / r^3 + 4.0 * M^2 / r^4)
ern_grav_potential(ell, M, r) = (1.0 - M / r)^2 * (ell * (ell + 1) / r^2 - M * (3.0 + sqrt(4.0 * ell * (ell + 1.0) + 1.0)) / r^3 + 4.0 * M^2 / r^4)

function write_potential(cfg::RunConfig, path::AbstractString)
    open(path, "w") do io
        println(io, "# x  V(x)")
        for x in cfg.integration.xmin:cfg.integration.dx:cfg.integration.xmax
            @printf(io, "% .12e % .12e\n", x, potential_value(cfg.physics, x))
        end
    end
    return path
end

function write_spectrum(result::QNMResult, path::AbstractString)
    open(path, "w") do io
        println(io, "# index  Re(E)  Im(E)  Re(omega)  Im(omega)")
        for i in eachindex(result.energy)
            E = result.energy[i]
            w = result.omega[i]
            @printf(io, "%5d % .12e % .12e % .12e % .12e\n", i, real(E), imag(E), real(w), imag(w))
        end
    end
    return path
end

function default_spectrum_filename(cfg::RunConfig)
    label = canonical_potential_label(cfg.physics)
    b = cfg.basis
    return "spectrum_$(label)_$(b.range)_poly_nmax$(b.nmax)_theta$(cfg.csm.the0)_N$(b.imax)_r$(b.r0)to$(b.rmax).dat"
end

function default_potential_filename(cfg::RunConfig)
    label = canonical_potential_label(cfg.physics)
    return "potential_$(label).dat"
end

end
