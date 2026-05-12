"Finite-element DVR backend on an exterior-complex-scaled contour."
Base.@kwdef struct FEDVRECSConfig
    theta_deg::Float64 = 60.0
    xmin::Float64 = -240.0
    xmax::Float64 = 240.0
    nelements::Int = 120
    order::Int = 10
    x0_left::Float64 = 25.0
    x0_right::Float64 = 25.0
    smoothing::Float64 = 2.5
    boundary::Symbol = :dirichlet
    element_breaks::Union{Nothing,Vector{Float64}} = nothing
end

"Full FEDVR-ECS run configuration. `physics` uses the same labels as `RunConfig.physics`."
Base.@kwdef struct FEDVRECSRunConfig
    ecs::FEDVRECSConfig = FEDVRECSConfig()
    physics::Dict{Symbol,Any} = Dict(:ipot => :schwarzschild_rw, :M => 1.0, :ell => 2, :s => 2, :lam => 0.0, :dim => 4, :c => :tensor)
    output::OutputConfig = OutputConfig()
end

struct FEDVRECSResult
    energy::Vector{ComplexF64}
    omega::Vector{ComplexF64}
    hamiltonian::Matrix{ComplexF64}
    overlap::Matrix{ComplexF64}
    x::Vector{Float64}
    z::Vector{ComplexF64}
    jacobian::Vector{ComplexF64}
    potential::Vector{ComplexF64}
    config::FEDVRECSRunConfig
end

function solve_qnm_ecs_fedvr_from_dict(; ecs=Dict{Symbol,Any}(), physics=Dict{Symbol,Any}(), output=Dict{Symbol,Any}())
    cfg = FEDVRECSRunConfig(
        ecs = FEDVRECSConfig(; symbol_keys(ecs)...),
        physics = symbol_keys(physics),
        output = OutputConfig(; symbol_keys(output)...),
    )
    return solve_qnm_ecs_fedvr(cfg)
end

function validate_fedvr_ecs_config(cfg::FEDVRECSRunConfig)
    e = cfg.ecs
    e.boundary === :dirichlet || error("Only boundary = :dirichlet is currently supported in the FEDVR-ECS solver.")
    e.xmax > e.xmin || error("ecs.xmax must be larger than ecs.xmin.")
    e.nelements >= 1 || error("ecs.nelements must be positive.")
    e.order >= 4 || error("ecs.order must be at least 4 for Gauss-Lobatto-Legendre FEDVR.")
    e.x0_left > 0 || error("ecs.x0_left must be positive.")
    e.x0_right > 0 || error("ecs.x0_right must be positive.")
    e.smoothing >= 0 || error("ecs.smoothing must be non-negative.")
    theta = e.theta_deg * DEG
    0 < theta < π || error("ecs.theta_deg must satisfy 0 < theta < 180 degrees.")
    e.xmin < -e.x0_left || error("ecs.xmin must be smaller than -ecs.x0_left, leaving a left ECS absorbing region.")
    e.xmax > e.x0_right || error("ecs.xmax must be larger than ecs.x0_right, leaving a right ECS absorbing region.")
    if e.element_breaks !== nothing
        b = e.element_breaks
        length(b) >= 2 || error("ecs.element_breaks must contain at least two points.")
        abs(b[1] - e.xmin) <= 100eps(Float64) * max(1.0, abs(e.xmin)) || error("ecs.element_breaks[1] must equal ecs.xmin.")
        abs(b[end] - e.xmax) <= 100eps(Float64) * max(1.0, abs(e.xmax)) || error("ecs.element_breaks[end] must equal ecs.xmax.")
        all(b[i + 1] > b[i] for i in 1:length(b)-1) || error("ecs.element_breaks must be strictly increasing.")
    end
    canonical_potential_label(cfg.physics)
    return nothing
end

"Gauss-Lobatto-Legendre nodes and weights on [-1,1]."
function legendre_gll(order::Integer)
    order >= 2 || error("GLL order must be at least 2.")
    n = order - 1
    ξ = zeros(Float64, order)
    w = zeros(Float64, order)
    ξ[1] = -1.0
    ξ[end] = 1.0
    if order > 2
        for i in 2:order-1
            # Chebyshev-Lobatto points are good initial guesses for the roots of P'_n.
            x = -cos((i - 1) * π / n)
            for _ in 1:80
                P, dP = legendre_p_and_dp(n, x)
                ddP = (2.0 * x * dP - n * (n + 1) * P) / (1.0 - x^2)
                step = dP / ddP
                xnew = x - step
                if abs(xnew - x) <= 10eps(Float64) * max(1.0, abs(x))
                    x = xnew
                    break
                end
                x = xnew
            end
            ξ[i] = x
        end
    end
    for i in 1:order
        P, _ = legendre_p_and_dp(n, ξ[i])
        w[i] = 2.0 / (n * (n + 1) * P^2)
    end
    return ξ, w
end

function legendre_p_and_dp(n::Integer, x::Real)
    if n == 0
        return 1.0, 0.0
    elseif n == 1
        return float(x), 1.0
    end
    pnm2 = 1.0
    pnm1 = float(x)
    pn = pnm1
    for k in 2:n
        pn = ((2k - 1) * x * pnm1 - (k - 1) * pnm2) / k
        pnm2, pnm1 = pnm1, pn
    end
    # Stable endpoint values.
    if abs(1.0 - x) < 100eps(Float64)
        dp = n * (n + 1) / 2.0
    elseif abs(-1.0 - x) < 100eps(Float64)
        dp = (-1)^(n + 1) * n * (n + 1) / 2.0
    else
        dp = n * (x * pn - pnm2) / (x^2 - 1.0)
    end
    return pn, dp
end

"Differentiation matrix for Lagrange polynomials at arbitrary distinct nodes."
function lagrange_derivative_matrix(nodes::AbstractVector{<:Real})
    n = length(nodes)
    λ = ones(Float64, n)
    for j in 1:n
        prod = 1.0
        xj = nodes[j]
        for k in 1:n
            k == j && continue
            prod *= xj - nodes[k]
        end
        λ[j] = 1.0 / prod
    end
    D = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            i == j && continue
            D[i, j] = λ[j] / (λ[i] * (nodes[i] - nodes[j]))
        end
        D[i, i] = -sum(D[i, j] for j in 1:n if j != i)
    end
    return D
end

function fedvr_element_breaks(e::FEDVRECSConfig)
    if e.element_breaks === nothing
        return collect(range(e.xmin, e.xmax; length=e.nelements + 1))
    else
        return copy(e.element_breaks)
    end
end

function fedvr_global_grid(e::FEDVRECSConfig)
    breaks = fedvr_element_breaks(e)
    ξ, wξ = legendre_gll(e.order)
    Dξ = lagrange_derivative_matrix(ξ)
    nloc = e.order
    ne = length(breaks) - 1
    nglobal = ne * (nloc - 1) + 1
    x = Vector{Float64}(undef, nglobal)
    for el in 1:ne
        a, b = breaks[el], breaks[el + 1]
        h = b - a
        offset = (el - 1) * (nloc - 1)
        for l in 1:nloc
            g = offset + l
            x[g] = 0.5 * (a + b) + 0.5 * h * ξ[l]
        end
    end
    return x, breaks, ξ, wξ, Dξ
end

function fedvr_ecs_matrices(cfg::FEDVRECSRunConfig)
    validate_fedvr_ecs_config(cfg)
    e = cfg.ecs
    x, breaks, ξ, wξ, Dξ = fedvr_global_grid(e)
    nloc = e.order
    ne = length(breaks) - 1
    nglobal = length(x)

    H = zeros(ComplexF64, nglobal, nglobal)
    N = zeros(ComplexF64, nglobal, nglobal)

    contour_cfg = ECSConfig(
        theta_deg=e.theta_deg,
        xmin=e.xmin,
        xmax=e.xmax,
        nx=5,
        x0_left=e.x0_left,
        x0_right=e.x0_right,
        smoothing=e.smoothing,
        boundary=e.boundary,
    )
    z = ComplexF64[ecs_contour(xi, contour_cfg) for xi in x]
    J = ComplexF64[ecs_jacobian(xi, contour_cfg) for xi in x]
    V = ComplexF64[potential_value(cfg.physics, zi) for zi in z]

    for el in 1:ne
        a, b = breaks[el], breaks[el + 1]
        h = b - a
        scale_dx = 2.0 / h
        jac_phys = h / 2.0
        offset = (el - 1) * (nloc - 1)
        inds = offset .+ (1:nloc)

        # Local weak form with GLL quadrature:
        # K_ab = ∫ dx J^{-1} dB_a/dx dB_b/dx,
        # N_ab = ∫ dx J B_a B_b,
        # V_ab = ∫ dx J V B_a B_b.
        for aidx in 1:nloc
            ga = inds[aidx]
            for bidx in 1:nloc
                gb = inds[bidx]
                kval = 0.0 + 0.0im
                for q in 1:nloc
                    gq = inds[q]
                    dBa = scale_dx * Dξ[q, aidx]
                    dBb = scale_dx * Dξ[q, bidx]
                    kval += jac_phys * wξ[q] * (1.0 / J[gq]) * dBa * dBb
                end
                H[ga, gb] += kval
            end
            # DVR diagonal mass and potential terms.
            g = ga
            weight = jac_phys * wξ[aidx]
            N[g, g] += weight * J[g]
            H[g, g] += weight * J[g] * V[g]
        end
    end

    # Dirichlet endpoints: remove first and last global nodes.
    interior = 2:nglobal-1
    return H[interior, interior], N[interior, interior], x[interior], z[interior], J[interior], V[interior]
end

function solve_qnm_ecs_fedvr(cfg::FEDVRECSRunConfig)
    H, N, x, z, J, V = fedvr_ecs_matrices(cfg)
    ene = eigvals(Matrix(H), Matrix(N))
    ene = ComplexF64[w for w in ene if isfinite(real(w)) && isfinite(imag(w))]
    omega = select_qnm_branch.(sqrt.(ene))
    result = FEDVRECSResult(
        ene,
        ComplexF64.(omega),
        Matrix{ComplexF64}(H),
        Matrix{ComplexF64}(N),
        x,
        z,
        J,
        V,
        cfg,
    )

    if cfg.output.write_potential
        path = something(cfg.output.potential_path, default_fedvr_ecs_potential_filename(cfg))
        write_fedvr_ecs_potential(result, path)
    end
    if cfg.output.write_spectrum
        path = something(cfg.output.spectrum_path, default_fedvr_ecs_spectrum_filename(cfg))
        write_fedvr_ecs_spectrum(result, path)
    end
    return result
end

function write_fedvr_ecs_spectrum(result::FEDVRECSResult, path::AbstractString)
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

function write_fedvr_ecs_potential(result::FEDVRECSResult, path::AbstractString)
    open(path, "w") do io
        println(io, "# Re(z)  Im(z)  Re(V)  Im(V)  Re(J)  Im(J)")
        for i in eachindex(result.z)
            z = result.z[i]
            V = result.potential[i]
            J = result.jacobian[i]
            @printf(io, "% .12e % .12e % .12e % .12e % .12e % .12e\n", real(z), imag(z), real(V), imag(V), real(J), imag(J))
        end
    end
    return path
end

function default_fedvr_ecs_spectrum_filename(cfg::FEDVRECSRunConfig)
    label = canonical_potential_label(cfg.physics)
    e = cfg.ecs
    return "spectrum_$(label)_fedvr_ecs_theta$(e.theta_deg)_Nel$(e.nelements)_p$(e.order)_x$(e.xmin)to$(e.xmax).dat"
end

function default_fedvr_ecs_potential_filename(cfg::FEDVRECSRunConfig)
    label = canonical_potential_label(cfg.physics)
    e = cfg.ecs
    return "potential_$(label)_fedvr_ecs_theta$(e.theta_deg)_Nel$(e.nelements)_p$(e.order).dat"
end
