"Gaussian-packet ECS backend with Takagi orthogonalization."
Base.@kwdef struct TakagiECSConfig
    theta_deg::Float64 = 60.0
    xmin::Float64 = -240.0
    xmax::Float64 = 240.0
    nbasis::Int = 240
    sigma_scale::Float64 = 1.5
    x0_left::Float64 = 30.0
    x0_right::Float64 = 60.0
    smoothing::Float64 = 5.0
    quadrature_panels::Int = 0
    quadrature_order::Int = 12
    overlap_cutoff::Float64 = 1e-10
    orthogonalization::Symbol = :auto   # :auto, :takagi, :hermitian, or :none
    envelope_dirichlet::Bool = true
    boundary::Symbol = :dirichlet
end

"Full Takagi-ECS run configuration. `physics` uses the same labels as `RunConfig.physics`."
Base.@kwdef struct TakagiECSRunConfig
    ecs::TakagiECSConfig = TakagiECSConfig()
    physics::Dict{Symbol,Any} = Dict(:ipot => :schwarzschild_rw, :M => 1.0, :ell => 2, :s => 2, :lam => 0.0, :dim => 4, :c => :tensor)
    output::OutputConfig = OutputConfig()
end

struct TakagiECSResult
    energy::Vector{ComplexF64}
    omega::Vector{ComplexF64}
    hamiltonian::Matrix{ComplexF64}
    overlap::Matrix{ComplexF64}
    overlap_values::Vector{Float64}
    orthogonalization_error::Float64
    orthogonalization_method::Symbol
    basis_size_before::Int
    basis_size_after::Int
    dropped_basis_vectors::Int
    centers::Vector{Float64}
    sigmas::Vector{Float64}
    config::TakagiECSRunConfig
end

function solve_qnm_ecs_takagi_from_dict(; ecs=Dict{Symbol,Any}(), physics=Dict{Symbol,Any}(), output=Dict{Symbol,Any}())
    cfg = TakagiECSRunConfig(
        ecs = TakagiECSConfig(; symbol_keys(ecs)...),
        physics = symbol_keys(physics),
        output = OutputConfig(; symbol_keys(output)...),
    )
    return solve_qnm_ecs_takagi(cfg)
end

function validate_takagi_ecs_config(cfg::TakagiECSRunConfig)
    e = cfg.ecs
    e.boundary === :dirichlet || error("Only boundary = :dirichlet is currently supported in the Takagi-ECS solver.")
    e.xmax > e.xmin || error("ecs.xmax must be larger than ecs.xmin.")
    e.nbasis >= 3 || error("ecs.nbasis must be at least 3.")
    e.sigma_scale > 0 || error("ecs.sigma_scale must be positive.")
    e.x0_left > 0 || error("ecs.x0_left must be positive.")
    e.x0_right > 0 || error("ecs.x0_right must be positive.")
    e.smoothing >= 0 || error("ecs.smoothing must be non-negative.")
    e.quadrature_order >= 2 || error("ecs.quadrature_order must be at least 2.")
    e.quadrature_panels >= 0 || error("ecs.quadrature_panels must be non-negative. Use 0 for automatic panel count.")
    e.overlap_cutoff > 0 || error("ecs.overlap_cutoff must be positive.")
    e.orthogonalization in (:auto, :takagi, :hermitian, :none) || error("ecs.orthogonalization must be :auto, :takagi, :hermitian, or :none.")
    theta = e.theta_deg * DEG
    0 < theta < π || error("ecs.theta_deg must satisfy 0 < theta < 180 degrees.")
    e.xmin < -e.x0_left || error("ecs.xmin must be smaller than -ecs.x0_left, leaving a left ECS absorbing region.")
    e.xmax > e.x0_right || error("ecs.xmax must be larger than ecs.x0_right, leaving a right ECS absorbing region.")
    canonical_potential_label(cfg.physics)
    return nothing
end

function takagi_as_ecs_config(e::TakagiECSConfig)
    return ECSConfig(
        theta_deg = e.theta_deg,
        xmin = e.xmin,
        xmax = e.xmax,
        nx = max(5, e.nbasis + 2),
        x0_left = e.x0_left,
        x0_right = e.x0_right,
        smoothing = e.smoothing,
        boundary = e.boundary,
    )
end

function takagi_ecs_centers_sigmas(e::TakagiECSConfig)
    h = (e.xmax - e.xmin) / (e.nbasis + 1)
    centers = collect(range(e.xmin + h, e.xmax - h; length=e.nbasis))
    sigmas = fill(e.sigma_scale * h, e.nbasis)
    return centers, sigmas
end

function gauss_legendre_rule(n::Integer)
    n >= 1 || error("Gauss-Legendre order must be positive.")
    if n == 1
        return [0.0], [2.0]
    end
    beta = [i / sqrt(4.0 * i^2 - 1.0) for i in 1:n-1]
    vals, vecs = eigen(SymTridiagonal(zeros(n), beta))
    weights = 2.0 .* abs2.(vecs[1, :])
    return collect(vals), collect(weights)
end

function takagi_ecs_panel_count(e::TakagiECSConfig)
    e.quadrature_panels > 0 && return e.quadrature_panels
    return max(e.nbasis, 32)
end

function gaussian_packet_values(x::Real, centers::AbstractVector{<:Real}, sigmas::AbstractVector{<:Real}, e::TakagiECSConfig)
    n = length(centers)
    phi = Vector{Float64}(undef, n)
    dphi = Vector{Float64}(undef, n)
    L = e.xmax - e.xmin
    if e.envelope_dirichlet
        env = 4.0 * (x - e.xmin) * (e.xmax - x) / L^2
        denv = 4.0 * (e.xmin + e.xmax - 2.0 * x) / L^2
    else
        env = 1.0
        denv = 0.0
    end
    @inbounds for i in 1:n
        σ = sigmas[i]
        ξ = (x - centers[i]) / σ
        normc = 1.0 / (π^(0.25) * sqrt(σ))
        g = normc * exp(-0.5 * ξ^2)
        dg = -(x - centers[i]) / σ^2 * g
        phi[i] = env * g
        dphi[i] = denv * g + env * dg
    end
    return phi, dphi
end

function gaussian_ecs_matrices(cfg::TakagiECSRunConfig)
    e = cfg.ecs
    ecs = takagi_as_ecs_config(e)
    centers, sigmas = takagi_ecs_centers_sigmas(e)
    n = length(centers)
    N = zeros(ComplexF64, n, n)
    K = zeros(ComplexF64, n, n)
    Vmat = zeros(ComplexF64, n, n)

    panels = takagi_ecs_panel_count(e)
    qnodes, qweights = gauss_legendre_rule(e.quadrature_order)
    edges = collect(range(e.xmin, e.xmax; length=panels + 1))

    for p in 1:panels
        a = edges[p]
        b = edges[p + 1]
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)
        for q in eachindex(qnodes)
            x = mid + half * qnodes[q]
            wq = half * qweights[q]
            z = ecs_contour(x, ecs)
            J = ecs_jacobian(x, ecs)
            Vx = potential_value(cfg.physics, z)
            phi, dphi = gaussian_packet_values(x, centers, sigmas, e)
            phic = ComplexF64.(phi)
            dphic = ComplexF64.(dphi)
            N .+= (wq * J) .* (phic * transpose(phic))
            K .+= (wq / J) .* (dphic * transpose(dphic))
            Vmat .+= (wq * J * Vx) .* (phic * transpose(phic))
        end
    end
    H = K + Vmat
    return H, N, centers, sigmas
end

function factor_values_vector(d)
    if d isa Diagonal
        return Float64.(real.(d.diag))
    elseif d isa AbstractVector
        return Float64.(real.(d))
    else
        return Float64.(real.(diag(Matrix(d))))
    end
end

function try_takagi_orthogonalizer(N::AbstractMatrix; cutoff::Real=1e-10)
    d, U = takagi_factor(Matrix(N), sort=-1)
    vals = factor_values_vector(d)
    keep = findall(vals .> cutoff)
    isempty(keep) && error("All Takagi singular values of the overlap matrix were below cutoff=$cutoff.")

    # TakagiFactorization.jl convention used here: N ≈ transpose(U) * Diagonal(vals) * U.
    # Build S and verify S^T N S ≈ I. If the convention changes, this diagnostic catches it.
    S = adjoint(U)[:, keep] * Diagonal(1.0 ./ sqrt.(vals[keep]))
    err = norm(transpose(S) * Matrix(N) * S - I(length(keep))) / max(1.0, length(keep))
    if !isfinite(err) || err > 1e-6
        error("Takagi orthogonalization did not satisfy S^T N S ≈ I; normalized error=$err")
    end
    return Matrix{ComplexF64}(S), keep, vals, Float64(real(err)), :takagi
end

function hermitian_metric_matrix(cfg::TakagiECSRunConfig, centers, sigmas)
    e = cfg.ecs
    ecs = takagi_as_ecs_config(e)
    n = length(centers)
    G = zeros(ComplexF64, n, n)
    panels = takagi_ecs_panel_count(e)
    qnodes, qweights = gauss_legendre_rule(e.quadrature_order)
    edges = collect(range(e.xmin, e.xmax; length=panels + 1))

    for p in 1:panels
        a = edges[p]
        b = edges[p + 1]
        mid = 0.5 * (a + b)
        half = 0.5 * (b - a)
        for q in eachindex(qnodes)
            x = mid + half * qnodes[q]
            wq = half * qweights[q]
            J = ecs_jacobian(x, ecs)
            phi, _ = gaussian_packet_values(x, centers, sigmas, e)
            phic = ComplexF64.(phi)
            # Positive auxiliary metric for basis conditioning only.
            G .+= (wq * abs(J)) .* (phic * adjoint(phic))
        end
    end
    return Hermitian(G)
end

function hermitian_conditioning_transform(cfg::TakagiECSRunConfig, centers, sigmas; cutoff::Real=1e-10)
    G = hermitian_metric_matrix(cfg, centers, sigmas)
    F = eigen(G)
    vals = Float64.(real.(F.values))
    scale = maximum(abs.(vals))
    thresh = cutoff * max(scale, 1.0)
    keep = findall(vals .> thresh)
    isempty(keep) && error("All Hermitian auxiliary overlap eigenvalues were below cutoff=$cutoff.")
    S = F.vectors[:, keep] * Diagonal(1.0 ./ sqrt.(vals[keep]))
    err = norm(adjoint(S) * Matrix(G) * S - I(length(keep))) / max(1.0, length(keep))
    return Matrix{ComplexF64}(S), keep, vals, Float64(real(err)), :hermitian
end

function direct_generalized_conditioning(N::AbstractMatrix; cutoff::Real=1e-10)
    n = size(N, 1)
    vals = ones(Float64, n)
    keep = collect(1:n)
    return Matrix{ComplexF64}(I, n, n), keep, vals, 0.0, :none
end

function takagi_ecs_conditioner(cfg::TakagiECSRunConfig, N::AbstractMatrix, centers, sigmas)
    method = cfg.ecs.orthogonalization
    if method === :takagi
        return try_takagi_orthogonalizer(N; cutoff=cfg.ecs.overlap_cutoff)
    elseif method === :hermitian
        return hermitian_conditioning_transform(cfg, centers, sigmas; cutoff=cfg.ecs.overlap_cutoff)
    elseif method === :none
        return direct_generalized_conditioning(N; cutoff=cfg.ecs.overlap_cutoff)
    else
        try
            return try_takagi_orthogonalizer(N; cutoff=cfg.ecs.overlap_cutoff)
        catch err
            @warn "Takagi orthogonalization failed; falling back to Hermitian auxiliary-metric conditioning." exception=(err, catch_backtrace())
            return hermitian_conditioning_transform(cfg, centers, sigmas; cutoff=cfg.ecs.overlap_cutoff)
        end
    end
end

function solve_qnm_ecs_takagi(cfg::TakagiECSRunConfig)
    validate_takagi_ecs_config(cfg)
    H, N, centers, sigmas = gaussian_ecs_matrices(cfg)
    basis_size_before = size(H, 1)
    S, keep, vals, err, method = takagi_ecs_conditioner(cfg, N, centers, sigmas)
    if method === :takagi
        Hortho = transpose(S) * H * S
        ene = eigvals(Matrix(Hortho))
    elseif method === :none
        Hortho = H
        ene = eigvals(Matrix(H), Matrix(N))
    else
        # Hermitian auxiliary conditioning is only a basis reduction/preconditioner.
        # The physical C-product generalized problem is preserved in the reduced subspace.
        Hred = adjoint(S) * H * S
        Nred = adjoint(S) * N * S
        Hortho = Hred
        ene = eigvals(Matrix(Hred), Matrix(Nred))
    end
    ene = ComplexF64[w for w in ene if isfinite(real(w)) && isfinite(imag(w))]
    omega = select_qnm_branch.(sqrt.(ene))
    result = TakagiECSResult(
        ene,
        ComplexF64.(omega),
        Matrix{ComplexF64}(Hortho),
        Matrix{ComplexF64}(N),
        vals,
        err,
        method,
        basis_size_before,
        length(keep),
        basis_size_before - length(keep),
        centers,
        sigmas,
        cfg,
    )

    if cfg.output.write_potential
        path = something(cfg.output.potential_path, default_takagi_ecs_potential_filename(cfg))
        write_takagi_ecs_potential(result, path)
    end
    if cfg.output.write_spectrum
        path = something(cfg.output.spectrum_path, default_takagi_ecs_spectrum_filename(cfg))
        write_takagi_ecs_spectrum(result, path)
    end
    return result
end

function write_takagi_ecs_spectrum(result::TakagiECSResult, path::AbstractString)
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

function write_takagi_ecs_potential(result::TakagiECSResult, path::AbstractString)
    e = result.config.ecs
    ecs = takagi_as_ecs_config(e)
    nplot = max(501, 2 * e.nbasis + 1)
    xs = collect(range(e.xmin, e.xmax; length=nplot))
    open(path, "w") do io
        println(io, "# x  Re(z)  Im(z)  Re(V)  Im(V)  Re(J)  Im(J)")
        for x in xs
            z = ecs_contour(x, ecs)
            J = ecs_jacobian(x, ecs)
            V = potential_value(result.config.physics, z)
            @printf(io, "% .12e % .12e % .12e % .12e % .12e % .12e % .12e\n", x, real(z), imag(z), real(V), imag(V), real(J), imag(J))
        end
    end
    return path
end

function default_takagi_ecs_spectrum_filename(cfg::TakagiECSRunConfig)
    label = canonical_potential_label(cfg.physics)
    e = cfg.ecs
    return "spectrum_$(label)_takagi_ecs_theta$(e.theta_deg)_Nb$(e.nbasis)_x$(e.xmin)to$(e.xmax).dat"
end

function default_takagi_ecs_potential_filename(cfg::TakagiECSRunConfig)
    label = canonical_potential_label(cfg.physics)
    e = cfg.ecs
    return "potential_$(label)_takagi_ecs_theta$(e.theta_deg)_Nb$(e.nbasis).dat"
end
