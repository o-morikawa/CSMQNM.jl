"Exterior complex scaling contour and finite-difference grid."
Base.@kwdef struct ECSConfig
    theta_deg::Float64 = 60.0
    xmin::Float64 = -120.0
    xmax::Float64 = 120.0
    nx::Int = 1201
    x0_left::Float64 = 40.0
    x0_right::Float64 = 40.0
    smoothing::Float64 = 0.0
    boundary::Symbol = :dirichlet
end

"Full ECS run configuration. `physics` uses the same labels as `RunConfig.physics`."
Base.@kwdef struct ECSRunConfig
    ecs::ECSConfig = ECSConfig()
    physics::Dict{Symbol,Any} = Dict(:ipot => :schwarzschild_rw, :M => 1.0, :ell => 2, :s => 2, :lam => 0.0, :dim => 4, :c => :tensor)
    output::OutputConfig = OutputConfig()
end

struct ECSResult
    energy::Vector{ComplexF64}
    omega::Vector{ComplexF64}
    hamiltonian::Matrix{ComplexF64}
    x::Vector{Float64}
    z::Vector{ComplexF64}
    jacobian::Vector{ComplexF64}
    potential::Vector{ComplexF64}
    config::ECSRunConfig
end

function solve_qnm_ecs_from_dict(; ecs=Dict{Symbol,Any}(), physics=Dict{Symbol,Any}(), output=Dict{Symbol,Any}())
    cfg = ECSRunConfig(
        ecs = ECSConfig(; symbol_keys(ecs)...),
        physics = symbol_keys(physics),
        output = OutputConfig(; symbol_keys(output)...),
    )
    return solve_qnm_ecs(cfg)
end

function validate_ecs_config(cfg::ECSRunConfig)
    e = cfg.ecs
    e.boundary === :dirichlet || error("Only boundary = :dirichlet is currently supported in the ECS finite-difference solver.")
    e.nx >= 5 || error("ecs.nx must be at least 5.")
    isodd(e.nx) || error("ecs.nx should be odd so that x=0 is included for symmetric domains.")
    e.xmax > e.xmin || error("ecs.xmax must be larger than ecs.xmin.")
    e.x0_left > 0 || error("ecs.x0_left must be positive.")
    e.x0_right > 0 || error("ecs.x0_right must be positive.")
    e.smoothing >= 0 || error("ecs.smoothing must be non-negative.")
    theta = e.theta_deg * DEG
    0 < theta < π || error("ecs.theta_deg must satisfy 0 < theta < 180 degrees.")
    e.xmin < -e.x0_left || error("ecs.xmin must be smaller than -ecs.x0_left, leaving a left ECS absorbing region.")
    e.xmax > e.x0_right || error("ecs.xmax must be larger than ecs.x0_right, leaving a right ECS absorbing region.")
    canonical_potential_label(cfg.physics)
    return nothing
end

"Return the ECS contour value z = g(x)."
function ecs_contour(x::Real, ecs::ECSConfig)
    θ = ecs.theta_deg * DEG
    q = exp(im * θ)
    if ecs.smoothing == 0
        if x < -ecs.x0_left
            return -ecs.x0_left + (x + ecs.x0_left) * q
        elseif x > ecs.x0_right
            return ecs.x0_right + (x - ecs.x0_right) * q
        else
            return complex(x)
        end
    else
        a = ecs.smoothing
        right = smooth_ramp(x - ecs.x0_right, a)
        left = smooth_ramp(-x - ecs.x0_left, a)
        return complex(x) + (q - 1) * right - (q - 1) * left
    end
end

"Return the ECS Jacobian J = dg/dx."
function ecs_jacobian(x::Real, ecs::ECSConfig)
    θ = ecs.theta_deg * DEG
    q = exp(im * θ)
    if ecs.smoothing == 0
        if x < -ecs.x0_left || x > ecs.x0_right
            return q
        else
            return 1.0 + 0im
        end
    else
        a = ecs.smoothing
        right = smooth_step(x - ecs.x0_right, a)
        left = smooth_step(-x - ecs.x0_left, a)
        return 1.0 + (q - 1) * right + (q - 1) * left
    end
end

smooth_step(y, a) = 0.5 * (1.0 + tanh(y / a))
smooth_ramp(y, a) = 0.5 * y + 0.5 * a * log(2.0 * cosh(y / a))

function solve_qnm_ecs(cfg::ECSRunConfig)
    validate_ecs_config(cfg)
    e = cfg.ecs
    x = collect(range(e.xmin, e.xmax; length=e.nx))
    dx = x[2] - x[1]
    z = ComplexF64[ecs_contour(xj, e) for xj in x]
    J = ComplexF64[ecs_jacobian(xj, e) for xj in x]
    Vfull = ComplexF64[potential_value(cfg.physics, zj) for zj in z]

    H = ecs_kinetic_matrix(x, e) + Diagonal(Vfull[2:end-1])
    ene = eigvals(Matrix(H))
    ene = ComplexF64[w for w in ene if isfinite(real(w)) && isfinite(imag(w))]
    omega = select_qnm_branch.(sqrt.(ene))

    result = ECSResult(
        ene,
        ComplexF64.(omega),
        Matrix{ComplexF64}(H),
        x[2:end-1],
        z[2:end-1],
        J[2:end-1],
        Vfull[2:end-1],
        cfg,
    )

    if cfg.output.write_potential
        path = something(cfg.output.potential_path, default_ecs_potential_filename(cfg))
        write_ecs_potential(result, path)
    end
    if cfg.output.write_spectrum
        path = something(cfg.output.spectrum_path, default_ecs_spectrum_filename(cfg))
        write_ecs_spectrum(result, path)
    end
    return result
end

"""
    ecs_kinetic_matrix(x, ecs)

Build the ECS kinetic-energy matrix on the interior grid points.  Dirichlet
boundary conditions are imposed by dropping the two endpoint degrees of freedom.

This is a conservative three-point flux-form discretization of

    -J(x)^(-1) d/dx [ J(x)^(-1) d/dx ],

where J(x) = dg/dx.  The half-grid coefficients are evaluated at the
midpoints x[j+1/2], which is more robust than composing two first-derivative
matrices and avoids even/odd grid decoupling artifacts.
"""
function ecs_kinetic_matrix(x::AbstractVector{<:Real}, ecs::ECSConfig)
    nx = length(x)
    nx >= 3 || error("At least three grid points are required.")

    dx = x[2] - x[1]
    all(abs((x[j + 1] - x[j]) - dx) <= 100eps(Float64) * max(1.0, abs(dx)) for j in 1:nx-1) ||
        error("ecs_kinetic_matrix currently requires a uniform grid.")

    n = nx - 2
    T = zeros(ComplexF64, n, n)

    # Node values on the interior points.
    Jinv_node = ComplexF64[1.0 / ecs_jacobian(x[j], ecs) for j in 2:nx-1]

    # Half-grid flux coefficients A_{j+1/2} = 1/J(x_{j+1/2}).
    # There are nx-1 intervals in the full grid.
    Jinv_half = ComplexF64[
        1.0 / ecs_jacobian(0.5 * (x[j] + x[j + 1]), ecs) for j in 1:nx-1
    ]

    invdx2 = 1.0 / dx^2
    for row in 1:n
        full_j = row + 1
        Aminus = Jinv_half[full_j - 1]
        Aplus = Jinv_half[full_j]
        pref = Jinv_node[row] * invdx2

        T[row, row] = pref * (Aminus + Aplus)
        if row > 1
            T[row, row - 1] = -pref * Aminus
        end
        if row < n
            T[row, row + 1] = -pref * Aplus
        end
    end
    return T
end

# Backward-compatible low-level constructor.  This keeps old internal calls
# working, but the preferred ECS path is `ecs_kinetic_matrix(x, ecs)` because it
# evaluates the flux coefficients at half-grid points.
function ecs_kinetic_matrix(J::AbstractVector{<:Complex}, dx::Real)
    nx = length(J)
    n = nx - 2
    Jinv = 1.0 ./ J
    T = zeros(ComplexF64, n, n)
    for row in 1:n
        j = row + 1
        Aplus = 0.5 * (Jinv[j] + Jinv[j + 1])
        Aminus = 0.5 * (Jinv[j] + Jinv[j - 1])
        pref = Jinv[j] / dx^2
        T[row, row] = pref * (Aplus + Aminus)
        if row < n
            T[row, row + 1] = -pref * Aplus
        end
        if row > 1
            T[row, row - 1] = -pref * Aminus
        end
    end
    return T
end

function write_ecs_spectrum(result::ECSResult, path::AbstractString)
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

function write_ecs_potential(result::ECSResult, path::AbstractString)
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

function default_ecs_spectrum_filename(cfg::ECSRunConfig)
    label = canonical_potential_label(cfg.physics)
    e = cfg.ecs
    return "spectrum_$(label)_ecs_theta$(e.theta_deg)_N$(e.nx)_x$(e.xmin)to$(e.xmax).dat"
end

function default_ecs_potential_filename(cfg::ECSRunConfig)
    label = canonical_potential_label(cfg.physics)
    e = cfg.ecs
    return "potential_$(label)_ecs_theta$(e.theta_deg)_N$(e.nx).dat"
end
