using CSMQNM

# Gaussian-packet ECS backend with Takagi orthogonalization.
# This is an experimental alternative to finite-difference ECS and FEDVR-ECS.

ell = 3
θ = 60.0

cfg = TakagiECSRunConfig(
    ecs = TakagiECSConfig(
        theta_deg = θ,
        xmin = -240.0,
        xmax = 240.0,
        nbasis = 240,
        sigma_scale = 1.5,
        x0_left = 30.0,
        x0_right = 60.0,
        smoothing = 5.0,
        quadrature_order = 12,
        overlap_cutoff = 1e-10,
        envelope_dirichlet = true,
    ),
    physics = Dict(
        :ipot => :schwarzschild_rw,
        :M => 0.5,
        :ell => ell,
        :s => 2,
    ),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm_ecs_takagi(cfg)
println("kept basis vectors: ", result.basis_size_after, " / ", result.basis_size_before)
println("Takagi orthogonalization residual: ", result.takagi_error)
println("first few omega values:")
for w in result.omega[1:min(end, 10)]
    println(w)
end
