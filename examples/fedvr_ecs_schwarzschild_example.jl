using CSMQNM

ell = 3
θ = 60.0

cfg = FEDVRECSRunConfig(
    ecs = FEDVRECSConfig(
        theta_deg = θ,
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
        :ell => ell,
        :s => 2,
    ),
    output = OutputConfig(write_spectrum=true, write_potential=true),
)

result = solve_qnm_ecs_fedvr(cfg)
println("Computed ", length(result.omega), " FEDVR-ECS eigenfrequencies.")
