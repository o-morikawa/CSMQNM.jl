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
    output = OutputConfig(
        write_spectrum = true,
        write_potential = true,
    ),
)

result = solve_qnm_ecs(cfg)
println(result.omega)
