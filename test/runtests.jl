using Test
using CSMQNM

@testset "configuration" begin
    @test CSMQNM.canonical_potential_label(Dict(:ipot=>1)) == :schwarzschild_rw
    @test CSMQNM.canonical_potential_label(Dict(:ipot=>:sds_dim_rw)) == :sds_dim_rw
    @test CSMQNM.basis_dimension(BasisConfig(imax=4, nmax=2, range=:real)) == 12
    @test CSMQNM.basis_dimension(BasisConfig(imax=4, nmax=2, range=:complex)) == 12
end

@testset "complex range validation" begin
    valid = RunConfig(
        basis = BasisConfig(imax=4, nmax=1, range=:complex, beta=π/2),
        csm = CSMConfig(the0=42.0),
        physics = Dict(:ipot=>:schwarzschild_rw),
    )
    @test CSMQNM.validate_config(valid) === nothing

    bad_beta = RunConfig(
        basis = BasisConfig(imax=4, nmax=1, range=:complex, beta=0.1),
        csm = CSMConfig(the0=80.0),
        physics = Dict(:ipot=>:schwarzschild_rw),
    )
    @test_throws ErrorException CSMQNM.validate_config(bad_beta)
end

@testset "Schwarzschild potential map" begin
    physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2)
    @test isfinite(CSMQNM.potential_value(physics, 0.0))
end

@testset "full small Schwarzschild solve" begin
    cfg = RunConfig(
        basis = BasisConfig(imax=4, nmax=1, r0=0.2, rmax=4.0, range=:real),
        csm = CSMConfig(the0=25.0),
        integration = IntegrationConfig(xmin=-5.0, xmax=5.0, dx=0.5),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm(cfg)
    @test length(result.energy) == result.basis_size_after
    @test all(imag.(result.omega) .<= 1e-12)
end

@testset "full small complex-range solve" begin
    cfg = RunConfig(
        basis = BasisConfig(imax=4, nmax=3, r0=0.2, rmax=4.0, range=:complex, beta=π/2),
        csm = CSMConfig(the0=42.0),
        integration = IntegrationConfig(xmin=-5.0, xmax=5.0, dx=0.5),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm(cfg)
    @test length(result.energy) == cfg.basis.imax * (cfg.basis.nmax + 1)
    @test result.dropped_basis_vectors == 0
    @test all(isfinite, real.(result.energy))
    @test all(isfinite, imag.(result.energy))
end

@testset "ECS configuration and contour" begin
    ecs = ECSConfig(theta_deg=60.0, xmin=-10.0, xmax=10.0, nx=101, x0_left=3.0, x0_right=3.0)
    @test CSMQNM.ecs_contour(0.0, ecs) == 0.0 + 0.0im
    @test CSMQNM.ecs_jacobian(0.0, ecs) == 1.0 + 0.0im
    @test CSMQNM.ecs_jacobian(5.0, ecs) ≈ exp(im * 60.0 * π / 180)
end

@testset "small Schwarzschild ECS solve" begin
    cfg = ECSRunConfig(
        ecs = ECSConfig(theta_deg=50.0, xmin=-8.0, xmax=8.0, nx=81, x0_left=3.0, x0_right=3.0),
        physics = Dict(:ipot=>:schwarzschild_rw, :M=>1.0, :ell=>2),
    )
    result = solve_qnm_ecs(cfg)
    @test length(result.energy) == cfg.ecs.nx - 2
    @test length(result.x) == cfg.ecs.nx - 2
    @test all(isfinite, real.(result.energy))
    @test all(isfinite, imag.(result.energy))
end
